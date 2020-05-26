export apply_gate!,
    tebd!, TEBD2,
    TEBD2Sweep, TEBD4

abstract type TEBDalg end

"""
    TEBD2 <: TEBDalg
TEBD using 2nd order Trotter decomposition of the time evolution
operator. Namely, if the Hamiltonian is given by ``\\sum_i h_{i,i+1}`` with
``h_{i,i+1}`` having support only on the sites ``i,i+1``, define
``H_{e} = \\sum_{i \\el \\mathrm{even}}, H_{o} = \\sum_{i \\el \\mathrm{odd}}``.
The time-evolution operator is then approximated as:

``U(\\Delta t) = \\exp(-i \\Delta t/2 H_{o}) \\exp(-i \\Delta t H_{e}) \\exp(-i \\Delta t/2 H_{o})``

In practice the ``\\Delta t/2 H_{o}`` terms from consecutive time steps are grouped together as long as
no measurement is performed. Hence the overhead with respect to 1st order Trotter decomposition is
minimal.

The Trotter error per time step is ``O(\\Delta t^2)``.
"""
struct TEBD2 <: TEBDalg
end

"""
    TEBD2Sweep <: TEBDalg
TEBD using a 2nd Trotter decomposition similar to [`TEBD2`](@ref), but
instead of an even-odd decomposition of the Hamiltonian the time evolution
operator is approximated as:

``U(\\Delta t) = e^{-i \\Delta t /2 h_{1,2} } e^{-i \\Delta t /2 h_{2,3}} ... e^{-i \\Delta t /2 h_{2,3}} e^{-i \\Delta t /2 h_{1,2} }``

That is the exponentials of the bond Hamiltonians are applied to each bond sweeping from left to right and then applied a second
time sweeping from right to left.

The Trotter error per time step is ``O(\\Delta t^2)``
"""
struct TEBD2Sweep <: TEBDalg
end


"""
    TEBD4 <: TEBDalg
TEBD using a 4th order Trotter decomposition.
The time-evolution operator is approximated as ``U(τ₁)U(τ₂)U(τ₃)U(τ₂)U(τ₁)``
with ``U(τᵢ) = exp(-iH_o τᵢ/2) exp(-i H_e τᵢ ) exp(-i H_o τᵢ/2)``,
and ``τ1=τ2=1/(4-4^(1/3)) dt, τ3=dt - 4 τ1``.

A reduction in the number of gate applications is achieved by grouping together
the ``exp(-iH_o τᵢ/2)`` terms from consecutive Us and time-steps (except at measurement times).
"""
struct TEBD4 <: TEBDalg
end

function time_evo_gates(dt, H::GateList, alg::TEBD2; ishermitian=false)
    Uhalf = exp.(-1im*dt/2 .* H[1:2:end]; ishermitian=ishermitian)
    Us = [exp.( -1im*dt .* H[2:2:end]; ishermitian=ishermitian),
          exp.(-1im*dt .* H[1:2:end]; ishermitian=ishermitian)]
    return [Uhalf ], Us, [Us[1], Uhalf]
end

function time_evo_gates(dt,H::GateList,alg::TEBD2Sweep; ishermitian=false)
    Us = [exp.(-1im*dt/2 .* H, ishermitian=ishermitian),
          exp.(-1im*dt/2 .* H; ishermitian=ishermitian)]
    return [], Us, []
end

function time_evo_gates(dt,H::GateList,alg::TEBD4; ishermitian=false)
    τ₁ = 1/(4-4^(1/3))*dt
    τ₂ = τ₁
    τ₃ = dt - 2*τ₁ - 2*τ₂

    e= 2:2:length(H)
    o = 1:2:length(H)
    Ustart = [exp.(-1im*τ₁/2 .* H[o])]

    sequence = [(τ₁, e), (τ₁, o ), (τ₂,e), ((τ₂+τ₃)/2, o), (τ₃, e),
                ((τ₂+τ₃)/2, o), (τ₂,e), (τ₂,o), (τ₁,e), (τ₁,o)]
    Us = map(x->exp.(-1im*x[1] .* H[x[2]]; ishermitian=ishermitian), sequence)

    end_sequence = [(τ₁, e), (τ₁, o ), (τ₂,e), ((τ₂+τ₃)/2, o), (τ₃, e),
                    ((τ₂+τ₃)/2, o), (τ₂,e), (τ₂,o), (τ₁,e), (τ₁/2,o)]
    Uend = map(x->exp.(-1im*x[1] .* H[x[2]]; ishermitian=ishermitian), end_sequence)

    return Ustart, Us , Uend
end

tebd!(psi::MPS,H::BondOperator, args...; kwargs...) = tebd!(psi,gates(H),args...;kwargs...)

function tebd!(psi::MPS, H::GateList, dt::Number, tf::Number, alg::TEBDalg = TEBD2() ; kwargs... )
    # TODO: think of the best way to avoid inexact error when dt is very small
    # one option would be to use round(tf/dt) and verify that abs(round(tf/dt)-tf/dt)
    # is smaller than some threshold. Another option would be to use big(Rational(dt)).
    nsteps = Int(tf/dt)

    # TODO: use ishermitian for imaginary time-evolution once exponentiation
    # of hermitian ITensor is fixed (see https://github.com/ITensor/NDTensors.jl/pull/15).
    # ishermitian = get(kwargs, :ishermitian, dt isa Complex && real(dt)==0)
    ishermitian = get(kwargs, :ishermitian, false)

    cb = get(kwargs,:callback, NoTEvoCallback())
    cb_func(dt,step,rev,se) = (psi; bond, spec) -> apply!(cb,psi;
                                                              t=dt*step,
                                                              bond=bond,
                                                              spec=spec,
                                                              sweepdir=rev ? "left" : "right",
                                                              sweepend= se,
                                                              alg = alg)
    orthogonalize_step = get(kwargs,:orthogonalize,0)

    #We can bunch together half-time steps, when we don't need to measure observables
    dtm = callback_dt(cb)
    if dtm > 0
        floor(Int64, dtm / dt) != dtm /dt && throw("Measurement time step $dtm incommensurate with time-evolution time step $dt")
        mstep = floor(Int64,dtm / dt)
        nbunch = gcd(mstep,nsteps)
    else
        nbunch = nsteps
    end
    orthogonalize_step > 0 && (nbunch = gcd(nbunch,orthogonalize_step))

    Ustart, Us, Uend = time_evo_gates(dt,H,alg; ishermitian=ishermitian)

    length(Ustart)==0 && length(Uend)==0 && (nbunch+=1)

    pbar = get(kwargs,:progress, true) ? Progress(nsteps, desc="Evolving state... ") : nothing

    step = 0
    switchdir(rev) = !(rev)
    rev = false
    while step < nsteps
        for U in Ustart
            apply_gates!(psi, U ; reversed = rev, kwargs...)
            rev = !rev
        end

        for i in 1:nbunch-1

            stime = @elapsed begin
            for U in Us
                apply_gates!(psi,U; reversed=rev,
                             cb=cb_func(dt,step,rev,false), kwargs...)
                rev = !rev
            end
            end
            step += 1

            !isnothing(pbar) && ProgressMeter.next!(pbar, showvalues=[("t", dt*step),
                                                                  ("dt step time", round(stime,digits=3)),
                                                                  ("Max bond-dim", maxlinkdim(psi))])
            checkdone!(cb,psi) && break
        end

        #finalize the last time step from the bunched steps
        stime = @elapsed begin
        for (i,U) in enumerate(Uend)
            apply_gates!(psi,U; reversed = rev, cb=cb_func(dt,step+1,rev,i>=length(Uend)-1), kwargs...)
            rev = !rev
        end
        end

        if length(Uend)>0
            step += 1
            if !isnothing(pbar)
                ProgressMeter.next!(pbar,
                                    showvalues=[("t", dt*step),
                                                ("dt step time", round(stime,digits=3)),
                                                ("Max bond-dim", maxlinkdim(psi))])
            end
        end

        # TODO: make this a callback
        (orthogonalize_step>0 && step % orthogonalize_step ==0) && reorthogonalize!(psi)

        checkdone!(cb,psi) && break
    end
    return psi
end



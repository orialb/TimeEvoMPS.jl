export apply_gate!,
       tebd!, TEBD2, TEBD1

abstract type TEBDalg end

struct TEBD1 <: TEBDalg
end

struct TEBD2 <: TEBDalg
end

function time_evo_gates(dt, H::GateList, alg::TEBD2)
    Uhalf = exp.(-1im*dt/2 .* H[1:2:end])
    Us = [exp.( -1im*dt .* H[2:2:end]), exp.(-1im*dt .* H[1:2:end])]
    return [Uhalf ], Us, [Us[1], Uhalf]
end

function time_evo_gates(dt,H::GateList,alg::TEBD1)
    Us = [exp.(-1im*dt/2 .* H), exp.(-1im*dt/2 .* H[end:-1:1])]
    return [], Us, []
end

tebd!(psi::MPS,H::BondOperator, args...; kwargs...) = tebd!(psi,gates(H),args...;kwargs...)

function tebd!(psi::MPS, H::GateList, dt::Number, tf::Number, alg::TEBDalg = TEBD2() ; kwargs... )
    # TODO: think of the best way to avoid inexact error when dt is very small
    # one option would be to use round(tf/dt) and verify that abs(round(tf/dt)-tf/dt)
    # is smaller than some threshold. Another option would be to use big(Rational(dt)).
    nsteps = Int(tf/dt)
    obs = get(kwargs,:observer, NoTEvoObserver())
    orthogonalize_step = get(kwargs,:orthogonalize,0)

    #We can bunch together half-time steps, when we don't need to measure observables
    nbunch =  measurement_step(obs) >0 ? gcd(measurement_step(obs),nsteps) : nsteps
    orthogonalize_step > 0 && (nbunch = gcd(nbunch,orthogonalize_step))

    Ustart, Us, Uend = time_evo_gates(dt,H,alg)

    step = 0
    switchdir(dir) =dir== "fromleft" ? "fromright" : "fromleft"
    dir = "fromleft"
    while step < nsteps
        for U in Ustart
            apply_gates!(psi, U ; dir = dir, kwargs...)
            dir = switchdir(dir)
        end

        for i in 1:nbunch-1
            for U in Us
                apply_gates!(psi,U; dir=dir, kwargs...)
                dir = switchdir(dir)
            end
            step += 1
            observe!(obs,psi)
            checkdone!(obs,psi) && break
        end

        #finalize the last time step from the bunched steps
        for U in Uend
            apply_gates!(psi,U; dir=dir, kwargs...)
            dir = switchdir(dir)
        end
        step += 1
        (orthogonalize_step>0 && step % orthogonalize_step ==0) && reorthogonalize!(psi)
        observe!(obs,psi)
        checkdone!(obs,psi) && break
    end
    return psi
end



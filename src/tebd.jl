export apply_gate!,
       tebd!

tebd!(psi::MPS,H::BondOperator, args...; kwargs...) = tebd!(psi,gates(H),args...;kwargs...)

function tebd!(psi::MPS, H::GateList, dt::Number, tf::Number ; kwargs... )
    # TODO: think of the best way to avoid inexact error when dt is very small
    # one option would be to use round(tf/dt) and verify that abs(round(tf/dt)-tf/dt)
    # is smaller than some threshold. Another option would be to use big(Rational(dt)).
    nsteps = Int(tf/dt)
    order = get(kwargs, :order,2)
    order == 2 || throw("Time evolution with trotter decomposition of order
                        $order is currently not supported")
    obs = get(kwargs,:observer, NoTEvoObserver())
    orthogonalize_step = get(kwargs,:orthogonalize,0)

    #We can bunch together half-time steps, when we don't need to measure observables
    nbunch =  measurement_step(obs) >0 ? gcd(measurement_step(obs),nsteps) : nsteps
    orthogonalize_step > 0 && (nbunch = gcd(nbunch,orthogonalize_step))

    L = length(psi)

    Uhalf = exp.(-1im*dt/2 .* H[1:2:L-1],)
    Us = [exp.( -1im*dt .* H[2:2:L-1]), exp.(-1im*dt .* H[1:2:L-1])]

    step = 0
    switchdir(dir) =dir== "fromleft" ? "fromright" : "fromleft"
    dir = "fromleft"
    while step < nsteps
        #TODO : optimize sweep direction!
        apply_gates!(psi, Uhalf ; dir = dir, kwargs...)
        dir = switchdir(dir)
        for i in 1:nbunch-1
            for (j,U) in enumerate(Us)
                apply_gates!(psi,U; dir=dir, kwargs...)
                dir = switchdir(dir)
            end
            step += 1
            observe!(obs,psi)
            checkdone!(obs,psi) && break
        end
        apply_gates!(psi,Us[1]; dir=dir, kwargs...)
        dir = switchdir(dir)
        apply_gates!(psi,Uhalf; dir = dir,kwargs...)
        dir = switchdir(dir)

        step += 1
        (orthogonalize_step>0 && step % orthogonalize_step ==0) && reorthogonalize!(psi)
        observe!(obs,psi)
        checkdone!(obs,psi) && break
    end
    return psi
end



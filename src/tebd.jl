struct BondEvoGates
end

function apply_bond_gates!(psi::MPS, U::BondEvoGates, which_bonds=:all; kwargs...)
end

function tdmrg!(psi::MPS, H::BondOperator, dt::Number, tf::Number ; kwargs... )
    nsteps = Int(tf/dt)
    order = get(kwargs, :order,2)
    obs = get(kwargs,:observer, TEvoObserver(dt,tf; kwargs...))
    orthogonalize_step = get(kwargs,:orthogonalize,0)

    #We can bunch together half-time steps, when we don't need to measure observables
    nbunch =  gcd(measurement_step(obs),nsteps)
    orthogonalize_step > 0 && nbunch = gcd(nbunch,orthogonalize_step)

    U = BondEvoGates(H)
    apply_bond_gates!(psi, U, dt/2, :odd; kwargs...)

    while currentstep(obs) < nsteps
        for i in 1:nbunch-1
            for which_bonds in [:even,:odd]
                apply_bond_gates!(psi,U,dt,which_bonds; kwargs...)
            end
            step!(obs)
            observe!(obs,psi)
            checkdone!(obs,psi) && break
        end
        apply_bond_gates!(psi,U,dt,:even; kwargs...)
        apply_bond_gates!(psi,U,dt/2,:odd; kwargs...)

        step!(obs)
        observe!(obs,psi)
        checkdone!(obs,psi) && break
    end
    return psi, obs
end

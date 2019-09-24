export BondGate,apply_gate!,BondOperator,
       tdmrg!

function apply_gate!(psi::MPS,G::BondGate ; kwargs...)
    b = bond(G)
    orthogonalize!(psi,b)
    wf = psi[b]*psi[b+1]
    wf = noprime( G*(psi[b]*psi[b+1]) )
    replaceBond!(psi, b, wf; kwargs...)
end

function apply_gates!(psi::MPS, Gs::Vector{BondGate} ; kwargs...)
    for G in Gs
        apply_gate!(psi,G ; kwargs...)
    end
end

function tdmrg!(psi::MPS, H::BondOperator, dt::Number, tf::Number ; kwargs... )
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

    Uhalf = exp.(H[1:2:L-1],-1im*dt/2)
    imag_evo = typeof(dt) <: Complex && real(dt) ==0
    if !(imag_evo)
        Us = [exp.(H[reverse(2:2:L-1)], -1im*dt), exp.(H[1:2:L-1], -1im*dt)]
    else
        Us = [exp.(H,-1im*dt)]
    end

    step = 0
    while step < nsteps
        apply_gates!(psi, Uhalf ; dir = "fromleft", kwargs...)
        for i in 1:nbunch-1
            for (j,U) in enumerate(Us)
                dir = (j==1 && !(imag_evo)) ? "fromright" : "fromleft"
                apply_gates!(psi,U; dir=dir, kwargs...)
            end
            step += 1
            observe!(obs,psi)
            checkdone!(obs,psi) && break
        end
        apply_gates!(psi,Us[1]; dir="fromright", kwargs...)
        apply_gates!(psi,Uhalf; dir = "fromleft",kwargs...)

        step += 1
        observe!(obs,psi)
        checkdone!(obs,psi) && break
    end
    return psi
end


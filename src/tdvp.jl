export tdvp!
using ITensors: position!

singlesite!(PH::ProjMPO) = (PH.nsite = 1)
twosite!(PH::ProjMPO) = (PH.nsite = 2)

function tdvp!(psi,H::MPO,dt,tf; kwargs...)
    nsteps = Int(tf/dt)
    obs = get(kwargs,:observer, NoTEvoObserver())
    orthogonalize_step = get(kwargs,:orthogonalize,0)
    hermitian = get(kwargs,:hermitian,false)

    N = length(psi)
    orthogonalize!(psi,1)
    PH = ProjMPO(H)
    position!(PH,psi,1)
    for s in 1:nsteps
        for (b,ha) in sweepnext(N)
            #evolve with two-site Hamiltonian
            twosite!(PH)
            ITensors.position!(PH,psi,b)
            wf = psi[b]*psi[b+1]
            wf, info = exponentiate(PH, -1im*dt/2, wf; ishermitian=hermitian )
            dir = ha==1 ? "fromleft" : "fromright"
            info.converged==0 && throw("exponentiate did not converge")
            replaceBond!(psi,b,wf; dir = dir, kwargs... )

            #evolve with single-site Hamiltonian backward in time
            singlesite!(PH)
            i = ha==1 ? b+1 : b
            ITensors.position!(PH,psi,i)
            psi[i], info = exponentiate(PH,1im*dt/2,psi[i]; ishermitian=hermitian)
            info.converged==0 && throw("exponentiate did not converge")
        end
        checkdone!(obs) && break
    end

end

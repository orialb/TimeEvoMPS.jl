export tdvp!
using ITensors: position!

singlesite!(PH::ProjMPO) = (PH.nsite = 1)
twosite!(PH::ProjMPO) = (PH.nsite = 2)

function tdvp!(psi,H::MPO,dt,tf; kwargs...)
    nsteps = Int(tf/dt)
    obs = get(kwargs,:observer, NoTEvoObserver())
    orthogonalize_step = get(kwargs,:orthogonalize,0)
    hermitian = get(kwargs,:hermitian,false)
    exp_tol = get(kwargs,:exp_tol, 1e-12/dt)
    verbose = get(kwargs,:verbose, false)

    N = length(psi)
    orthogonalize!(psi,1)
    PH = ProjMPO(H)
    position!(PH,psi,1)
    for s in 1:nsteps
        stime = @elapsed begin
        for (b,ha) in sweepnext(N)
            #evolve with two-site Hamiltonian
            twosite!(PH)
            ITensors.position!(PH,psi,b)
            wf = psi[b]*psi[b+1]
            # @assert scalar(dag(wf)*wf) ≈ 1
            wf, info = exponentiate(PH, -1im*dt/2, wf; ishermitian=hermitian , tol=exp_tol)
            dir = ha==1 ? "fromleft" : "fromright"
            info.converged==0 && throw("exponentiate did not converge")
            replaceBond!(psi,b,wf; dir = dir, kwargs... )

            #evolve with single-site Hamiltonian backward in time
            singlesite!(PH)
            i = ha==1 ? b+1 : b
            if 1<i<N
                ITensors.position!(PH,psi,i)
                psi[i], info = exponentiate(PH,1im*dt/2,psi[i]; ishermitian=hermitian, tol=exp_tol)
                info.converged==0 && throw("exponentiate did not converge")
            end
            @assert b==1 || isleftortho(psi,i-1)
            @assert i==N || isrightortho(psi,i+1)
            # @assert scalar(dag(psi[i])*psi[i]) ≈ 1
        end
        (orthogonalize_step>0 && s % orthogonalize_step ==0) && reorthogonalize!(psi)
        end
        if verbose
            @printf("Step %d : maxLinkDim= %d, step time= %.3f \n", s,maxLinkDim(psi), stime)
        end
        checkdone!(obs) && break
        # println("step $s, $(inner(psi,psi))")
    end

end

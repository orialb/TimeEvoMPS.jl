using LinearAlgebra

function init_neel(L)
    N = floor(Int64,L/2)
    U0 = zeros(L,N)

    # init in Naeel state
    for i in 1:N
        U0[i*2-1,i]=1.
    end

    return U0
end

function nn_hopping_hamiltonian(L; periodic = false)
    h= diagm(1=> ones(L-1), -1=> ones(L-1))
    if periodic
        h[1,end] = 1
        h[end,1] = 1
    end
    return h
end

function free_fermions_densities(N,dt,tf)
    U = complex(init_neel(N))
    Unext = zeros(ComplexF64,size(U))

    expH = exp(-1im*dt*nn_hopping_hamiltonian(N))
    ns = [real.(diag(U*U'))]

    for i in 1:floor(Int64,tf/dt)
        #evolution
        Unext .=  expH*U
        # make sure U is normalized
        Q,R = qr(Unext)
        U .= Matrix(Q)
        push!(ns, real.(diag(U*U')))
    end
    return ns
end



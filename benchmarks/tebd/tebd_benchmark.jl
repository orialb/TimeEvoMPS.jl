using ITensors, TimeEvoMPS
using BenchmarkTools

function run_tebd_tfi(N,mindim,dt,tf,alg=TEBD2(); kwargs...)
    J = 1.
    h = 0.5
    sites = spinHalfSites(N)
    psi = productMPS(sites,ones(Int,N))
    H = tfi_bondop(sites,J,h)
    Hgates = gates(H)

    tebd!(psi,H,dt,tf,alg; mindim = mindim , kwargs...)
    println("max bond dim = $(maxLinkDim(psi))")
    return nothing
end

function run_tebd_xxx(N,mindim,dt,tf,alg=TEBD2(); kwargs...)
    sites = spinHalfSites(N)
    psi = productMPS(sites,ones(Int,N))
    H = BondOperator(sites)
    for b in 1:N-1
        for o in ["Sz","Sx","Sy"]
            add!(H,1.,o,o,b)
        end
    end
    Hgates = gates(H)
    tebd!(psi,H,dt,tf,alg; mindim=mindim, kwargs...)
    println("max bond dim = $(maxLinkDim(psi))")
    return nothing
end

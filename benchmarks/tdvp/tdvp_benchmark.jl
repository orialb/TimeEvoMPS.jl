using ITensors, TimeEvoMPS


function run_tdvp(N,maxdim,tf; kwargs...)
    J = 1.
    h = 0.5
    sites = spinHalfSites(N)
    H = tfi_mpo(J,h,sites)
    psi = complex!(productMPS(sites, ones(Int,N)))

    tdvp!(psi,H,0.1,tf; maxdim=maxdim, kwargs...)
    return maxLinkDim(psi)
end


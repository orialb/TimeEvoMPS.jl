using Test, ITensors, TimeEvoMPS

@testset "Imaginary time-evolution TFI model" begin
    N=10
    J = 1.
    h = 0.5
    sites = spinHalfSites(N)
    psi = productMPS(sites,ones(Int,N))
    #hack to make psi complex
    for b in 1:length(psi)
        psi[b] *= (1.0 + 0.0im)
    end

    ampo = AutoMPO()
    for j=1:length(sites)-1
        add!(ampo,-1.,"Sz",j,"Sz",j+1)
        add!(ampo,-h,"Sx",j)
    end
    add!(ampo,-h,"Sx",length(sites))
    H = MPO(ampo,sites)

    Es = []
    for dt in [0.1,0.01,1e-3,1e-4]
        tdvp!(psi,H,-1im*dt,-1im*500*dt ; maxdim=50, cutoff=1e-12, hermitian=true)
        push!(Es, inner(psi,H,psi))
    end

    # exact expression for ground-state energy
    # at criticality (ref? took it from ITensors.jl tests)

    eexact = 0.25 -0.25/sin(π/(4*N + 2))
    @test Es[end] ≈ eexact atol=1e-4
end


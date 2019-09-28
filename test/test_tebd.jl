using Test, TimeEvoMPS, ITensors
using TimeEvoMPS: isleftortho, isrightortho

@testset "Basic TEBD tests" begin
    # check that evolving with identity doesn't change the state
    N=10
    sites = spinHalfSites(N)
    trivialH = BondOperator(sites)
    psi0 = randomMPS(sites)
    orthogonalize!(psi0,1)
    psi = deepcopy(psi0)
    tebd!(psi,trivialH,0.01,1.;)

    @test inner(psi0,psi) ≈ 1

    #check that evolving forward and backward in time brings us back to the same state
    #(up to numerical errors of course)
    J,h = 0.3, -0.7
    H = tfi_bondop(sites,J,h)

    psi0 = randomMPS(sites)
    psi = deepcopy(psi0)
    tebd!(psi,H,0.01,1.)
    tebd!(psi,H,-0.01,-1.)

    @test inner(psi0,psi) ≈ 1

    #check that bond dimension is growing to maximum during evolution
    psi = productMPS(sites,ones(Int,N))
    tebd!(psi,H,0.01,5.)
    @test maxLinkDim(psi) == 2^5
end

@testset "Imaginary time-evolution TFI model" begin
    N=10
    J = 1.
    h = 0.5
    sites = spinHalfSites(N)
    psi = productMPS(sites,ones(Int,N))
    H = tfi_bondop(sites,J,h)
    Hgates = gates(H)

    Es = []
    for dt in [0.1,0.01,1e-3,1e-4]
        tebd!(psi,H,-1im*dt,-1im*500*dt ; maxdim=50,orthogonalize=10)
        push!(Es, measure(gates(H),psi))
    end

    # exact expression for ground-state energy
    # at criticality (ref? took it from ITensors.jl tests)

    eexact = 0.25 -0.25/sin(π/(4*N + 2))
    @test Es[end] ≈ eexact atol=1e-4
end

using Test, TimeEvoMPS, ITensors
using TimeEvoMPS: isleftortho, isrightortho

"create a BondOperator for the transverse-field ising model"
function tfi_bondop(sites,J,h)
    N = length(sites)
    H = BondOperator(sites)
    for b in 1:N-1
        add!(H,J,"Sz","Sz",b)
        add!(H,h,"Sx",b)
    end
    add!(H,h,"Sx",N)
    return H
end

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
    tebd!(psi,H,0.01,5.)
    tebd!(psi,H,-0.01,-5.)

    @test inner(psi0,psi) ≈ 1
    @test maxLinkDim(psi) == 2^5
end

"get ground state of transverse-field Ising model"
function TFIgs(sites,h)
    ampo = AutoMPO()
    for j=1:length(sites)-1
        add!(ampo,-1.,"Sz",j,"Sz",j+1)
        add!(ampo,-h,"Sx",j)
    end
    add!(ampo,-h,"Sx",length(sites))
    H = MPO(ampo,sites)

    psi0 = randomMPS(sites)
    sweeps = Sweeps(15)
    maxdim!(sweeps, 10,20,100,100,200)
    cutoff!(sweeps, 1E-10)
    energy, psi = dmrg(H,psi0,sweeps,quiet=true)
    return psi,energy
end

@testset "Imaginary time-evolution TFI model" begin
    N=10
    J = 1.
    h = 0.5
    sites = spinHalfSites(N)
    psi = randomMPS(sites)
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

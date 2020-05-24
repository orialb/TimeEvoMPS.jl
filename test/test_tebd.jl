using Test, TimeEvoMPS, ITensors
using TimeEvoMPS: isleftortho, isrightortho, measure!

@testset "Basic TEBD tests" begin
    for alg in [TEBD2(), TEBD2Sweep(), TEBD4()]
        @testset "$alg" begin
            # check that evolving with identity doesn't change the state
            N=10
            sites = siteinds("S=1/2",N)
            trivialH = BondOperator(sites)
            psi0 = randomMPS(sites)
            orthogonalize!(psi0,1)
            psi = deepcopy(psi0)
            tebd!(psi,trivialH,0.01,1., alg, progress=false)

            @test inner(psi0,psi) ≈ 1

            #check that evolving forward and backward in time brings us back to the same state
            #(up to numerical errors of course)
            J,h = 0.3, -0.7
            H = tfi_bondop(sites,J,h)

            psi0 = randomMPS(sites)
            psi = deepcopy(psi0)
            tebd!(psi,H,0.01,1.,alg, progress=false)
            tebd!(psi,H,-0.01,-1.,alg, progress=false)

            @test inner(psi0,psi) ≈ 1

            #check that bond dimension is growing to maximum during evolution
            psi = productMPS(sites,ones(Int,N))
            tebd!(psi,H,0.01,5.,alg, progress=false)
            @test maxlinkdim(psi) == 2^5

        end
    end

    #check that the state stays normalized with truncation
    N = 20
    J,h = 0.3, -0.7
    sites = siteinds("S=1/2",N)
    H = tfi_bondop(sites,J,h)
    psi = productMPS(sites, fill("↑",N))
    tebd!(psi,H,0.1,20,TEBD2(),maxdim=5, progress=false)

    @test isapprox(inner(psi,psi), 1., atol=1e-10)

end

@testset "Imaginary time-evolution TFI model" begin
    N=10
    J = 1.
    h = 0.5
    sites = siteinds("S=1/2",N)
    for alg in [TEBD2(), TEBD2Sweep(), TEBD4()]
        @testset "$alg" begin
            psi = productMPS(sites,ones(Int,N))
            H = tfi_bondop(sites,J,h)
            Hgates = gates(H)

            Es = []
            for dt in [0.1,0.01]
                tebd!(psi,H,-1im*dt,-1im*500*dt, alg ; maxdim=50, cutoff=1e-8, orthogonalize=10, progress=false)
                push!(Es, measure(gates(H),psi))
            end
            # exact expression for ground-state energy
            # at criticality (ref? took it from ITensors.jl tests)

            eexact = 0.25 -0.25/sin(π/(4*N + 2))
            @test Es[end] ≈ eexact atol=1e-4
        end
    end
end

@testset "Compare TEBD2 and TEBD4" begin
    N=10
    J,h = 1., 5.87
    dt, tf = 0.01,2
    sites= siteinds("S=1/2",N)
    psi2 = productMPS(sites,ones(Int,N))
    psi4 = complex!(deepcopy(psi2))
    H = tfi_bondop(sites,J,h)

    tebd!(psi2,H,dt,tf,cutoff=1e-12, progress=false)
    tebd!(psi4,H,dt,tf,TEBD4(),cutoff=1e-12, progress=false)

    # since error in TEBD2 is O(t*dt²) while
    # error in TEBD4 is O(t*dt⁴) the most we can
    # hope for is an agreement up to O(t*dt²)
    @test inner(psi2,psi4) ≈ 1 atol= dt^2

    #measure magnetizations and make sure they are the same
    sz2 = measure!(psi2,"Sz")
    sz4 = measure!(psi4,"Sz")
    @test maximum(abs.(sz2 - sz4)) < dt^2

    # we can also check that the Trotter error for TEBD4 scales like dt^4
    # by comparing TEBD4 with step size dt and TEBD2 with step size dt^2
    # (since N is very small we don't expect truncation errors kicking in here)
    psi4 = productMPS(sites,ones(Int,N))
    dt4 = 0.1
    tebd!(psi4,H,dt4,tf,TEBD4(), progress=false)
    @test inner(psi2,psi4) ≈ 1 atol=dt^2
    sz4 = measure!(psi4,"Sz")
    @test maximum(abs.(sz2 - sz4)) < dt^2
end

using Test, ITensors, TimeEvoMPS
using TimeEvoMPS: measure!

@testset "Basic TDVP tests" begin
    N=10
    sites= siteinds("S=1/2",N)

    #check that evolving forward and backward in time brings us back to the same state
    #(up to numerical errors of course)
    J,h = 1., 0.5
    H = tfi_mpo(J,h,sites)

    psi0 = productMPS(ComplexF64,sites, [rand()>0.5 ? "↑" : "↓" for i in 1:length(sites)])
    psi = deepcopy(psi0)
    nsteps = 20
    dt = 1e-1
    tdvp!(psi,H,dt,nsteps*dt, hermitian=true,exp_tol=1e-12/dt,cutoff=1e-14, progress=false)
    tdvp!(psi,H,-dt,-nsteps*dt, hermitian=true, exp_tol=1e-12/dt, cutoff=1e-14, progress=false)

    @test inner(psi,psi) ≈ 1
    @test inner(psi0,psi) ≈ 1

    #check that energy is approximately conserved up to time-step error
    #for short time evolution with no truncation
    psi0 = productMPS(ComplexF64,sites, [rand()>0.5 ? "↑" : "↓" for i in 1:length(sites)])

    nsteps = 20
    dt = 1e-1
    E0 = inner(psi,H,psi)

    tdvp!(psi,H,dt,nsteps*dt, hermitian=true, progress=false)
    E1 = inner(psi,H,psi)
    @test E0 ≈ E1

    #check that the state stays normalized with truncation
    N = 20
    sites = siteinds("S=1/2",N)
    H = tfi_mpo(J,h,sites)
    psi = productMPS(ComplexF64,sites, fill("↑",N))
    tdvp!(psi,H,1,10,maxdim=5,hermitian=true, progress=false)

    @test isapprox(inner(psi,psi), 1., atol=1e-10)

end

@testset "compare short-time evolution TDVP and TEBD" begin
    N=10
    J,h = 1., 5.87
    dt, tf = 0.01,0.5
    sites= siteinds("S=1/2",N)
    psi_tebd = productMPS(sites,ones(Int,N))
    psi_tdvp = complex!(deepcopy(psi_tebd))

    H = tfi_mpo(J,h,sites)
    tdvp!(psi_tdvp,H,dt,tf, cutoff=1e-12, hermitian=true, progress=false)

    H = tfi_bondop(sites,J,h)
    tebd!(psi_tebd,H,dt,tf,cutoff=1e-12,hermitian=true, progress=false)

    # problem is that the error in TDVP scales as
    # Δt^4 while in TEBD2 it scales as Δt^2
    # so I guess we should be ok with difference which is
    # O(Δt^2) ?
    @test inner(psi_tdvp,psi_tebd) ≈ 1 atol= dt^2

    #measure magnetizations and make sure they are the same
    sz_tebd = measure!(psi_tebd,"Sz")
    sz_tdvp = measure!(psi_tdvp,"Sz")
    @test maximum(abs.(sz_tebd - sz_tdvp)) < dt^2
end

# TODO: compare some time-dependent observable with some known results
# TODO: find some analytical results to compare to

@testset "Imaginary time-evolution TFI model" begin
    N=10
    J = 1.
    h = 0.5
    sites = siteinds("S=1/2",N)
    H = tfi_mpo(J,h,sites)
    psi = productMPS(sites,ones(Int,N))

    Es = []
    nsteps = 100
    for dt in [0.1]
        tdvp!(psi,H,-1im*dt,-1im*nsteps*dt ; maxdim=50, hermitian=true, exp_tol = 1e-14/dt, progress= false)
        push!(Es, inner(psi,H,psi))
    end

    # exact expression for ground-state energy
    # at criticality (ref? took it from ITensors.jl tests)
    eexact = 0.25 -0.25/sin(π/(4*N + 2))
    @test Es[end] ≈ eexact atol=1e-4
end

@testset "compare free-fermions" begin
    N = 20
    dt = 0.1
    tf = 2
    sites = siteinds("Fermion",N, conserve_nf = true)
    ampo = AutoMPO()
    for b in 1:N-1
        add!(ampo,1,"Cdag",b, "C",b+1)
        add!(ampo,-1,"C",b,"Cdag",b+1)
    end
    H =MPO(ampo,sites)
    psi = productMPS(sites, [i%2==0 ? "Emp" : "Occ" for i in 1:N])
    cb = LocalMeasurementCallback(["N"], sites,0.1)
    tdvp!(psi,H,dt,tf, cutoff=1e-10,callback=cb, progress=false)

    ns_exec = free_fermions_densities(N,0.1,tf)
    ns = measurements(cb)["N"]

    for i in 1:length(ns)
        @test all( abs.(ns[i] - ns_exec[i+1]).< 5e-5)
    end
end

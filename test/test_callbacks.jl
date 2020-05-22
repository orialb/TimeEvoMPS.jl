using Test, ITensors, TimeEvoMPS
const te = TimeEvoMPS

@testset "NoTEvoCallback" begin
    obs = NoTEvoCallback()
    sites = siteinds("S=1/2",10)
    psi = randomMPS(sites)
    @test te.apply!(obs,psi; t=0,bond=1,sweeend=true) == nothing
    @test te.checkdone!(obs,psi) == false
    @test te.callback_dt(obs) == 0
end


@testset "LocalMeasurementCallback" begin
    sites = siteinds("S=1/2",10)
    psi = randomMPS(sites)
    cb = te.LocalMeasurementCallback(["Sz","Sx"], sites,0.1)
    @test te.callback_dt(cb) ==0.1
    @test length(te.measurement_ts(cb))==0

    # make sure that no measurements are performed at wrong times
    for (t,b,se) in [(0.03,1,true), (0.1,5,true), (0.1,1,false)]
        te.apply!(cb, psi, t=t, bond = b, sweepend=se)
        @test length(te.measurement_ts(cb))==0
    end

    te.apply!(cb, psi, t=0.1, bond = 1, sweepend=true)
    @test length(te.measurement_ts(cb))==1
    @test psi.llim==0
    @test psi.rlim==2
    for o in ["Sz","Sx"]
        @test length(te.measurements(cb)[o][end])==length(psi)
    end

    N=10
    H = tfi_bondop(sites,1.0,1.0)
    psi = productMPS(sites, fill("↑",N))
    cb = te.LocalMeasurementCallback(["Sz","Sx"], sites,0.5)
    tebd!(psi,H,0.1,5,TEBD2(),maxdim=5,callback=cb)

    @test length(te.measurement_ts(cb))==length(0.5:0.5:5)
    for i in 1:length(psi)
        orthogonalize!(psi,i)
        m = scalar(dag(psi[i])*noprime(op(sites, "Sz", i)*psi[i]))
        @test te.measurements(cb)["Sz"][end][i] ≈ m
    end

    H = tfi_mpo(1.0,1.0,sites)
    psi = productMPS(sites, fill("↑",N))
    cb2 = te.LocalMeasurementCallback(["Sz","Sx"], sites,0.5)
    tdvp!(psi,H,0.1,5,maxdim=5,hermitian=true,callback=cb2)

    @test length(te.measurement_ts(cb2))==length(0.5:0.5:5)
    for i in 1:length(psi)
        orthogonalize!(psi,i)
        m = scalar(dag(psi[i])*noprime(op(sites, "Sz", i)*psi[i]))
        @test te.measurements(cb2)["Sz"][end][i] ≈ m
    end

end



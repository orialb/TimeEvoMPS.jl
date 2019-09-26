using Test, ITensors, TimeEvoMPS
const te = TimeEvoMPS

@testset "NoTEvoObserver" begin
    obs = NoTEvoObserver()
    sites = spinHalfSites(10)
    psi = randomMPS(sites)
    @test te.observe!(obs,psi) == nothing
    @test te.checkdone!(obs,psi) == false
    @test te.measurement_step(obs) == 0
end

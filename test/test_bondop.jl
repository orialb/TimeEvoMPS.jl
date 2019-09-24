using Test, TimeEvoMPS, ITensors
const te = TimeEvoMPS

@testset "bond operator" begin
    sites = spinHalfSites(10)
    H = BondOperator(sites)

    add!(H,1im,"Sx",1)
    @test op(sites,siteterms(H,1)[1]) == 1im*op(sites,"Sx",1)

    add!(H,1.,"Sz","Sx",1)
    bt = bondterms(H,1)[1]
    @test te.leftop(sites,bt) == op(sites,"Sz",1)
    @test te.rightop(sites,bt) == op(sites,"Sx",2)

end

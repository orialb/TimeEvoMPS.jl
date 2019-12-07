using Test, TimeEvoMPS, ITensors
const te = TimeEvoMPS


@testset "SiteTerm and BondTerm" begin
    sites = siteinds("S=1/2",10)
    st = SiteTerm(1,1,"Sx*Sx")
    @test op(sites,st)*4 == op(sites,"Id",1)
    st = SiteTerm(1,4,"Sx*Sz*Sx")
    @test op(sites,st) == -1*op(sites,"Sz",1)
end

"""
utility function for testing gates created from BondOperator
"""
function TFIsingBondGate(sites, b, J,h)
    fac = i -> (i==1 || i==length(sites) ? 1 : 1/2)
    G = J*op(sites,"Sz",b)*op(sites,"Sz",b+1) +
            h*(fac(b)*op(sites,"Sx",b)*op(sites,"Id",b+1) + fac(b+1)*op(sites,"Id",b)*op(sites,"Sx",b+1))
    return BondGate(G,b)
end

@testset "bond operator" begin
    sites = siteinds("S=1/2",10)
    H = BondOperator(sites)

    add!(H,1im,"Sx",1)
    @test op(sites,siteterms(H,1)[1]) == 1im*op(sites,"Sx",1)

    add!(H,1.,"Sz","Sx",1)
    bt = bondterms(H,1)[1]
    @test te.leftop(sites,bt) == op(sites,"Sz",1)
    @test te.rightop(sites,bt) == op(sites,"Sx",2)

    #compare against manually-built operator for TF-Ising model
    J = -1
    h = 0.7
    H = BondOperator(sites)
    for b in 1:length(H)-1
        add!(H,J,"Sz","Sz",b)
        add!(H,h,"Sx",b)
    end
    add!(H,h,"Sx",length(H))

    for b in 1:length(H)-1
        @test bondgate(H,b) == TFIsingBondGate(sites,b,J,h)
    end
end

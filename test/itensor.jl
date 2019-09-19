using Test, TimeEvoMPS, ITensors
const te = TimeEvoMPS

@testset "exponentiate" begin
    s1 = Index(2,"s1")
    s2 = Index(2,"s2")
    Amat = rand(2,2,2,2)
    A = ITensor(Amat, prime(s1),prime(s2),s1,s2)

    Aexp = te.exp(A,te.findprimeinds(A,1))
    Amatexp = reshape( exp(reshape(Amat,4,4)), 2,2,2,2)
    Aexp_from_mat = ITensor(Amatexp, prime(s1),prime(s2),s1,s2)
    @test Aexp ≈ Aexp_from_mat



    #test that exponentiation works when indices need to be permuted
    Aexp = te.exp(A,te.findprimeinds(A,0))
    Amatexp = Array( exp(  reshape(Amat,4,4))' )
    Aexp_from_mat = ITensor(reshape(Amatexp,2,2,2,2), s1,s2,prime(s1),prime(s2))
    @test Aexp ≈ Aexp_from_mat

    @test_throws DimensionMismatch te.exp(A,IndexSet(s1))

    #TODO: test a more general case (e.g. indices dimensions differ from each other)
end



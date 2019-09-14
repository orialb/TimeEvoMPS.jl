using Test, TimeEvoMPS, ITensors
const te = TimeEvoMPS

@testset "exponentiate" begin
    s1 = Index(2,"s1")
    s2 = Index(2,"s2")
    Amat = rand(2,2,2,2)
    A = ITensor(Amat, prime(s1),prime(s2),s1,s2)

    Aexp = te.exp(A,te.findprimeinds(A,1), te.findprimeinds(A,0))
    Amatexp = exp(reshape(Amat,4,4))
    Aexp_from_mat = ITensor(reshape(Amatexp,2,2,2,2), prime(s1),prime(s2),s1,s2)
    @test Aexp ≈ Aexp_from_mat
end

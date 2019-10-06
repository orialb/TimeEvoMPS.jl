using ITensors, BenchmarkTools, KrylovKit

struct LocalOp
    A::ITensor
    B::ITensor
end

(L::LocalOp)(t::ITensor) = noprime(L.A*t*L.B)

m=100
i,j = Index(m,"i"), Index(m,"j")
v = randomITensor(i,j)
A= randomITensor(i',i)
B = randomITensor(j,j')
M =LocalOp(A/norm(A),B/norm(B) )

@btime exponentiate($M,0.1,$v,tol=1e-12)

using ITensors
using KrylovKit: exponentiate
using BenchmarkTools


struct LocalOp
    A::ITensor
    B::ITensor
end

function product(o::LocalOp,
                 v::ITensor)::ITensor
    Hv = v;
    Hv *= o.A;
    Hv *= o.B;
    return noprime(Hv)
end

(L::LocalOp)(t::ITensor) = product(L,t)


function exp_random_op()
    m=100
    i,j = Index(m,"i"), Index(m,"j")
    v = randomITensor(i,j)
    A= randomITensor(i',i)
    B = randomITensor(j,j')
    M =LocalOp(A/norm(A),B/norm(B) )

    @btime exponentiate($M,0.1,$v,tol=1e-10)
end


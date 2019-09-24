# some extensions for ITensor functionality
# the goal is to eventually contribute these upstream
import ITensors: _permute_for_factorize
import LinearAlgebra.exp

function findprimeinds(is::IndexSet, plevel::Int=-1)
    if plevel>=0
        return IndexSet(filter(x->plev(x)==plevel, is.inds))
    else
        return IndexSet(filter(x->plev(x)>0, is.inds) )
    end
end

findprimeinds(A::ITensor,args...) = findprimeinds(A.inds,args...)

function isleftortho(M,i)
    i==length(M) && return true
    L = M[i]*prime(dag(M[i]), i==1 ? "Link" : commonindex(M[i],M[i+1]))
    l = linkindex(M,i)
    return norm(L-delta(l,l')) < 1E-12
end

function isrightortho(M,i)
    i==1 && return true
    R = M[i]*prime(dag(M[i]),i==length(M) ? "Link" : commonindex(M[i],M[i-1]))
    r = linkindex(M,i-1)
    return norm(R-delta(r,r')) < 1E-12
end


function storage_exp(As::Dense{T}, Lis,Ris) where {T}
    expAdata = LinearAlgebra.exp( reshape(data(As),dim(Lis),dim(Ris)) )
    return Dense{T}(vec(expAdata))
end

"""
    exp(A::ITensor, Lis::IndexSet)
compute the exponent of the tensor `T` by treating as a matrix ``T_{lr}`` with
the left index `l` running over all indices in `Lis` and `r` running over all
indices not in `Lis`. Must have `dim(Lis) == dim(inds(A))/dim(Lis)` for the exponentiation to
be defined.
"""
function exp(A::ITensor, Lis::IndexSet)
    (dim(Lis) == dim(inds(A))/dim(Lis)) || throw(DimensionMismatch("dimension of the left index set `Lis` must be
                                                                     equal to `dim(inds(A))/dim(Lis)`"))
    A, Lis, Ris = _permute_for_factorize(A,Lis)
    expAs = storage_exp(store(A), Lis,Ris)
    return ITensor(inds(A),expAs)
end


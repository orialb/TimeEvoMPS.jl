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


function storage_exp(As::Dense{T}, Lis,Ris) where {T}
    expAdata = LinearAlgebra.exp( reshape(data(As),dim(Lis),dim(Ris)) )
    return Dense{T}(vec(expAdata))
end

"""
    exp(A::ITensor, Lis::IndexSet, Ris::IndexSet)
compute the exponent of the tensor `T` by treating as a matrix ``T_{lr}`` with
the left index `l` running over all indices in `Lis` and `r` running over all
indices in `Ris`. Must have `dim(Lis) == dim(Ris)` for the exponentiation to
be defined.
"""
function exp(A::ITensor, Lis::IndexSet, Ris::IndexSet)
    dim(Lis) == dim(Ris) || throw("dimension of the left index set `Lis` must be
                                    equal to dimension of the right index set `Ris`")
    A, Lis, Ris = _permute_for_factorize(A,Lis)
    expAs = storage_exp(store(A), Lis,Ris)
    return ITensor(inds(A),expAs)
end


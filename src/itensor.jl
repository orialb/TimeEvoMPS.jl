# some extensions for ITensor functionality
# the goal is to eventually contribute these upstream

function findprimeinds(is::IndexSet, plevel::Int=-1)
    if plevel>=0
        return IndexSet(filter(x->plev(x)==plevel, is.inds))
    else
        return IndexSet(filter(x->plev(x)>0, is.inds) )
    end
end

findprimeinds(A::ITensor,args...) = findprimeinds(A.inds,args...)

"""
    exp(T::ITensor, Lis::IndexSet, Ris::IndexSet)
compute the exponent of the tensor `T` by treating as a matrix ``T_{lr}`` with
the left index `l` running over all indices in `Lis` and `r` running over all
indices in `Ris`. Must have `size(Lis) == size(Ris)` for the exponentiation to
be defined.
"""
function exp(T::ITensor, Lis::IndexSet, Ris::IndexSet)
    size(Lis) == size(Ris) || throw("size of the left index set `Lis` must be
                                    equal to size of the right index set `Ris`")
end

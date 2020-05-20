# some extensions for ITensors functionality
# the goal is to eventually contribute these upstream if found appropriate

function findprimeinds(is::IndexSet, plevel::Int=-1)
    if plevel>=0
        return filter(x->plev(x)==plevel, is)
    else
        return filter(x->plev(x)>0, is)
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

function reorthogonalize!(psi::MPS)
    ITensors.setleftlim!(psi,-1)
    ITensors.setrightlim!(psi,length(psi)+2)
    orthogonalize!(psi,1)
    psi[1] /= sqrt(inner(psi,psi))
end



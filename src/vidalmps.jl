export VidalMPS, Smat

"""
    VidalMPS
A wrapper around `ITensors.MPS` to represent
MPS in Vidal's form where we keep track of the singular values.
"""
struct VidalMPS
    psi::MPS
    S_::Vector{ITensor}
end

function VidalMPS(psi::MPS)
    S_ = Vector{ITensor}()
    for b in reverse(1:length(psi)-1)
        orthogonalize!(psi,b)
        li = linkindex(psi,b)
        if dim(linkindex(psi,b)) > 1
            _,S,_ = svd(psi[b]*psi[b+1],inds(psi[b]) )
            #TODO : using storage interface here, should be a way
            # to do that in a more robust way (maybe Vector(S) for DiagITensor)
            push!(S_,diagITensor(data(store(S)), IndexSet(li',li)))
        else
            push!(S_,diagITensor([1.], IndexSet(li',li)))
        end
    end
    # orthogonalize!(psi,1)
    return VidalMPS(psi,reverse(S_))
end

@forward VidalMPS.psi Base.getindex,Base.setindex!, ITensors.maxLinkDim,ITensors.siteindex, ITensors.linkindex, Base.length
@forward VidalMPS.psi Base.lastindex

ITensors.orthogonalize!(psi::VidalMPS,b) = throw("VidalMPS does not support shifting orthogonality center")

Smat(psi::VidalMPS,b) = psi.S_[b]
setSmat!(psi::VidalMPS,b,S::ITensor) = (psi.S_[b] = S)

function debugortho(M,i)
    i==1 && return 0
    R = M[i]*prime(dag(M[i]),i==length(M) ? "Link" : commonindex(M[i],M[i-1]))
    r = linkindex(M,i-1)
    return norm(R-delta(r,r'))
end


function apply_gate!(psi::VidalMPS,G::BondGate ; kwargs...)
    b = bond(G)
    wf_ = noprime( G*(psi[b]*psi[b+1]) )
    wf = b>1 ? noprime(Smat(psi,b-1)*wf_ ) : wf_
    U,S,V,u,v= svd(wf, inds(psi[b]); kwargs...)
    linktag = tags(linkindex(psi,b))
    replacetags!(V,tags(v) , linktag)
    replacetags!(U,tags(u) , linktag)

    normS = sqrt(scalar(S*S))

    psi[b] = wf_*dag(V) / normS
    psi[b+1] = V
    li = linkindex(psi,b)
    setSmat!(psi,b,diagITensor(data(store(S)) ./ normS, IndexSet(li',li)))
    return nothing
end

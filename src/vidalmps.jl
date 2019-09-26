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
        if dim(linkindex(psi,b)) > 1
            _,S,_ = svd(psi[b]*psi[b+1],inds(psi[b]) )
            #TODO : using storage interface here, should be a way
            # to do that in a more robust way (maybe Vector(S) for DiagITensor)
            push!(S_,diagITensor(data(store(S)), linkindex(psi,b)))
        else
            push!(S_,diagITensor([1.], linkindex(psi,b)))
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

function apply_gate!(psi::VidalMPS,G::BondGate ; kwargs...)
    b = bond(G)
    wf_ = noprime( G*(psi[b]*psi[b+1]) )
    wf = b>1 ? Smat(psi,b-1)*wf_ : wf_
    U,S,V = svd(wf, inds(psi[b]); kwargs...)
    linktag = tags(linkindex(psi,b))
    replacetags!(V,tags(findinds(V,"Link,v")[1]) , linktag)
    #FIX: this is not working currently, need to renormalize
    psi[b] = wf_*V
    psi[b+1] = V
    @assert isrightortho(psi,b)
    @assert isrightortho(psi,b+1)
    setSmat!(psi,b,ITensor(data(store(S)), linkindex(psi,b)))
    return nothing
end

export BondGate, measure, apply_gate!, apply_gates!, exp

import Base: * , ==
import LinearAlgebra.exp

#TODO: need to think if this type really adds something
# or it would be better to just work with an ITensor and extract
# the bond from the site indices.
struct BondGate
    G::ITensor
    b::Int
end

bond(G::BondGate) = G.b

*(G::BondGate, wf::ITensor) = G.G*wf
*(wf::ITensor,G::BondGate) = wf*G.G
*(c::Number, G::BondGate) = BondGate(c*G.G,bond(G))
*(G::BondGate,c::Number) = c*G
==(G1::BondGate,G2::BondGate) = (G1.G == G2.G && G1.b == G2.b)

const GateList = Vector{BondGate}

exp(G::BondGate; ishermitian=false) = BondGate(exp(G.G, findprimeinds(G.G), findprimeinds(G.G,0);
                                                    ishermitian=ishermitian), G.b)

"""
    apply_gate!(psi::MPS,G::BondGate ; kwargs...)
Apply a two-site bond gate `G` to the MPS `psi` and update the
MPS tensors `psi[b],psi[b+1]` (with `b` the bond on which `G` acts).

# Keyword arguments:
All keyword accepted by `ITensors.replacebond!`.
The important ones are:
- `ortho::String`: Options are ["left", "right"].
    If "left" orthogonality center will be at `b` after update, if "right"
    orthogonality center will be at `b+1`
    and `psi[b+1]= sqrt(S)*V'`.
- `maxdim::Int`: If specified, keep only `maxdim` largest singular values after applying gate.
- `mindim::Int`: Minimal number of singular values to keep if truncation is performed according to
    value specified by `cutoff`.
- `cutoff::Float`: If specified, keep the minimal number of singular-values such that the discarded weight is
    smaller than `cutoff` (but bond dimension will be kept smaller than `maxdim`).
- `use_absolute_cutoff::Bool = false`: If `true` truncate all singular-values whose square is smaller than `cutoff`.
"""
function apply_gate!(psi::MPS,G::BondGate ; kwargs...)
    b = bond(G)
    orthogonalize!(psi,b)
    wf = psi[b]*psi[b+1]
    wf = noprime( G*(psi[b]*psi[b+1]) )
    spec = replacebond!(psi, b, wf; normalize=get(kwargs, :normalize,true), kwargs...)
    return spec
end

function apply_gates!(psi::MPS, Gs::Vector{BondGate}; kwargs...)
    rev = get(kwargs,:reversed,false)
    cb_func = get(kwargs,:cb, nothing)
    Gs_ordered = rev ? (@view Gs[end:-1:1]) : Gs
    for g in Gs_ordered
        spec=apply_gate!(psi,g; kwargs...)
        !isnothing(cb_func) && cb_func(psi; bond=bond(g),spec=spec)
    end
end

"""
    measure(H::GateList,psi::MPS)
measure the expectation value of the operator
defined by ``\\sum_i H_i`` with respect to an
MPS `psi`.
"""
function measure(H::GateList,psi::MPS)
    energy = 0
    for (b,h) in enumerate(H)
        orthogonalize!(psi,b)
        wf = psi[b]*psi[b+1]
        energy += scalar(dag(wf)*noprime( h*wf))
    end
    return energy
end


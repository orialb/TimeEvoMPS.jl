struct BondGate
    G::ITensor
    b::Int
end

bond(G::BondGate) = G.b

import Base: *
*(G::BondGate, wf::ITensor) = G.G*wf

const BondOperator = Vector{BondGate}

exp(G::BondGate, dt; hermitian=false) = BondGate( exp(dt*G.G, findprimeinds(G.G); hermitian=hermitian), G.b)

####
## BondOperator
###
export SiteTerm, coeff

struct SiteTerm
    i::Int
    coeff::Number
    opname::String
end

coeff(sop::SiteTerm) = sop.coeff

function ITensors.op(sites, sop::SiteTerm)
    ops = split(sop.opname,"*")
    o = op(sites,String(ops[1]),sop.i)
    for j in 2:length(ops)
        o *= prime(op(sites,String(ops[j]),sop.i))
        o = mapprime(o,2,1)
    end
    return coeff(sop)*o
end


struct BondTerm
    b::Int
    coeff::Number
    leftop::SiteTerm
    rightop::SiteTerm
end

coeff(bop::BondTerm) = bop.coeff
leftop(sites,bop::BondTerm) = op(sites,bop.leftop,b)
rightop(bop::BondTerm) = op(sites,bop.rightop,b)

"""
    BondHamiltonian
A representation of a generic operator which is composed only of one- and two-site terms
"""
struct BondHamiltonian
    sites::SiteSet
    bondops::Dict{Int,Vector{BondTerm}}
    siteops::Dict{Int,Vector{SiteTerm}}
end


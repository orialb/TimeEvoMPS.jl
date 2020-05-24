export add!, BondOperator, bondgate,
       gates, siteterms, bondterms

####
## BondOperator
###
export SiteTerm, coeff, BondTerm

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

    BondTerm(b::Int,c::Number,op1::String,op2::String) =
        new(b,c,SiteTerm(b,1,op1),SiteTerm(b+1,1,op2))
end

coeff(bop::BondTerm) = bop.coeff
leftop(sites,bop::BondTerm) = op(sites,bop.leftop)
rightop(sites,bop::BondTerm) = op(sites,bop.rightop)

"""
    BondOperator
A representation of a generic operator which is composed of a sum of one- and two- site terms
"""
struct BondOperator
    sites::Vector{<: Index}
    bondterms::Dict{Int,Vector{BondTerm}}
    siteterms::Dict{Int,Vector{SiteTerm}}

    BondOperator(sites::Vector{<: Index}) =
        new(sites, Dict(b=> Vector{BondTerm}() for b in 1:length(sites)-1),
            Dict(i => Vector{SiteTerm}() for i in 1:length(sites)))
end

function Base.show(io::IO, bo::BondOperator)
    println(io,"BondOperator")
    print(io,"Site indices: ")
    show(io,bo.sites)
end

Base.length(bo::BondOperator) = length(bo.sites)

#TODO : add some error messages when trying to access site or bond out of range
siteterms(bo::BondOperator,i) = bo.siteterms[i]
bondterms(bo::BondOperator,b) = bo.bondterms[b]

ITensors.add!(bo::BondOperator, op::String, i::Int) = add!(bo,1.,op,i)
ITensors.add!(bo::BondOperator, c::Number, op::String, i::Int) = (push!(bo.siteterms[i], SiteTerm(i,c,op)))
ITensors.add!(bo::BondOperator, op1::String, op2::String, b::Int) = add!(bo,1.,op1,op2,b)
ITensors.add!(bo::BondOperator,c::Number, op1::String, op2::String, b::Int) = (push!(bo.bondterms[b], BondTerm(b,c,op1,op2)))

function bondgate(bo::BondOperator, b::Int)
    sites = bo.sites
    gate = ITensor(sites[b],sites[b]',sites[b+1], sites[b+1]')
    for t in bondterms(bo,b)
        gate += (coeff(t)*leftop(sites,t)*rightop(sites,t))
    end
    for i in [b,b+1]
        fac = i==length(bo) || i==1 ? 1 : 1/2
        f = (i,o) -> fac*(i==b ? op(sites,o)*op(sites,"Id",b+1) : op(sites,"Id",b)*op(sites,o) )
        for st in siteterms(bo,i)
            gate += f(i,st)
        end
    end
    return BondGate(gate,b)
end

gates(bo::BondOperator) = (x->bondgate(bo,x)).(1:length(bo)-1)




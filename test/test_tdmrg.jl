using Test, TimeEvoMPS, ITensors
using TimeEvoMPS: isleftortho, isrightortho

function TFIsingBondGate(sites::SiteSet, b, J,h)
    fac = i -> (i==1 || i==length(sites) ? 1 : 1/2)
    G = J*op(sites,"Sz",b)*op(sites,"Sz",b+1) +
        h*(fac(b)*op(sites,"Sx",b)*op(sites,"Id",b+1) + fac(b+1)*op(sites,"Id",b)*op(sites,"Sx",b+1))
    return BondGate(G,b)
end

TFIsingBondOperator(sites::SiteSet, J,h) = map(x->TFIsingBondGate(sites,x,J,h), 1:length(sites)-1)

function noevolution_bondop(sites::SiteSet)
    gates = Vector{BondGate}()
    for b in 1:length(sites)-1
        s1,s2 = sites[b], sites[b+1]
        d = dim(s1)
        push!(gates, BondGate( ITensor(zeros(d,d,d,d), s1',s2',s1,s2), b))
    end
    return gates
end

function tfimMPO(sites::SiteSet,
                 h::Float64)
  # Input operator terms which define a Hamiltonian
  ampo = AutoMPO(sites)
  for j=1:length(sites)-1
    add!(ampo,-1.,"Sz",j,"Sz",j+1)
    add!(ampo,h,"Sx",j)
  end
  add!(ampo,h,"Sx",length(sites))
  # Convert these terms to an MPO tensor network
  return toMPO(ampo)
end

function TFIgs(sites,h)
    psi0 = randomMPS(sites)

    # define parameters for DMRG sweeps
    sweeps = Sweeps(15)
    maxdim!(sweeps, 10,20,100,100,200)
    cutoff!(sweeps, 1E-10)
    H = tfimMPO(sites,h)
    energy, psi = dmrg(H,psi0,sweeps)
    return psi
end


function measure(H::BondOperator,psi)
    energy = 0
    for (b,h) in enumerate(H)
        orthogonalize!(psi,b)
        wf = psi[b]*psi[b+1]
        energy += scalar(dag(wf)*noprime( h*wf))
    end
    return energy
end

N=10
sites = spinHalfSites(N)
psi = randomMPS(sites)
H = TFIsingBondOperator(sites,1,1)

# Exact energy for transverse field Ising model
# with open boundary conditions at criticality
energy_exact = 0.25 - 0.25/sin(π/(2*(2*N+1)))

Es = []
for dt in [0.01,1e-3,1e-4]
    tdmrg!(psi,H,1im*dt,1im*250*dt ; maxdim=50)
    orthogonalize!(psi,length(psi))
    orthogonalize!(psi,1)
    psi[1] /= sqrt(inner(psi,psi))
    push!(Es, measure(H,psi))
end

inner(psi,psi)

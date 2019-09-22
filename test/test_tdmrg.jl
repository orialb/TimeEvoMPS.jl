using Test, TimeEvoMPS, ITensors

function TFIsingBondGate(sites::SiteSet, b, J,h)
    fac = i -> (i==1 || i==length(sites) ? 1 : 1/2)
    G = J*op(sites,"Sz",b)*op(sites,"Sz",b+1) +
        h*(fac(b)*op(sites,"Sx",b)*op(sites,"Id",b+1) + fac(b+1)*op(sites,"Id",b)*op(sites,"Sx",b+1))
    return BondGate(G,b)
end

TFIsingBondOperator(sites::SiteSet, J,h) = map(x->TFIsingBondGate(sites,x,J,h), 1:length(sites)-1)

N=10
sites = spinHalfSites(N)
psi = randomMPS(sites)
H = TFIsingBondOperator(sites,-1,1)
# tdmrg!(psi,H,-0.01*1im,-10im ; maxdim=50)
# tdmrg!(psi,H,-0.001*1im,-1im ; maxdim=50)
tdmrg!(psi,H,0.01,0.02 ; maxdim=50)



# Exact energy for transverse field Ising model
# with open boundary conditions at criticality
energy_exact = 0.25 - 0.25/sin(π/(2*(2*N+1)))

energy = 0
for (b,h) in enumerate(H)
    orthogonalize!(psi,b)
    wf = psi[b]*psi[b+1]
    global energy += scalar(wf*noprime( h*wf))
end
energy

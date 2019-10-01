# TimeEvoMPS (WIP)
[![Build Status](https://travis-ci.org/orialb/TimeEvoMPS.jl.svg?branch=master)](https://travis-ci.org/orialb/TimeEvoMPS.jl)
[![codecov](https://codecov.io/gh/orialb/TimeEvoMPS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/orialb/TimeEvoMPS.jl)

The goal of this package is to provide implementations of time-evolution algorithms for matrix-product states using 
[ITensors.jl](https://github.com/ITensor/ITensors.jl). 

The package is currently at a very initial stage. Contributions and suggestions are very welcome. 

Algorithms currently implemented:
- TEBD with 2nd order Trotter decomposition [1]
- TDVP (two-site variant) [2]


## Installation
Since ITensors.jl is not yet a registered package you will have to make sure it is installed before installing TimeEvoMPS.jl.
In the Julia REPL:
```
] add https://github.com/ITensor/ITensors.jl 
```
After you have installed ITensors.jl you can go ahead and install TimeEvoMPS:
```
] add https://github.com/orialb/TimeEvoMPS.jl
```

# Usage
The following code example shows how to evolve an MPS for a spin-half chain with the transverse-field Ising Hamiltonian, starting from a fully polarized state (functionality to perform measurements during time evolution is still missing, but will be added very soon).

```julia
using ITensors, TimeEvoMPS

N=10
J = 1.
h = 0.5

# Use ITensors to define initial fully-polarized MPS
sites = spinHalfSites(N)
psi = productMPS(sites, fill("↑",N))

# Build Hamiltonian
H = BondOperator(sites)
for b in 1:N-1
    #add two-site term at bond b
    add!(H,J,"Sz","Sz",b)
    #add single-site term at site b
    add!(H,h,"Sx",b)
end
add!(H,h,"Sx",N)

# set maximal bond dimension during evolution
maxdim = 10
#time step and total evolution time 
dt = 0.01
tf =1.
#evolve
tebd!(psi,H,dt,tf, maxdim=maxdim)
``` 

## References
[1] Vidal, G. (2004). Efficient Simulation of One-Dimensional Quantum Many-Body
Systems. Physical Review Letters, 93(4), 040502.
https://doi.org/10.1103/PhysRevLett.93.040502

[2] Haegeman, J., Lubich, C., Oseledets, I., Vandereycken, B., & Verstraete,
F. (2016). Unifying time evolution and optimization with matrix product states.
Physical Review B, 94(16). https://doi.org/10.1103/PhysRevB.94.165116

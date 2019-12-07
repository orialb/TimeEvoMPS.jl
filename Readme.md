# TimeEvoMPS 
[![Build Status](https://travis-ci.org/orialb/TimeEvoMPS.jl.svg?branch=master)](https://travis-ci.org/orialb/TimeEvoMPS.jl)
[![codecov](https://codecov.io/gh/orialb/TimeEvoMPS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/orialb/TimeEvoMPS.jl)

The goal of this package is to provide implementations of time-evolution algorithms for matrix-product states using 
[ITensors.jl](https://github.com/ITensor/ITensors.jl). 

The package is currently at an initial stage. Contributions and suggestions are very welcome. 

Algorithms currently implemented:
- TEBD (with 2nd and 4th order Trotter decomposition) [1]
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

First we setup the initial state 
```julia
using ITensors, TimeEvoMPS

N=10
J = 1.
h = 0.5

# Use ITensors to define initial fully-polarized MPS
sites = siteinds("S=1/2",N)
psi = productMPS(sites, fill("↑",N))
```

## TEBD

Define our Hamiltonian as a `BondOperator` which is an object that represents 
an operator composed of a sum of two-site gates H = Σᵢ h_{i,i+1} (calling `gates(bo::BondOperator)` builds
the actual gates).

```julia
# Build Hamiltonian
H = BondOperator(sites)
for b in 1:N-1
    #add two-site term at bond b
    add!(H,-J,"Sz","Sz",b)
    #add single-site term at site b
    add!(H,-h,"Sx",b)
end
add!(H,-h,"Sx",N)
```

Now we can run TEBD, by default 2nd order Trotter decomposition is used.
```julia
# Set maximal bond dimension during evolution.
# Other truncation parmeters supported by ITensors.jl are possible, see 
# documentation of `apply_gate!`
maxdim = 10
#time step and total evolution time 
dt = 0.01
tf =1.
#evolve
tebd!(psi,H,dt,tf, maxdim=maxdim)
``` 

Alternatively we could use 4th order Trotter decomposition
```julia
tebd!(psi,H,dt,tf,TEBD4(), maxdim=maxdim)
```

## TDVP

We could also use TDVP to evolve our state. For this we first need to build an MPO representation
of the Hamiltonian. This is easily done using the AutoMPO functionality from ITensors.

```julia
ampo = AutoMPO()
for j=1:length(sites)-1
    add!(ampo,-J,"Sz",j,"Sz",j+1)
    add!(ampo,-h,"Sx",j)
end
add!(ampo,-h,"Sx",length(sites))
H= MPO(ampo,sites)
```

Now we can run time-evolution using TDVP. 
```julia
psi = complex!(psi)
tdvp!(psi,H,dt,tf,maxdim=maxdim)
```

( Note that currently for real-time evolution the ITensors stored 
in the MPS `psi` must be of type `Complex` (this will hopefully not be necessary in the future).
I couldn't find a way to initialize an ITensors.MPS with `ComplexF64` (probably this will be available in the future)
so I wrote a quick hack (the `complex!` function) to convert a `Float` MPS to a `Complex` MPS just by multiplying each 
ITensor by complex unity. This is just a temporary thing.)

## References
[1] Vidal, G. (2004). Efficient Simulation of One-Dimensional Quantum Many-Body
Systems. Physical Review Letters, 93(4), 040502.
https://doi.org/10.1103/PhysRevLett.93.040502

[2] Haegeman, J., Lubich, C., Oseledets, I., Vandereycken, B., & Verstraete,
F. (2016). Unifying time evolution and optimization with matrix product states.
Physical Review B, 94(16). https://doi.org/10.1103/PhysRevB.94.165116

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
In order to use this package you will need to have ITensors.jl installed, run the following in the Julia REPL:
```
] add ITensors
```
This package is not yet registered, in order to install it run the following in the Julia REPL:
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
psi = productMPS(ComplexF64,sites, fill("↑",N))
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
tdvp!(psi,H,dt,tf,maxdim=maxdim)
```

## Callbacks (i.e. performing measurements during time-evolution) (WIP)

Probably you will be interested in measuring some observables at different time
points during the evolution. This is possible using the callback mechanism which
allows you to access the state of the system at each point of the evolution
(including steps in the middle of the TEBD/TDVP sweep). Callbacks can also be
used to stop the time-evolution if certain criteria have been met (e.g.
convergence of some observables). You could use one of the implemented callbacks
or implement a callback of your own (the callback mechanism is still work in
progress, so the API might change).

### Local measurements callback

Currently a `LocalMeasurmentCallback` is implemented which allows measuring single-site operators. 
Here is how you would use it to measure the expectation values of `Sz,Sx,Sy`, with measurement interval `dt=0.2` 
(time evolving with the MPO `H` defined in the previous section):
```julia
psi = productMPS(sites, fill("↑",N))
cb = LocalMeasurementCallback(["Sz","Sx","Sy"], sites,0.2)
tdvp!(psi,H,0.1,5,maxdim=30,callback=cb)
```
note that the measurment interval must be commensurate with the time-evolution step (here `0.1`).

Now you can extract the measurement results from the callback calling `measurements(cb)`. For example here is how you could plot
the measurements, using e.g. PyPlot.jl: 
```julia
using PyPlot
ts = measurement_ts(cb)
for o in ["x","y","z"]
    S5 = getindex.(measurements(cb)["S$o"],5)
    plot(ts,S5,"-o",label="\$S^$(o)_5\$")
end
legend()
xlabel("t")
```

### Implementing a custom callback (TODO)

## References
[1] Vidal, G. (2004). Efficient Simulation of One-Dimensional Quantum Many-Body
Systems. Physical Review Letters, 93(4), 040502.
https://doi.org/10.1103/PhysRevLett.93.040502

[2] Haegeman, J., Lubich, C., Oseledets, I., Vandereycken, B., & Verstraete,
F. (2016). Unifying time evolution and optimization with matrix product states.
Physical Review B, 94(16). https://doi.org/10.1103/PhysRevB.94.165116

# TimeEvoMPS (WIP)
The goal of this package is to provide implementations of time-evolution algorithms for matrix-product states using 
[ITensors.jl](https://github.com/ITensor/ITensors.jl). 

The package is currently at a very initial stage. Contributions and suggestions are very welcome. 
Currently only TEBD using 2nd order Trotter decomposition is implemented. 

Note: since ITensors.jl is not yet a registered package you will have to make sure it is installed before installing TimeEvoMPS.jl. 

# TODO
- [ ] implement observer to perform measurements during time evolution
- [ ] TEBD with 4th order Trotter
- [ ] TDVP

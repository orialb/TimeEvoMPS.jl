These benchmarks are currently only for the purpose of sanity check, to verify
that there aren't any big mistakes in the implementations and that performance is 
reasonable. 

I compared timings between TimeEvoMPS, implementations from the TenPy library, and the c++ version 
of iTensor. So far it seems that performance of all three is similar when bond dimension is large enough.
Let me stress again that the benchmarks I did so far are not very rigorous 
as I am not 100% sure that I managed to reproduce exactly the same computation in all of the above mentioned libraries so you should 
not interpret the fact that there are a few seconds difference as one package being faster than the other.


For the Julia benchmarks, I used Julia 1.2 compiled with MKL. 
The iTensor (c++) was also compiled with MKL, and TenPy was run with numpy compiled with MKL as well.

I am using MKL_NUM_THREADS=2 for the benchmark.

```
Platform Info:
  OS: macOS (x86_64-apple-darwin17.7.0)
  CPU: Intel(R) Core(TM) i7-5557U CPU @ 3.10GHz
```

# TEBD

Evolve with critical transverse-field Ising model, starting with fully polarized state. 
Require `maxdim=100, cutoff=1e-14` (hopefully I managed to choose the equivelent parameters in TenPy) 
for a chain of size `N=20`. Run with `dt=0.1` up to `t=10`.

For a more robust benchmark I should probably prepare first a state with maximal desired bond dimension 
(using TEBD or random unitaries) and than benchmark a few sweeps on it. 

- TenPy (2nd order trotter) :  26.5 s
- TimeEvoMPS (2nd order trotter) : 22.621 s

Same as above but with a different type of 2nd order Trotter decomposition (TEBD2Sweep):
- iTensor (2nd order trotter left-right sweep) : 22.114 s
- TimeEvoMPS (2nd order trotter left-right sweep) :  36.77 s

This seems like too big of a difference might be a hint to some performance issue in TEBD2Sweep, although it makes
sense that TEBD2Sweep will take roughly twice than TEBD2 because there are roughly twice as many gate applications. 

Same as above but with `maxdim=100` and `t=5`:
- iTensor : 56.895s
- TimeEvoMPS : 46.385s

This hints that there might be some differences in the way truncation is done (or I am using some parameters wrong) in the two examples. 
In general I would expect the c++ code to be a bit faster. 

# TDVP 
Here I used `MKL_NUM_THREADS=1`.

Times are per sweep after saturation to maximal bond dimension , `N=20, maxdim=100`. I restricted 
the maximal Krylov dimension to `krylovdim=20` and used `tol=1e-10` (which was the default tolerance setting in the
ITensor version). I don't know if those are good defaults but I was just trying to have a fair comparison.
After performing the benchmarks I realized that a possible mistake is that maybe the `tol` parameter in `KrylovKit.exponentiate` does
not have exactly the same meaning as `ErrGoal` in `applyExp` in ITensor, so it might be that after all I did not do a fair comparison.

- ITensor (tdvp branch) ~ 14.4
- TimeEvoMPS ~ 8.7

When using `krylovdim=30`:

- ITensor ~  20.2
- TimeEvoMPS ~ 13.9

I Also looked at timings for the Krylov exponentiation methods themselves, using a random projected local
operator acting on a two-legged tensor with both indices having dimension 100, I get (caveat above still applies)

- ITensor ~ 32.762ms
- KrylovKit ~ 22.525ms

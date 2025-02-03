# In-Medium Similarity Renormalization Group (IMSRG)

Functions for the In-Medium Similarity Renormalization Group (IMSRG) method.

Under the current implementation, the IMSRG is the most time-consuming method among the ones available in the package.
One should consider to run a script using multiple threads to speed up the calculations like the following example:

```bash
julia -t 12 path_to_script.jl
```

Some functions in the package are already parallelized, but it may be still insufficient for large systems.

```@autodocs
Modules = [PairingHamiltonian]
Pages = ["imsrg.jl"]
``` 


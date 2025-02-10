# How to start 

The main API of the package is the `PairingHamiltonian` function in `src/hamiltonian.jl`.

## Example

One may use the following script to run the Hamiltonian for different methods and values of the pairing strength `gval`.

```julia
# using release version
using PairingHamiltonian

# # using a local version for development
# include("src/PairingHamiltonians.jl")
# using .PairingHamiltonians

methods = ["Full-CI(2-fold)", "HF", "BCS", "CCD" , "IMSRG(2)"]
Norb = 20
Nocc = 10 

debug_mode = 0
gvals = collect(-1.0:0.1:1.0)

for gval in gvals
    println("gval = $gval")
    for method in methods
        @timeit to "$method" begin
            Eret = PairingHamiltonian(;Norb_in=Norb, Nocc_in=Nocc, gval=gval,
                                        debug_mode=debug_mode, solver=method)
        end
    end
end
```

For `Norb=20` and `Nocc=10`, the results will be like:

![](Energies_Norb20_Nocc10.png)



## Optional arguments

One can specify the following optional arguments:

| Argument | Description | Default |
|:---------|:------------|:--------|
| `Norb_in::Int` | number of orbitals | 8 |
| `Nocc_in::Int` | number of occupied states | 4 |
| `gval::Float64` | pairing strength | 0.33 |
| `delta_eps::Float64` | energy difference between orbitals | 1.0 |
| `debug_mode::Int` | specify the debug mode (can be 0, 1, 2) | 0 |
| `solver::String` | method to solve the Hamiltonian | "FCI(2-fold)" |
| `save_Exact_wf::Bool` | save the full-CI wave functions as HDF5 files | false |
| `to` | TimerOutput object, if defined in a user script | nothing |

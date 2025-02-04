# Many-Body Perturbation Theory

One can use the Many-Body Perturbation Theory (MBPT) to calculate the ground state energy of a system.
In the package, PT2 and PT3 are implemented, and can be used to calculate the ground state energy of a system.
Users do not have to specify "MBPT2" or "MBPT3" in the input file, as the package will automatically use the methods when the Hartree-Fock calculation is called.


```@autodocs
Modules = [PairingHamiltonians]
Pages = ["mbpt.jl"]
``` 


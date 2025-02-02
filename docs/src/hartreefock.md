# Hartree-Fock

The Hartree-Fock (HF) method is based on the assumption that the many-body wave function can be approximated by a single Slater determinant. 
The HF is nothing but a naive filling for the global pairing Hamiltonian, but it is a good starting point for more advanced methods such as MBPT, post-HF methods (CC, IMSRG, etc.).
For this reason, some functions are not necessary for the HF method, but they are implemented for future use.

```@autodocs
Modules = [PairingHamiltonian]
Pages = ["hartreefock.jl"]
``` 


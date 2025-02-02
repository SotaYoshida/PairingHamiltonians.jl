# Hamiltonian

Functions for constructing paring Hamiltonian and some utilities for handling it.

```math
H_\mathrm{pair} = \sum_{i=1}^{N} \epsilon_i a_i^\dagger a_i - \frac{g}{4} \sum_{ij} a_i^\dagger a^\dagger_{\bar{i}} a_{\bar{j}} a_j
```

One can use the following methods to solve the pairing Hamiltonian with the package:
- "Full-CI": Full-CI/Exact diagonalization
- "Full-CI(2-fold)": Full-CI/Exact diagonalization assuming that particles are paired
- "HF": Hartree-Fock (nothing but naive filling due to the pairing Hamiltonian)  
   Under this method, MBPT calculation with the HF reference state are carried out.
- "BCS": Bardeen-Cooper-Schrieffer
- "CCD": Coupled Cluster, especially CCD
- "IMSRG": In-Medium Similarity Renormalization Group, especially IM-SRG(2)


```@autodocs
Modules = [PairingHamiltonian]
Pages = ["hamiltonian.jl"]
``` 


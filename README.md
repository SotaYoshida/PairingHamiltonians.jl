# PairingHamiltonian.jl

Julia package solving pairing hamiltonian

The package covers the following methods:

- Full Configuration Interaction (FCI)/Exact Diagonalization
- Hartree-Fock (HF)
- Bardeen-Cooper-Schrieffer (BCS)
- Many-Body Perturbation Theory (MBPT) a.k.a. Moeller-Plesset method
- Coupled Cluster (CC) 
- In-Medium Similarity Renormalization Group (IM-SRG)
- Eigenvector Continuation (EC)

## Installation

Assuming you have Julia installed, you can install the package by running the following commands in the Julia REPL:
```julia
] add PairingHamiltonian
``` 

## How to start

If you cloned the repository, you are ready to run a simple example in the `examples` folder. Just run the following commands in terminal:

```bash
```julia -t 12 examples/sample_script.jl
```

where `12` is the number of threads. Note that the only IM-SRG part has some parallelization.


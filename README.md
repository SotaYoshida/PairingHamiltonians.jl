# PairingHamiltonian.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://SotaYoshida.github.io/NuclearToolkit.jl/dev)
[![Build Status](https://github.com/SotaYoshida/PairingHamiltonian.jl/actions/workflows/CI.yml/badge.svg?branch=dev)](https://github.com/SotaYoshida/PairingHamiltonian.jl/actions/workflows/CI.yml/badge.svg?branch=dev)

<img src="https://github.com/SotaYoshida/PairingHamiltonian.jl/blob/main/logo/logo_PairingHamiltonian.png?raw=true" width=60%>

Julia package to solve pairing Hamiltonian.

The package covers the following methods:

- Full Configuration Interaction (FCI)/Exact Diagonalization
- Hartree-Fock (HF)
- Bardeen-Cooper-Schrieffer (BCS)
- Many-Body Perturbation Theory (MBPT) a.k.a.  Møller–Plesset method
- Coupled Cluster (CC) 
- In-Medium Similarity Renormalization Group (IM-SRG)
- Eigenvector Continuation (EC)

## Installation

Note: This package is still under development. The following instructions are not working yet.

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


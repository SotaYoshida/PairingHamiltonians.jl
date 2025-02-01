# PairingHamiltonian

Julia package to solve pairing Hamiltonian with various methods.

Upon the request to give a lecture in a summer school in Japan, I decided to develop a Julia package to solve the pairing Hamiltonian with various methods.
The package itself is not covering quantum computing methods, but this can be used as a reference to play with the quantum computing methods.

Regarding the quantum computing methods, I have developed and opened many materials to solve the pairing Hamiltonian with quantum algorithms, VQE on PennyLane/Qiskit/Pytket, Hadamard test, QPE, Quantum Lanczos, and so on.



## Installation and example

First, prepare Julia environment.

Second, add the package in Pkg mode of the Julia REPL.
```julia
julia>]add PairingHamiltonian
```

Then, you can use the package as follows.
```julia
using PairingHamiltonian
```

That's it!

## Package features 

PairingHamiltonian.jl is a Julia package to solve the pairing Hamiltonian with a bunch of many-body machinery.
The package is covering the following methods:

- Full-CI/Exact diagonalization 
- Hartree-Fock
- Bardeen-Cooper-Schrieffer (BCS)
- Many-Body Perturbation Theory
- Coupled Cluster
- In-Medium Similarity Renormalization Group
- Eigenvector Continuation

To author's knowledge, this is the first open-source software that provides all these methods in a single package.
Once can use the package to dive into the many-body physics and compare the results of different methods, learn the methods, and develop new ones.

## Issues/Pull requests

This package is designed to be an open-source software and to guarantee reproducibility and transparency of the future works.
Making issues and pull requests are welcome. More details are in the [CONTRIBUTING.md]().

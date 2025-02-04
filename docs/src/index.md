# PairingHamiltonians.jl

Julia package to solve pairing Hamiltonian with a bunch of many-body methods used in nuclear physics.

## Why I Created PairingHamiltonians.jl

I started to use Julia in 2018, and I was fascinated by the performance, the simplicity of the language, and 
the philosophy of the Julia language, such as the following:
[Why We Created Julia](https://julialang.org/blog/2012/02/why-we-created-julia)

Upon the request to give a lecture in a summer school in Japan, I decided to develop a Julia package and open it to the public.
The pairing Hamiltonian is a simple model, but it is very important in nuclear physics, condensed matter physics, and provides a good starting point to learn methods to solve many-body problems, quantum computing, and so on.

Nuclear many-body methods today are very diverse and rich. This may originate from a long history of nuclear physics and rich phenomena in nuclei, which have led to the development of many methods to describe the nuclear structure and reactions from different perspectives.
This is a good thing and very fascinating. I really love the diversity of the methods and the rich phenomena in nuclei.
However, it may be difficult for beginners like Bachelor students or first-year graduate students to understand the whole picture of the nuclear many-body methods.

I was one such student, to be honest. I was struggling to understand the many-body methods in nuclear physics in the first year of the graduate school. (I am still struggling, but I am enjoying it now.)
This package is designed to help a student (like old me) to dive into the many-body methods in nuclear physics.

I would be very happy if this package provides a good starting point for students to learn the many-body methods in nuclear physics, and to develop something new.

The package itself is not covering quantum computing methods, but this can be used as a reference to play with the quantum computing methods.
Regarding the quantum computing methods, I have developed and opened some materials to solve the pairing Hamiltonian with quantum algorithms, VQE, Hadamard test, QPE, Quantum Lanczos with either PennyLane/Qiskit/Pytket.
Those materials will be available through the following repository: [Lecture_SummerSchool2025](https://github.com/SotaYoshida/Lecture_SummerSchool2025)
This is for the summer school in Japan, so the materials are written in Japanese, but I am planning to make an English version in the future.

I have been developing another Julia package, [NuclearToolkit.jl](https://github.com/SotaYoshida/NuclearToolkit.jl) covering chiral EFT potentials, and methods like ones in this package.
After getting used to methods provided in this package, you can use the NuclearToolkit.jl to solve more realistic nuclear physics problems.

## Installation and example

First, prepare Julia environment.

Second, add the package in Pkg mode of the Julia REPL.
```julia
julia>]add PairingHamiltonians
```

Then, you can use the package as follows.
```julia
using PairingHamiltonians
```

That's it. You are ready to use the package.

## Package features 

PairingHamiltonians.jl is a Julia package to solve the pairing Hamiltonian with a bunch of many-body machinery.
The package is covering the following methods:

- Full-CI/Exact diagonalization 
- Hartree-Fock
- Bardeen-Cooper-Schrieffer (BCS)
- Many-Body Perturbation Theory
- Coupled Cluster
- In-Medium Similarity Renormalization Group
- Eigenvector Continuation

To author's knowledge, this is the first open-source software that provides all these methods in a single package.
One can use the package to get used to the many-body methods, and to build a new method by combining or modifying the existing methods.

## Issues/Pull requests

This package is designed to be an open-source software, and the author is welcoming any contributions.

Especially, the author would like to get feedback from students who are using this package to learn the many-body methods in nuclear physics.
Of course, contributions from other fields are highly welcome.
Making issues and pull requests are welcome. More details are in the [Contributing to PairingHamiltonians.jl]().

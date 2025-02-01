module PairingHamiltonian

using Base.Threads
using Combinatorics
using Glob
using HDF5
using LinearAlgebra
using Random
using Printf
using StatsBase
using TimerOutputs

include("hamiltonian.jl")
export main_pairHamil
export eval_Hamil
include("fullCI.jl")
include("bcs.jl")
include("hartreefock.jl")
include("mbpt.jl")
include("coupled_cluster.jl")
include("imsrg.jl")
include("eigenvector_continuation.jl")
export main_EC_from_FCI

end # module ParingHamiltonian

using Documenter
using PairingHamiltonians
## For development in local machine
#include("../src/PairingHamiltonians.jl")
#push!(LOAD_PATH,"../src/")

DocMeta.setdocmeta!(PairingHamiltonians, :DocTestSetup, :(using PairingHamiltonians); recursive=true)
makedocs(;
         modules=[PairingHamiltonians],
         authors="SotaYoshida <syoshida@cc.utsunomiya-u.ac.jp>",
         repo="https://github.com/SotaYoshida/PairingHamiltonians.jl/blob/{commit}{path}#{line}",
         sitename="PairingHamiltonians.jl",
         format=Documenter.HTML(;
                                prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://SotaYoshida.github.io/PairingHamiltonians.jl",
                                assets=String[],
                                ),
         warnonly = [:missing_docs],
         pages=[
            "Home" => "index.md",
            "How to start" => "howtostart.md",
            "Contributing to PairingHamiltonians" => "contributing.md",
            "Hamiltonian" => "hamiltonian.md",
            "Many-body methods" => [
                "Full CI" => "fullCI.md",
                "Hartree-Fock" => "hartreefock.md",
                "MBPT" => "mbpt.md",
                "BCS" => "bcs.md",                
                "Coupled Cluster" => "coupledcluster.md",
                "IM-SRG" => "imsrg.md",
                "Eigenvector Continuation" => "eigenvectorcontinuation.md",
            ]
    ],
)

deploydocs(;
    repo="github.com/SotaYoshida/PairingHamiltonians.jl.git",
    devbranch="dev",
)

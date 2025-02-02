using Documenter
#using PairingHamiltonian
include("../src/PairingHamiltonian.jl")
push!(LOAD_PATH,"../src/")

DocMeta.setdocmeta!(PairingHamiltonian, :DocTestSetup, :(using PairingHamiltonian); recursive=true)
makedocs(;
         modules=[PairingHamiltonian],
         authors="SotaYoshida <syoshida@cc.utsunomiya-u.ac.jp>",
         repo="https://github.com/SotaYoshida/PairingHamiltonian.jl/blob/{commit}{path}#{line}",
         sitename="PairingHamiltonian.jl",
         format=Documenter.HTML(;
                                prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://SotaYoshida.github.io/PairingHamiltonian.jl",
                                assets=String[],
                                ),
         warnonly = [:missing_docs],
         pages=[
            "Home" => "index.md",
            "Contributing to PairingHamiltonian" => "contributing.md",
            "Hamiltonian" => "hamiltonian.md",
            "Many-body methods" => [
                "Full-CI" => "fullCI.md",
                "HartreeFock" => "hartreefock.md",
                "MBPT" => "mbpt.md",
                "BCS" => "bcs.md",                
                "Coupled Cluster" => "coupledcluster.md",
                "IM-SRG" => "imsrg.md",
                "Eigenvector Continuation" => "eigenvectorcontinuation.md",
            ]
    ],
)

deploydocs(;
    repo="github.com/SotaYoshida/PairingHamiltonian.jl.git",
    devbranch="dev",
)

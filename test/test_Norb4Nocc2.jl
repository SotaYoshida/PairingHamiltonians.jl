using PairingHamiltonians
using Test

function run_test()
    Norb = 4; Nocc = 2
    gvals = collect(-0.5:0.1:0.5)
    methods = ["Full-CI(2-fold)", "HF", "BCS", "CCD", "IMSRG(2)"]

    for gval in gvals
        println("--------------- gval = $gval ------------")        
        for method in methods             
            Eret = PairingHamiltonian(Norb_in=Norb, Nocc_in=Nocc, gval=gval, solver=method)
        end
        println("\n\n")
    end
    return true
end

@test run_test()
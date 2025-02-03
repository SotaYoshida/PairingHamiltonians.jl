using PairingHamiltonian
using Test

function run_test()
    Norb = 8; Nocc = 4 
    gvals = collect(-0.5:0.05:0.5)
    methods = ["Full-CI(2-fold)", "HF", "BCS", "CCD", "IMSRG(2)"]

    for gval in gvals
        println("--------------- gval = $gval ------------")        
        for method in methods             
            Eret = Pairing_Hamiltonian(Norb_in=Norb, Nocc_in=Nocc, gval=gval, solver=method)
        end
        println("\n\n")
    end
    return true
end

@test run_test()
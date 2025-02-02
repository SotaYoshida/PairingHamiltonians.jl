#using PairingHamiltonian

# for development
include("src/PairingHamiltonian.jl")
using .PairingHamiltonian

using Plots
using TimerOutputs
using LaTeXStrings

function plot_Es(methods, x, data, Norb, Nocc)
    markers = Dict("Full-CI" => :star5, "Full-CI(2-fold)" => :star5,
                    "HF" => :diamond, "MBPT2" => :utriangle, "MBPT3" => :dtriangle,
                    "CCD" => :circle,
                    "BCS" => :hexagon, "IMSRG(2)" => :octagon
    )
    cols = palette(:seaborn_muted)

    # Energy plot
    yref  = nothing
    p = plot(size =(400,300), legend = :bottom)
    for method in methods
        y = data[method]
        if occursin("Full-CI", method)
            plot!(p, x, y, label= "Full-CI", c=:black)
            ymin = minimum(y); ymax = maximum(y)
            ylims!(p, ( ymin - 0.1*abs(ymin), ymax + 0.1*abs(ymax)))
            yref = y
        else
            plot!(p, x, y, label=method, marker=markers[method],
            markersize=2, alpha=0.8,
            palette=:seaborn_muted
            )
        end
    end
    xlabel!(p, L"g")
    ylabel!(p, "Energy")
    savefig(p, "Energies_Norb$(Norb)_Nocc$(Nocc).pdf")
    
    # Diff. plot
    p = plot(size =(400,300), legend = :bottom)  
    ylims!(p, ( -5.0, 5.0))
    for method in methods
        y = data[method]
        if occursin("Full-CI", method)
            continue
        end
        println("method: $method y $y")
        diff = y .- yref
        plot!(p, x, diff, label=method, marker=markers[method],
        markersize=2, alpha=0.8,
        palette=:seaborn_muted,
        )
    end
    xlabel!(p, L"g")
    ylabel!(p, "Energy difference from Full-CI")
    savefig(p, "Diff_Energies_Norb$(Norb)_Nocc$(Nocc).pdf")

end

function run_test(;debug_mode=0)
    gvals = collect(-0.85:0.01:0.85)

    methods = ["Full-CI(2-fold)", "HF", "BCS", "CCD", "IMSRG(2)"]
    methods_plot = [ "HF", "BCS", "MBPT2", "MBPT3", "CCD", "IMSRG(2)", "Full-CI(2-fold)"]

    Norb = 8; Nocc = 4 
    write_wf = !true 

    data = Dict{String, Vector{Float64}}()
    to = TimerOutput()
    for gval in gvals
        println("--------------- gval = $gval ------------")
        for method in methods 
            tf_write_wf = false            
            if method == "Full-CI(2-fold)" && write_wf
                tf_write_wf = true
            end
            @timeit to "$method" Eret = main_pairHamil(Norb_in=Norb, Nocc_in=Nocc, gval=gval, 
                                                      debug_mode=debug_mode, solver=method,  save_Exact_wf=tf_write_wf,
                                                      to_in=to)
            for key in keys(Eret)
                if !haskey(data, key)
                    data[key] = Float64[]
                end
                push!(data[key], Eret[key])
            end
        end
        println("\n\n")
    end
    @timeit to "plot" plot_Es(methods_plot, gvals, data, Norb, Nocc)

    show(to, allocations=true, compact=true); println("")
end
#run_test()


function for_developing_a_method()
    methods = ["Full-CI(2-fold)", "HF", "BCS", "CCD" , "IMSRG(2)"]

    Norb = 20; Nocc = 10 

    debug_mode = 0
    gvals = [0.33]

    to = TimerOutput()
    for gval in gvals
        println("gval = $gval")
        for method in methods
            @timeit to "$method" Eret = main_pairHamil(;Norb_in=Norb, Nocc_in=Nocc, gval=gval,
                                                       debug_mode=debug_mode, solver=method, save_Exact_wf=false,
                                                       to_in=to)
        end
    end
    show(to, allocations=true, compact=true); println("")
end
for_developing_a_method( )

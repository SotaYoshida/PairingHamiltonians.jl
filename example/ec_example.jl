#using PairingHamiltonian

# for development
include("src/PairingHamiltonian.jl")
using .PairingHamiltonian

using Plots
using LaTeXStrings

function plot_ec(sets, labels, results, Norb, Nocc)
    p = plot(size=(400, 300),
        palette=:seaborn_deep, 
        xlabel=L"g",
        ylabel=L"$E_{\mathrm{EC}} - E_{\mathrm{exact}}$"
    )
    markers = [:circle, :diamond, :utriangle]
    for idx = 1:length(sets)
        label = labels[idx]
        Data = results[idx]

        x = sort([gval for gval in keys(Data)])
        y = zeros(Float64, length(x))
        ydiff = zeros(Float64, length(x))
        for (i, gval) in enumerate(x)
            Eec, Eexact = Data[gval]
            y[i] = Eec
            ydiff[i] = abs(Eec - Eexact)
        end
        plot!(p, x, ydiff, label=label, 
              marker=markers[idx], markersize=3,
              alpha=0.8,
              line=:solid, legend=:topleft)
    end
    savefig(p, "EC_Norb$(Norb)_Nocc$(Nocc).pdf")

end

function example_EC()
    Norb = 8
    Nocc = 4
    set_A = collect(-1.0:0.5:0.0)
    set_B = collect(-1.0:1.0:1.0)
    set_C = collect(0.0:0.5:1.0)

    sets = [set_A, set_B, set_C]
    labels = ["set-A", "set-B", "set-C"]

    results = [ ]
    for target_set in sets
        data = main_EC_from_FCI(Norb, Nocc; gvals_specified=target_set)
        push!(results, data)
    end

    plot_ec(sets, labels, results, Norb, Nocc)

end
example_EC()

struct snapshot
    Norb::Int64
    Nocc::Int64
    gval::Float64
    evec::Vector{Float64}
end

"""
    read_snapshot(fn::String)

Read the snapshot from the file `fn` and return the data as a tuple of `Norb`, `Nocc`, `gval`, and `evec`.
Note that the snapshot is assumed to be written out by `write_Exact_wf_hdf5` function in `src/fullCI.jl`.
"""
function read_snapshot(fn)
    io = h5open(fn, "r")
    Norb = read(io, "Norb")
    Nocc = read(io, "Nocc")
    gval = read(io, "gval")
    evec = read(io, "evec")
    close(io)
    return Norb, Nocc, gval, evec
end

"""
    solve_gen_eig(tildeH, tildeN)

Solve the generalized eigenvalue problem ``\\tilde{H} v = E N v`` and return the minimum eigenvalue ``E``.
Now it is implemented by the `eigen` function in Julia, but it can be ill-conditioned for large systems,
so it is better to use e.g. Cholesky decomposition and iterative solver like `Arpack` or `KrylovKit` for more general cases.
"""
function solve_gen_eig(tildeH, tildeN)
    evals, evecs = eigen(tildeH, tildeN)
    return minimum(evals)
end

"""
    EC_from_FCI(Norb_in, Nocc_in; target_dirs="eigenstates_fullCI", gvals_specified::Vector{Float64} = Vector{Float64}[], gvals_target=collect(-2.0:0.1:2.0), degenerate = true, delta_eps::Float64=1.0)

Main function to compute the energy curve by EC method from the FCI wavefunctions.
Assuming that the FCI wavefunctions are already evaluated and saved in a directory.

# Arguments
- `Norb_in::Int64`: Number of orbitals. 
- `Nocc_in::Int64`: Number of electrons.
- `target_dirs::String`: Directory name where the FCI wavefunctions are saved.
- `gvals_specified::Vector{Float64}`: List of gvals to be used as snapshots. If empty, all the snapshots extracted via `glob` will be used.
- 
"""
function EC_from_FCI(Norb_in, Nocc_in; 
    target_dirs="eigenstates_fullCI", 
    gvals_specified::Vector{Float64} = Vector{Float64}[],
    gvals_target=collect(-2.0:0.1:2.0),
    degenerate = true,
    delta_eps::Float64=1.0
    )
    solver = "Full-CI"


    fns = glob("*Norb$(Norb_in)*Nocc$(Nocc_in)*.h5", target_dirs)
    if length(fns) == 0
        println("No corresponding files found in $target_dirs")
        println("Please run the fullCI calculation first to prepare w.f. with Norb=$Norb_in, Nocc=$Nocc_in")
        println("Aborted.")
        exit(1)
    end

    # Load the FCI wavefunctions
    if length(gvals_specified) == 0
        gvals_specified = [ float(split(split(fn, "_g")[end], ".h5")[1]) for fn in fns ]
        @warn "Since gvals_specified is not specified, all the gvals in the snapshots are to be used. This may not be what you want."
    end
    println(gvals_specified)
    snapshots = snapshot[ ]
    for i = 1:length(fns)
        Norb, Nocc, gval, evec = read_snapshot(fns[i])
        if (gval in gvals_specified) 
            push!(snapshots, snapshot(Norb, Nocc, gval, evec))
        end
    end

    # Sort the snapshots by gval
    sort!(snapshots, by=x->x.gval)
    for i = 1:length(snapshots)
        println("gval = ", snapshots[i].gval)
    end

    # Set up basis, epsilon
    Norb_calc = div(Norb_in, 2)
    Nocc_calc = div(Nocc_in, 2)
    epsilon = eval_epsilon(Norb_calc, delta_eps, degenerate)
    basis = prepare_basis(Norb_calc, Nocc_calc)
    dim_vec = length(basis)

    # To evaluate overlap and Hamiltonian matrix elements
    dim = length(snapshots)
    Hv = zeros(Float64, dim_vec)
    tildeH = zeros(Float64, dim, dim)
    tildeN = zeros(Float64, dim, dim)

    Data = Dict{Float64, Tuple{Float64, Float64}}()
    for g in gvals_target
        # construct Hamiltonian under the specified gval
        dummy, dummy_, Hamil_target = eval_Hamil(solver, basis, epsilon, g, Norb_calc, degenerate)
        # evaluate N and H
        for i = 1:dim
            for j = i:dim
                evec_i = snapshots[i].evec
                evec_j = snapshots[j].evec
                Nij = dot(evec_i, evec_j)
                operate_H_on_vec!(Hv, Hamil_target, evec_j)
                Hij = dot(evec_i, Hv)
                tildeH[i, j] = tildeH[j, i] = Hij
                tildeN[i, j] = tildeN[j, i] = Nij
            end
        end
        E_ec = solve_gen_eig(tildeH, tildeN)

        print(@sprintf("g = %6.2f, E_ec = %12.5f", g, E_ec))
        # comparing to the exact solution if available
        fn_exact = target_dirs * "/eigenstate_Norb$(Norb_in)_Nocc$(Nocc_in)_g$(g).h5"
        if isfile(fn_exact)
            _, __, ___, evec_exact = read_snapshot(fn_exact)
            operate_H_on_vec!(Hv, Hamil_target, evec_exact)
            E_exact = dot(evec_exact, Hv)
            print(@sprintf(", E_exact = %12.5f Diff. = %9.2e", E_exact, E_ec - E_exact), "\n")
        end
        Data[g] = (E_ec, E_exact)
    end
    return Data
end
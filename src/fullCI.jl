"""
    operate_H_on_vec!(w, Hamil_mat::Array{Float64, 2}, v)

Function to compute the matrix-vector product of the Hamiltonian matrix and a vector.
"""
function operate_H_on_vec!(w, Hamil_mat::Array{Float64, 2}, v)
    mul!(w, Hamil_mat, v)
    return nothing
end

"""
    operate_H_on_vec!(w, Hamil_mat::Dict{UInt64, Float64}, v)

Sparse version of the function to compute the matrix-vector product of the Hamiltonian matrix in the form of `Dict{UInt64, Float64}` and a vector.
"""
function operate_H_on_vec!(w, Hamil_mat::Dict{UInt64, Float64}, v)
    w .= 0.0
    for nkey in keys(Hamil_mat)
        val = Hamil_mat[nkey]
        i, j = unhash_key2(nkey)
        w[i] += val * v[j]
        if i != j
            w[j] += val * v[i]
        end
    end
    return nothing
end

"""
    lanczos(Hamil_mat, dim, save_Exact_wf, to; itnum=300, tol=1e-9, debug_mode=0)

Function to compute the lowest eigenvalue of the Hamiltonian using the Lanczos method.

Constructing a Krylov subspace ``\\mathcal{K}_m(H,v) = \\mathrm{span}\\{v, Hv, H^2v, \\cdots, H^{m-1}v\\}``,
and the tridiagonal matrix ``T_m = V_m^T H V_m`` where ``V_m = [v_1, v_2, \\cdots, v_m]`` and ``v_{m+1} = H v_m - \\alpha_m v_m - \\beta_m v_{m-1}``, the Lanczos method iteratively constructs the matrix ``T_m`` and diagonalizes it to obtain the smallest eigenvalue of ``H``.

# Arguments
- `Hamil_mat`: Hamiltonian matrix, either a `Matrix{Float64}` or a `Dict{UInt64, Float64}`.
- `dim`: Dimension of the basis, i.e. number of configurations
- `save_Exact_wf`: If `true`, save the exact wave function to a HDF5 file
- `to`: TimerOutput object

# Optional arguments
- `itnum`: Number of Lanczos iterations
- `tol`: Tolerance for convergence
- `debug_mode`: the level of debug information
"""
function lanczos(Hamil_mat, dim, save_Exact_wf, to; itnum=300, tol=1e-9, debug_mode=0)
    println("Starting Lanczos iteration...")
    Random.seed!(1234)
    v1 = rand(dim)
    v1 ./= norm(v1)
    alpha = beta = 0.0
    T = zeros(Float64, itnum, itnum)
    w = zeros(Float64, dim)
    Es = [1e5, 1e5]
    tE = 1e5
    vks = [zeros(Float64, dim) for _ in 1:itnum]
    vks[2] = v1
    i = 0
    @timeit to "iter" for i = 1:itnum
        v = vks[i+1]
        operate_H_on_vec!(w, Hamil_mat, v)
        alpha = dot(w, v)
        T[i, i] = alpha
        if i > 1
            tE = minimum(eigvals(T[1:i, 1:i]))
            Es[2] = tE
            if debug_mode > 0
                println("iter = ", @sprintf("%6i", i), " tE = ", @sprintf("%12.5f", tE))
            end
            if abs(Es[1] - tE) < tol
                break
            end
        end
        w .-= alpha .* v
        w .-= beta .* vks[i]
        reOrthogonalize!(w, vks, i)
        beta = norm(w)
        vks[i+2] = w ./ beta
        T[i, i+1] = beta
        T[i+1, i] = beta
        if i > 1
            Es[1] = Es[2]
        end
    end
    Ritz_vector = [0.0]
    if save_Exact_wf
        Ritz_vector = make_gs_wf(vks, T, i, 1)
    end
    evec = ifelse(save_Exact_wf, Ritz_vector, [0.0])
    return tE, evec
end

"""
    make_gs_wf(vks, Tmat, Lan_itnum, num_ev)

Function to construct the ground state wave function from the Lanczos vectors.

```math
|\\psi_{\\rm GS}\\rangle = \\sum_{k=1}^{D} c_k |v_k\\rangle
```
where `D` is the number of Lanczos iterations, and ``c_k`` is the eigenvector of the tridiagonal matrix ``T`` corresponding to the smallest eigenvalue.

# Arguments
- `vks`: Lanczos vectors
- `Tmat`: Tridiagonal matrix
- `Lan_itnum`: Number of Lanczos iterations
- `num_ev`: Number of the eigenvector to be used
"""
function make_gs_wf(vks, Tmat, Lan_itnum, num_ev)
    evec = zeros(Float64, length(vks[1]))
    vals, vecs = eigen(@views Tmat[1:Lan_itnum,1:Lan_itnum])
    for k = 1:length(vals)
        evec .+= vecs[k, num_ev] .* vks[k]
    end
    return evec
end

"""
    write_Exact_wf_hdf5(evec, Norb, Nocc, gval)

Function to save the exact wave function to a HDF5 file.

# Arguments
- `evec`: Exact wave function
- `Norb`: Number of orbitals
- `Nocc`: Number of occupied orbitals
- `gval`: Interaction strength
"""
function write_Exact_wf_hdf5(evec, Norb, Nocc, gval)
    if !isdir("eigenstates_fullCI")
        mkdir("eigenstates_fullCI")
    end
    fname = "eigenstates_fullCI/eigenstate_Norb$(Norb*2)_Nocc$(Nocc*2)_g$(gval).h5"
    h5open(fname, "w") do file
        write(file, "Norb", Norb)
        write(file, "Nocc", Nocc)
        write(file, "gval", gval)
        write(file, "evec", evec)
    end
    return nothing
end

"""
    _main_FCI(Hamil_mat, dim_basis, Norb, Nocc, gval, save_Exact_wf, to, debug_mode=0; return_only_E0=true)

Main function to compute the ground state energy with the full CI method.
If the system size is small, the Hamiltonian matrix is explicitly constructed and diagonalized using the `eigen` function in LinearAlgebra.jl.
If the system size is large, the Hamiltonian matrix is sparsely represented and the Lanczos method is used to compute the smallest eigenvalue.

# Arguments
- `Hamil_mat`: Hamiltonian matrix, either a `Matrix{Float64}` or a `Dict{UInt64, Float64}`.
- `dim_basis`: Dimension of the basis, i.e. number of configurations
- `Norb`: Number of orbitals
- `Nocc`: Number of occupied orbitals
- `gval`: Interaction strength
- `save_Exact_wf`: If `true`, save the exact wave function to a HDF5 file
- `to`: TimerOutput object
- `debug_mode`: the level of debug information
- `return_only_E0`: If `true`, return only the ground state energy, otherwise return the eigenvalues and eigenvectors
"""
function _main_FCI(Hamil_mat, dim_basis, Norb, Nocc, gval, save_Exact_wf, to, debug_mode=0; return_only_E0=true)
    sparce_rep = (typeof(Hamil_mat) != Matrix{Float64})
    if !sparce_rep
        @timeit to "solver: eigen" evals, evecs = eigen(Hamil_mat)
        if save_Exact_wf
            write_Exact_wf_hdf5(evecs[:, 1], Norb, Nocc, gval)
        end
        if return_only_E0
            return minimum(evals)
        else
            return evals, evecs
        end
    else
        @timeit to "solver: Lanczos" E0, evec = lanczos(Hamil_mat, dim_basis, save_Exact_wf, to;debug_mode=debug_mode)
        if save_Exact_wf
            write_Exact_wf_hdf5(evec, Norb, Nocc, gval)
        end
        if !return_only_E0
            @error "return_only_E0 must be `true`` when using the Lanczos method"
        end
        return E0
    end
end

"""
    reOrthogonalize!(w, vks, i)

Re-orthogonalize the vector w with respect to the previous vectors:

```math
w := w - \\sum_{j=1}^{i} \\langle w, v_j \\rangle v_j
```
"""
function reOrthogonalize!(w, vks, i)
    for j in 1:i
        w .-= dot(w, vks[j]) .* vks[j]
    end
    return w
end
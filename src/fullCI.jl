function operate_H_on_vec!(w, Hamil_mat::Array{Float64, 2}, v)
    mul!(w, Hamil_mat, v)
    return nothing
end

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
Use the Lanczos method to compute the smallest eigenvalue of the Hamiltonian.
"""
function Lanczos(Hamil_mat, dim, save_Exact_wf, to; itnum=300, tol=1e-9, debug_mode=0)
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

function make_gs_wf(vks, Tmat, Lan_itnum, num_ev)
    evec = zeros(Float64, length(vks[1]))
    vals, vecs = eigen(@views Tmat[1:Lan_itnum,1:Lan_itnum])
    for k = 1:length(vals)
        evec .+= vecs[k, num_ev] .* vks[k]
    end
    return evec
end

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
Compute the eigenvalues of the Hamiltonian by diagonalization.
"""
function get_Egs_diagonalization(Hamil_mat, dim_basis, Norb, Nocc, gval, save_Exact_wf, to, debug_mode=0; return_only_E0=true)
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
        @timeit to "solver: Lanczos" E0, evec = Lanczos(Hamil_mat, dim_basis, save_Exact_wf, to;debug_mode=debug_mode)
        if save_Exact_wf
            write_Exact_wf_hdf5(evec, Norb, Nocc, gval)
        end
    end
    return E0
end

"""
Reorthogonalize for Lanczos methods.
"""
function reOrthogonalize!(w, vks, i)
    for j in 1:i
        w .-= dot(w, vks[j]) .* vks[j]
    end
    return w
end
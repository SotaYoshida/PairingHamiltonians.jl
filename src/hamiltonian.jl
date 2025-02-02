"""
struct for pairing Hamiltonian

# Fields
- `Norb::Int`: Number of orbitals
- `Nocc::Int`: Number of occupied orbitals
- `degenerate::Bool`: Whether assume degenerate orbitals or not. It will be True only when Norb is even and the method is "Full-CI(2-fold)".
- `gval::Float64`: Pairing strength
- `dim_basis::Int`: Dimension of the basis
- `basis::Vector{Int}`: List of configurations
- `delta_eps::Float64`: Energy difference between orbitals
- `epsilon::Vector{Float64}`: Single-particle energies
- `H::Matrix{Float64}`: Hamiltonian matrix
"""
struct paringHamiltonian
    Norb::Int
    Nocc::Int
    degenerate::Bool
    gval::Float64
    dim_basis::Int
    basis::Array{Int,1}
    delta_eps::Float64
    epsilon::Array{Float64,1}
    H::Array{Float64,2}
end

function show_vector(label, vec::Vector{Float64}, step=1)
    println("$label:")
    for i = 1:step:length(vec)
        print(@sprintf("%9.4f", vec[i])," ")
    end
    println()
    return nothing
end

function show_matrix(label, mat::Array{Float64, 2})
    println("$label:")
    Nq = size(mat, 1)
    for i = 1:Nq
        for j = 1:Nq
            print(@sprintf("%9.4f", mat[i, j])," ")
        end
        println()
    end
    return nothing
end

function eval_epsilon(Norb::Int, delta_eps::Float64, degenerate::Bool)
    epsilon = zeros(Float64, Norb)
    for i = 1:Norb
        if degenerate
            epsilon[i] = (i-1) * delta_eps * 2
        else
            epsilon[i] = div(i-1,2) * delta_eps
        end
    end
    return epsilon
end

"""
Create a list of configurations where n out of N are occupied.

"""
function prepare_basis(N, n)
    basis = Int64[]
    for i in 1:2^N
        if count_ones(i) == n
            push!(basis, i)
        end
    end
    @assert length(basis) == binomial(N, n) "# of basis $(length(basis)) != # of comb. $(binomial(N, n))"
    return basis
end

"""
    bitstr2int(bitstr::String)::Int

Convert a bitstring representation (e.g., "00110") to an integer (e.g., 6).
"""
function bitstr2int(bitstr::String)
    N = length(bitstr)
    return parse(Int, bitstr, base=2)
end

"""
    int2bitstr(x::Int, N::Int)::String

Convert an integer to a bitstring of length `N`.
Similar thing can be done with string function in Julia, but this is more proper way for "0" padding. 
"""
function int2bitstr(x, N)
    return string(bitstring(x)[end-N+1:end] |> s -> lpad(s, N, '0'))
end

"""
    check_connectivity(psi_l::Int, psi_r::Int, p::Int, q::Int, Nq::Int)::Bool

Check if the configurations (in bit representation) are connected by the Hamiltonian.
"""
function check_connectivity(psi_l::Int, psi_r::Int, p::Int, q::Int, Nq::Int)
    if p == q
        if int2bitstr(psi_l, Nq)[p] == '1'
            return true
        end
    else
        if int2bitstr(psi_r, Nq)[p] == '0' && int2bitstr(psi_r, Nq)[q] == '1'
            if int2bitstr(psi_l, Nq)[p] == '1' && int2bitstr(psi_l, Nq)[q] == '0'
                return true
            end
        end
    end
    return false
end

"""
    unhash_key2(key::UInt64)::Tuple{Int, Int}
    
Unhash the key to the pair of indices.
""" 
function unhash_key2(key::UInt64)
    i = Int(key >> 32)    
    j = Int(key & 0xffffffff)
    return i, j
end

"""
    hash_key2(i::Int, j::Int)::UInt64

Hash the pair of indices to a key.
"""
function hash_key2(i,j)::UInt64
    key = (UInt64(i) << 32) +  UInt64(j)
    return key
end

"""
    eval_Hamil(method::String, basis::Vector{Int}, eps::Vector{Float64}, g::Float64, Nq::Int, degenerate::Bool; debug_mode::Int=0, sparce_rep::Bool=false)::Tuple{Matrix{Float64}, Matrix{Float64}, Union{Matrix{Float64}, Dict{UInt64, Float64}}}

Function to evaluate the Hamiltonian matrix elements.
Once the number of orbitals and the number of occupied states are given,
one can prepare the Hamiltonian matrix elements in the full matrix representation or sparse representation.
For larger systems, sparse representation using `Dict{UInt64, Float64}` is used instead of storing the full matrix.
The keys of the dictionary correspond to the integer representation of many-body configurations.

Note that `"1"` and `'1'` are different in Julia. The former is a string, while the latter is a character.
"""
function eval_Hamil(method::String, basis::Vector{Int}, eps::Vector{Float64}, g::Float64, 
                    Nq::Int, degenerate::Bool; debug_mode::Int=0, sparce_rep::Bool=false)
    dim = length(basis)
    req = dim^2 * 8 / 1024^3
    mem = Sys.total_memory() / 1024^3
    if req > min(1.0, 0.5 * mem)
        sparce_rep = true
        if debug_mode > 0
            println("Since the memory usage is large $(@sprintf("%5.1e", req)) GB, we use sparse representation.")
        end
    end
   
    Hamil_mat = nothing
    if sparce_rep
        Hamil_mat = Dict{UInt64, Float64}()
    else
        Hamil_mat = zeros(Float64, dim, dim)
    end
    if debug_mode > 0
        println("method $method dim $dim sparse_rep $sparce_rep")
    end
    sum_step = ifelse(degenerate, 1, 2)

    if occursin("Full-CI", method)
        num_nonzero = 0
        for i = 1:dim
            psi = int2bitstr(basis[i], Nq)
            nkey = hash_key2(i,i)
            if sparce_rep
                Hamil_mat[nkey] = 0.0
            end
            for p = 1:sum_step:Nq
                if degenerate 
                    if psi[p] == '1'
                        if sparce_rep
                            Hamil_mat[nkey] += eps[p] - g
                        else
                            Hamil_mat[i,i] += eps[p] - g
                        end
                    end
                else # partner should be occupied to gain pairing energy
                    te = 0.0
                    partner = ifelse(p%2==1, p+1, p-1)
                    occp = (psi[p] == '1')
                    occpp = (psi[partner] == '1')
                    if occp
                        te += eps[p] 
                    end
                    if occpp
                        te += eps[partner]
                    end
                    if occp && occpp
                        te += - g
                    end
                    if sparce_rep
                        Hamil_mat[nkey]  += te
                    else
                        Hamil_mat[i,i] += te
                    end
                end
            end
            num_nonzero += 1

            for j = i+1:dim
                psi_l = basis[i]
                psi_r = basis[j]
                diff_bit = psi_l ⊻ psi_r   
                num_bit_diff = count_ones(diff_bit)
                if (degenerate && num_bit_diff != 2) || (!degenerate && num_bit_diff != 4)
                    continue
                end
                nkey = hash_key2(i,j)
                # pqrs: |0011> → |1100> can contribute
                tmp = 0.0
                if degenerate 
                    for q = 1:Nq
                        for p = q+1:Nq
                            if !check_connectivity(psi_l, psi_r, p, q, Nq)
                                continue
                            end
                            tmp += -g 
                        end
                    end
                else 
                    # a^†_q a^†_qbar a_pbar a_p 
                    bitstr_ket = int2bitstr(psi_l, Nq)
                    bitstr_bra = int2bitstr(psi_r, Nq)
                    diff_bit_str = int2bitstr(diff_bit, Nq)
                    for q = 1:sum_step:Nq
                        qbar = q + 1
                        if diff_bit_str[q] == '0' || diff_bit_str[qbar] == '0'
                            continue
                        end

                        if !(bitstr_ket[q] == bitstr_ket[qbar] == '1')
                            continue
                        end
                        if !(bitstr_bra[q] == bitstr_bra[qbar] == '0')
                            continue
                        end

                        for p = 1:sum_step:q-2# $q-2  
                            pbar = p + 1
                            if !(bitstr_ket[p] == bitstr_ket[pbar] == '0')
                                continue
                            end
                            if !(bitstr_bra[p] == bitstr_bra[pbar] == '1')
                                continue
                            end
                            tmp += -g 
                        end
                    end
                end
                if tmp != 0.0
                    num_nonzero += 2
                    if sparce_rep
                        Hamil_mat[nkey] = tmp
                    else
                        Hamil_mat[i, j] = tmp
                        Hamil_mat[j, i] = tmp
                    end
                end
            end
        end
        if debug_mode > 0
            println("# of nonzero Hamiltonian elements: $num_nonzero / 10^$(@sprintf("%5.2f",log10(dim^2))),  sparcity = $(@sprintf("%6.1e", num_nonzero / dim^2))")
        end
    end
    h1b, h2b, = eval_h1_h2(Nq, eps, g, sum_step)
    return h1b, h2b, Hamil_mat
end

function eval_h1_h2(Nq, eps::Vector{Float64}, g::Float64, sum_step)
    h1b = zeros(Float64, Nq, Nq)
    h2b_dim = div(Nq * (Nq - 1), 2)
    h2b = zeros(Float64, h2b_dim)
    idx = 0 
    for i = 1:Nq
        h1b[i,i] = eps[i]
        for j = i+1:Nq
            idx += 1
            if (i % 2 == 1 && j == i + 1) || (i % 2 == 0 && j == i - 1)
                h2b[idx] = -g
            end
        end
    end
    return h1b, h2b
end

"""
    main_pairHamil(to; Norb_in::Int=8, Nocc_in::Int=4, gval::Float64=0.33, delta_eps::Float64=1.0,
                        debug_mode::Int=0, solver::String="FCI(2-fold)", save_Exact_wf::Bool=false)

Main function to evaluate the ground state energy of the pairing Hamiltonian.

# Arguments
- `to::TimerOutput`: timer object to measure the elapsed time

# Optional arguments
- `Norb_in::Int(8)`: number of orbitals
- `Nocc_in::Int(4)`: number of occupied states
- `gval::Float64(0.33)`: pairing strength
- `delta_eps::Float64(1.0)`: energy difference between orbitals
- `debug_mode::Int(0)`: specify the debug mode
- `solver::String("FCI(2-fold)")`: method to solve the Hamiltonian. This can be one of "Full-CI", "Full-CI(2-fold)", "HF", "BCS", "CCD", "IMSRG(2)". If "HF" is chosen, PT2/PT3 energies are also calculated.
- `save_Exact_wf::Bool(false)`: save the full-CI wave functions as HDF5 files, which can be used for e.g. analysis of the wave functions and constructing surrogate models like eigenvector continuation.
"""
function main_pairHamil(to; Norb_in::Int=8, Nocc_in::Int=4, gval::Float64=0.33, delta_eps::Float64=1.0,
                        debug_mode::Int=0, solver::String="Full-CI(2-fold)", save_Exact_wf::Bool=false)
    # working on single particle basis or pair basis
    degenerate = false
    if occursin("2-fold", solver)
        degenerate = true
    end

    # evaluate epsilon (single-particle energies)
    epsilon = eval_epsilon(Norb_in, delta_eps, degenerate)

    # In what follows, Norb&Nocc are the number of pairs if degenerate is true. 
    Norb = Norb_in; Nocc = Nocc_in
    if degenerate
        Norb = div(Norb_in, 2)
        Nocc = div(Nocc_in, 2)
    end

    # prepare base
    basis = prepare_basis(Norb, Nocc)
    if debug_mode > 1
        println("basis for $Norb, $Nocc:")
        for bit_int in basis
            bit_str = int2bitstr(bit_int, Norb)
            println(bit_str)
        end
    end

    # evaluate dim_basis
    dim_basis = length(basis)

    # evaluate H
    @timeit to "eval-H" begin
        h1b, h2b, Hmat = eval_Hamil(solver, basis, epsilon, gval, Norb, degenerate; debug_mode=debug_mode)        
        println("Norb $Norb Nocc $Nocc dim_basis $dim_basis")
    end

    if debug_mode > 1
        print("\tChecking bitstr2int for given basis")
        for bit_int in basis
            bit_str = int2bitstr(bit_int, Norb)
            bitstr2int(bit_str) == bit_int || error("bitstr2int failed")
        end
        println("... PASSED!")
    end

    E0 = 0.0
    Eret = Dict{String, Float64}()
    @timeit to "solve" begin
        if occursin("Full-CI", solver)
            E0 = _main_FCI(Hmat, dim_basis, Norb, Nocc, gval, save_Exact_wf, to, debug_mode)
            println("E(Full-CI) = $E0")
            Eret[solver] = E0
        elseif solver == "HF" || solver == "CCD" || solver == "IMSRG(2)"
            E0, EPT2, EPT3, HNO, holes, particles = _main_HF(Nocc, h1b, h2b, gval, to, debug_mode)
            F = HNO.f
            if solver == "HF"
                println("E(HF) = $E0, E(HF+PT2) = $(E0 + EPT2) E(HF+PT2+PT3) = $(E0 + EPT2 + EPT3)")
                Eret[solver] = E0
                Eret["MBPT2"] = E0 + EPT2
                Eret["MBPT3"] = E0 + EPT2 + EPT3
            elseif solver == "CCD"
                ECCD = _main_CC(F, gval, Nocc, to, debug_mode)                 
                Eret[solver] = E0 + ECCD
                println("E(CCD) = $(E0+ECCD)")
            elseif solver == "IMSRG(2)"
                E_IMSRG = _main_IMSRG(HNO, holes, particles, gval, Nocc, to, debug_mode) 
                Eret[solver] = E_IMSRG
                println("E(IMSRG(2)) = $(E_IMSRG)")
            end
        elseif solver == "BCS"
            E0 = NaN 
            if gval > 0.2
                E0 = _main_BCS(Norb, Nocc, epsilon, gval, to, debug_mode)
                println("E(BCS) = $E0")
            end
            Eret[solver] = E0
        else
            @error "Unknown solver: $solver"
        end
    end
    return Eret
end



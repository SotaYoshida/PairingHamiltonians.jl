struct ket2b
    i::Int64
    j::Int64
    pp::Bool
    hh::Bool
    ph::Bool
end

struct Operator
    zerobody::Vector{Float64}
    onebody::Matrix{Float64}
    twobody::Matrix{Float64}
end 

struct Mat222ph
    Mat_Nab::Matrix{Float64}
    Mat_Npp::Matrix{Float64}
    Mat_Nhh::Matrix{Float64}
    eta_ph::Matrix{Float64}
    Γ_ph::Matrix{Float64}
    dummy1::Matrix{Float64}
    dummy2::Matrix{Float64}
    dΓ_ph::Matrix{Float64}
end

function generate_2b_kets(dim1b, holes, particles)
    kets = ket2b[ ]
    idx2b = Dict{UInt64, Int64}()
    count = 0
    for i = 1:dim1b
        for j = 1:dim1b
            count += 1
            pp = hh = ph = false
            if i in holes && j in holes
                hh = true
            elseif i in particles && j in particles
                pp = true
            else
                ph = true
            end
            push!(kets, ket2b(i, j, pp, hh, ph))
            tkey = hash_key2(i,j)
            idx2b[tkey] = count
        end
    end
    return idx2b, kets
end

function set_zero_Operator(dim1b, kets)::Operator
    dim2b = size(kets, 1)
    return Operator(zeros(Float64, 1), zeros(Float64, dim1b, dim1b), zeros(dim2b, dim2b))
end

function check_pair_braket(ket::ket2b)
    i = ket.i; j = ket.j
    if i % 2 == 1
        if j == i + 1
            return true
        end
    else
        if j == i - 1
            return true
        end
    end
    return false
end

function set_HNO!(Hamiltonian, F, gval, kets)
    # H1b is to be Fock matrix
    Hamiltonian.onebody .= F

    # H2b is to be interaction matrix
    H2b = Hamiltonian.twobody
    count = 0
    for n = 1:length(kets)
        bra = kets[n]
        if bra.ph; continue; end
        if !check_pair_braket(bra); continue;end
        ordered_bra = ifelse(bra.i < bra.j, 1, 2)
        factor_bra = (-1)^(ordered_bra)
        for m = n:length(kets)
            ket = kets[m]
            if ket.ph; continue; end
            if !check_pair_braket(ket); continue;end
            ordered_ket = ifelse(ket.i < ket.j, 1, 2)
            factor_ket = (-1)^(ordered_ket + 1) # because of Vijkl associated with a^+_i a^+_j a_l a_k
            count += ifelse(m==n, 1, 2)
            H2b[n, m] = H2b[m, n] = factor_bra * factor_ket * gval
        end
    end
    return nothing
end

function get_idx2b(idx2b, i, j)
    tkey = hash_key2(i,j)
    return idx2b[tkey]
end

"""
    eta_white_atan!(eta_1b, eta_2b, f, Γ, idx2b, kets, holes, particles)

Destructive function to evaluate the generator of IMSRG flow equation, η.

```math
\\begin{align}
\\eta & = \\eta^{(1)} + \\eta^{(2)} \\nonumber \\\\
& = \\frac{1}{2} \\sum_{ph} \\arctan\\left(\\frac{f_{ph}(s)}{\\Delta_{ph}(s)}\\right) \\{  a^\\dagger_p a_h \\}
+ \\frac{1}{8} \\sum_{pp'hh'} \\arctan\\left(\\frac{2 \\Gamma_{pp'h'h}(s)}{\\Delta_{pp'hh'}(s)}\\right) \\{  a^\\dagger_p a^\\dagger_{p'} a_h a_{h'} \\} - H.c.. \\nonumber
\\end{align}
```
"""
function eta_white_atan!(eta_1b, eta_2b, f, Γ, idx2b, kets, holes, particles)
    # one-body
    for a in particles
        for i in holes
            idx_ai = get_idx2b(idx2b, a, i)            
            denominator = f[a, a] - f[i, i] + Γ[idx_ai, idx_ai]
            tmp =  0.5 * atan( 2 * f[a,i] / denominator)
            eta_1b[a, i] =  tmp
            eta_1b[i, a] = -tmp
        end
    end
    # two-body
    for a in particles
        for b in particles
            idx_ab = get_idx2b(idx2b, a, b)
            for i in holes
                idx_ai = get_idx2b(idx2b, a, i)
                idx_bi = get_idx2b(idx2b, b, i)
                for j in holes
                    delta = f[a, a] + f[b, b] - f[i, i] - f[j, j] 
                    idx_ij = get_idx2b(idx2b, i, j)
                    idx_aj = get_idx2b(idx2b, a, j)
                    idx_bj = get_idx2b(idx2b, b, j)
                    delta += Γ[idx_ab, idx_ab]
                    delta += Γ[idx_ij, idx_ij]
                    delta -= Γ[idx_ai, idx_ai]
                    delta -= Γ[idx_aj, idx_aj]
                    delta -= Γ[idx_bi, idx_bi]
                    delta -= Γ[idx_bj, idx_bj]
                    tmp = 0.5 * atan( 2 * Γ[idx_ab, idx_ij] / delta)
                    eta_2b[idx_ab, idx_ij] =  tmp
                    eta_2b[idx_ij, idx_ab] = -tmp
                end
            end
        end
    end
    return nothing
end

function show_status(s, Hamiltonian, PT2, eta_1b, eta_2b)
    println(@sprintf("%7.3f  %12.6f  %12.6f  %12.6f  %12.6f", 
            s, 
            Hamiltonian.zerobody[1],
            PT2,
            norm(eta_1b), norm(eta_2b)))
    return nothing
end

"""
Evaluating ``A_{i\\bar{j}k\\bar{l}} = - A_{ijkl}``
"""
function eval_2b_ph_trans!(eta, eta_ph, Γ, Γ_ph, idx2b, kets)
    dim = size(eta, 1)
    for idx_bra = 1:dim
        bra = kets[idx_bra]; i = bra.i; j = bra.j
        for idx_ket = 1:dim
            ket = kets[idx_ket]; k = ket.i; l = ket.j
            org_idx_bra = get_idx2b(idx2b, i, l)
            org_idx_ket = get_idx2b(idx2b, k, j)
            eta_ph[idx_bra, idx_ket] = -eta[org_idx_bra, org_idx_ket]
            Γ_ph[idx_bra, idx_ket] = -Γ[org_idx_bra, org_idx_ket]
        end
    end
    return nothing
end

function eval_2b_ph_trans_inv!(A, A_ph, idx2b, kets)
    dim = length(kets)
    for idx_bra = 1:dim
        bra = kets[idx_bra]; i = bra.i; j = bra.j
        for idx_ket = 1:dim
            ket = kets[idx_ket]; k = ket.i; l = ket.j            
            org_idx_bra = get_idx2b(idx2b, i, l)
            org_idx_ket = get_idx2b(idx2b, k, j)
            A[org_idx_bra, org_idx_ket] = - A_ph[idx_bra, idx_ket]
        end
    end
    return nothing
end

function prepare_idx_kets_ph(kets)
    idx_kets_ph = Int64[]
    for i = 1:length(kets)
        ket = kets[i]
        if ket.ph
            push!(idx_kets_ph, i)
        end
    end
    return idx_kets_ph
end

function eval_derivatives!(dH::Operator, eta_1b, eta_2b, f, Γ, idx2b, coppy_matrices, kets, kets_ph,                     
                           holes, particles, holes_2b, particles_2b, to, debug_mode)
    dim1b = size(f, 1)
    ##################
    ## 0-body dE/ds ##
    ##################
    @timeit to "0-body"  begin
        dE = 0
        # 1st term: this does not contribute to dE in pairing Hamiltonian (diagonal eta/f)
        for i in holes
            for a in particles
                dE += eta_1b[i,a] * f[a,i] - eta_1b[a,i] * f[i, a]
            end
        end
        # 2nd term
        for a in holes
            for b in holes
                for c in particles
                    for d in particles
                        idx_ab = get_idx2b(idx2b, a, b)
                        idx_cd = get_idx2b(idx2b, c, d)
                        dE += 0.5 * eta_2b[idx_ab, idx_cd] * Γ[idx_cd, idx_ab]
                    end
                end
            end
        end
        dH.zerobody[1] = dE
    end

    ##################
    ## 1-body df/ds ##
    ##################
    @timeit to "1-body" begin
        df = dH.onebody
        # 1st term: doesn't contribute
        df .= eta_1b * f - f * eta_1b 
        # 2nd term: doesn't contribute
        for i = 1:dim1b
            for j = 1:dim1b
                tmp = 0.0
                for h in holes
                    idx_hi = get_idx2b(idx2b, h, i)
                    idx_hj = get_idx2b(idx2b, h, j)
                    for p in holes
                        idx_pi = get_idx2b(idx2b, p, i)
                        idx_pj = get_idx2b(idx2b, p, j)                    
                        tmp += eta_1b[h,p] * Γ[idx_pi, idx_hj]
                        tmp -= f[h,p] * eta_2b[idx_pi, idx_hj]
                        tmp -= eta_1b[p,h] * Γ[idx_hi, idx_pj]
                        tmp += f[p,h] * eta_2b[idx_hi, idx_pj]
                    end
                end
                df[i,j] += tmp
            end
        end
        # 3rd term: contributing
        for i = 1:dim1b
            for j =1:dim1b
                tmp = 0.0
                for c in holes
                    idx_ci = get_idx2b(idx2b,c,i)
                    idx_cj = get_idx2b(idx2b,c,j)
                    eta_2b_ci_pp = @view eta_2b[idx_ci, particles_2b]
                    eta_2b_cj_pp = @view eta_2b[idx_cj, particles_2b]
                    Γ_pp_cj = @view Γ[particles_2b, idx_cj]
                    Γ_pp_ci = @view Γ[particles_2b, idx_ci]
                    tmp += dot(eta_2b_ci_pp, Γ_pp_cj)
                    tmp += dot(eta_2b_cj_pp, Γ_pp_ci)
                end
                for c in particles
                    idx_ci = get_idx2b(idx2b,c,i)
                    idx_cj = get_idx2b(idx2b,c,j)
                    eta_2b_ci_hh = @view eta_2b[idx_ci, holes_2b]
                    eta_2b_cj_hh = @view eta_2b[idx_cj, holes_2b]
                    Γ_hh_cj = @view Γ[holes_2b, idx_cj]
                    Γ_hh_ci = @view Γ[holes_2b, idx_ci]
                    tmp += dot(eta_2b_ci_hh, Γ_hh_cj)
                    tmp += dot(eta_2b_cj_hh, Γ_hh_ci)
                end 
                df[i,j] += 0.5 * tmp
            end
        end
    end

    ##################
    ## 2-body dΓ/ds ##
    ################## 
    @timeit to "2-body" begin
        dΓ = dH.twobody; dΓ .= 0.0
        @timeit to "122" @threads for i_pq = 1:dim1b^2
            p = div(i_pq-1, dim1b) + 1
            q = mod(i_pq-1, dim1b) + 1
            idx_bra = get_idx2b(idx2b, p, q)
            for r = 1:dim1b
                for s = 1:dim1b
                    idx_ket = get_idx2b(idx2b, r, s)
                    tmp = 0.0
                    for a = 1:dim1b
                        idx_aq = get_idx2b(idx2b, a, q)
                        idx_ap = get_idx2b(idx2b, a, p)
                        idx_as = get_idx2b(idx2b, a, s)
                        idx_ar = get_idx2b(idx2b, a, r)

                        tmp += eta_1b[p,a] * Γ[idx_aq, idx_ket]
                        tmp -= f[p,a] * eta_2b[idx_aq, idx_ket]
                        tmp -= eta_1b[q,a] * Γ[idx_ap, idx_ket]
                        tmp += f[q,a] * eta_2b[idx_ap, idx_ket]

                        tmp -= eta_1b[a,r] * Γ[idx_bra, idx_as]
                        tmp += f[a,r] * eta_2b[idx_bra, idx_as]
                        tmp += eta_1b[a,s] * Γ[idx_bra, idx_ar]
                        tmp -= f[a,s] * eta_2b[idx_bra, idx_ar]

                    end
                    dΓ[idx_bra,idx_ket] += tmp
                end
            end
        end
        if debug_mode > 2
            println("norm(dΓ) @1  ", norm(dΓ), " ", norm(dΓ,1))
        end
        dummy1 = coppy_matrices.dummy1
        Mat_Npp = coppy_matrices.Mat_Npp
        Mat_Nhh = coppy_matrices.Mat_Nhh
        @timeit to "222pphh" begin  
            BLAS.gemm!('N', 'N', 0.5, Mat_Npp, Γ, 0.0, dummy1); BLAS.gemm!('N', 'N', 1.0, eta_2b, dummy1, 1.0, dΓ)
            BLAS.gemm!('N', 'N',-0.5, Γ, Mat_Npp, 0.0, dummy1); BLAS.gemm!('N', 'N', 1.0, dummy1, eta_2b, 1.0, dΓ)
            BLAS.gemm!('N', 'N',-0.5, Mat_Nhh, Γ, 0.0, dummy1); BLAS.gemm!('N', 'N', 1.0, eta_2b, dummy1, 1.0, dΓ)
            BLAS.gemm!('N', 'N', 0.5, Γ, Mat_Nhh, 0.0, dummy1); BLAS.gemm!('N', 'N', 1.0, dummy1, eta_2b, 1.0, dΓ)     
            if debug_mode > 2
                println("norm(dΓ) @2  ", norm(dΓ), " ", norm(dΓ,1))
            end
        end
        # 3rd term: this part is usually the most expensive.
        # The following code (commented out) is naive implementation, and the better way below
        # is to replace loop by matrix multiplication via particle-hole transformation

        # @timeit to "222ph" @threads for p = 1:dim1b
        #     for q = 1:dim1b
        #         idx_bra = get_idx2b(idx2b, p, q)
        #         for r = 1:dim1b
        #             for s = 1:dim1b
        #                 idx_ket = get_idx2b(idx2b, r, s)
        #                 tmp = 0.0
        #                 for hp in holes
        #                     for pp in particles
        #                         idx_pp_p = get_idx2b(idx2b, pp, p); idx_pp_q = get_idx2b(idx2b, pp, q)
        #                         idx_pp_r = get_idx2b(idx2b, pp, r); idx_pp_s = get_idx2b(idx2b, pp, s)
        #                         idx_hp_p = get_idx2b(idx2b, hp, p); idx_hp_q = get_idx2b(idx2b, hp, q)
        #                         idx_hp_r = get_idx2b(idx2b, hp, r); idx_hp_s = get_idx2b(idx2b, hp, s)

        #                         tmp -= eta_2b[idx_pp_q, idx_hp_s] * Γ[idx_hp_p, idx_pp_r]
        #                         tmp += eta_2b[idx_pp_p, idx_hp_s] * Γ[idx_hp_q, idx_pp_r]
        #                         tmp += eta_2b[idx_pp_q, idx_hp_r] * Γ[idx_hp_p, idx_pp_s]
        #                         tmp -= eta_2b[idx_pp_p, idx_hp_r] * Γ[idx_hp_q, idx_pp_s]

        #                         tmp += eta_2b[idx_hp_q, idx_pp_s] * Γ[idx_pp_p, idx_hp_r]
        #                         tmp -= eta_2b[idx_hp_p, idx_pp_s] * Γ[idx_pp_q, idx_hp_r]
        #                         tmp -= eta_2b[idx_hp_q, idx_pp_r] * Γ[idx_pp_p, idx_hp_s]
        #                         tmp += eta_2b[idx_hp_p, idx_pp_r] * Γ[idx_pp_q, idx_hp_s]
        #                     end
        #                 end
        #                 dΓ[idx_bra, idx_ket] += tmp    
        #             end
        #         end
        #     end
        # end
    
        @timeit to "222ph w/ ph trans." begin
            dim = size(dΓ, 1)
            Mat_Nab = coppy_matrices.Mat_Nab
            eta_ph = coppy_matrices.eta_ph
            Γ_ph = coppy_matrices.Γ_ph
            dΓ_ph = coppy_matrices.dΓ_ph
            dummy1 = coppy_matrices.dummy1
            dummy2 = coppy_matrices.dummy2
            BLAS.gemm!('N', 'N', 1.0, Mat_Nab, Γ_ph, 0.0, dummy1)
            BLAS.gemm!('N', 'N', 1.0, eta_ph, dummy1, 0.0, dummy2)
            eval_2b_ph_trans_inv!(dΓ_ph, dummy2, idx2b, kets)

            # Completing permutation to take ph transformation
            dummy1 .= dΓ_ph
            for ibra = 1:dim
                bra = kets[ibra]; i = bra.i; j = bra.j
                idx_ji = get_idx2b(idx2b, j, i)
                for iket = 1:dim
                    ket = kets[iket]; k = ket.i; l = ket.j
                    idx_lk = get_idx2b(idx2b, l, k)
                    dummy1[ibra,iket] -= dΓ_ph[idx_ji, iket]
                    dummy1[ibra,iket] -= dΓ_ph[ibra, idx_lk]
                    dummy1[ibra,iket] += dΓ_ph[idx_ji, idx_lk]
                end
            end
            dΓ_ph .= dummy1
            dΓ .+= dΓ_ph

        end
        if debug_mode > 2
            println("norm(dΓ) @END ", norm(dΓ), " ", norm(dΓ,1))
        end
    end
    return nothing
end

function get_2b_holeparticle_idx(kets, holes, particles)
    holes_2b = Int64[]
    particles_2b = Int64[]
    for i = 1:length(kets)
        ket = kets[i]
        if ket.hh
            push!(holes_2b, i)
        elseif ket.pp
            push!(particles_2b, i)
        end
    end
    return holes_2b, particles_2b
end

function flatten_operator!(Op::Operator, y::Vector{Float64})
    y[1:length(Op.zerobody)] .= Op.zerobody
    y[length(Op.zerobody)+1:length(Op.zerobody)+length(Op.onebody)] .= Op.onebody[:]
    y[length(Op.zerobody)+length(Op.onebody)+1:end] .= Op.twobody[:]
    return nothing
end

function update_ys!(ys, dyds, ds)
    axpy!(ds, dyds, ys)
    return nothing
end

function update_Op_from_fvec!(Hamiltonian, ys)
    dim0b = length(Hamiltonian.zerobody)
    dim1b = size(Hamiltonian.onebody, 1)
    dim2b = size(Hamiltonian.twobody, 1)
    Hamiltonian.zerobody .= ys[1:dim0b]
    Hamiltonian.onebody .= reshape(ys[dim0b+1:dim0b+dim1b*dim1b], dim1b, dim1b)
    Hamiltonian.twobody .= reshape(ys[dim0b+dim1b*dim1b+1:end], dim2b, dim2b)
    return nothing
end

function PT2_from_HNO(Hamiltonian, holes, particles, idx2b)
    E_PT2 = 0.0
    F = Hamiltonian.onebody
    Γ = Hamiltonian.twobody
    for i in holes
        for j in holes
            for a in particles
                for b in particles
                    Delta = F[i,i] + F[j,j] - F[a,a] - F[b,b] 
                    ME = Γ[get_idx2b(idx2b, a, b), get_idx2b(idx2b, i, j)]
                    E_PT2 += 0.5 * ME^2 / Delta
                end
            end
        end
    end
    return E_PT2
end

function making_matrices_for_222ph(Γ, kets, holes)
    Mat_Nab = zeros(Float64,  size(Γ))
    Mat_Npp = zeros(Float64,  size(Γ))
    Mat_Nhh = zeros(Float64,  size(Γ))
    for i = 1:size(Γ,1)
        ket = kets[i]
        b = ket.i; a = ket.j
        na = ifelse(a in holes, 1.0, 0.0)
        nb = ifelse(b in holes, 1.0, 0.0)
        Mat_Nab[i,i] = (na - nb)
        Mat_Npp[i,i] = ifelse(ket.pp, 1.0, 0.0)
        Mat_Nhh[i,i] = ifelse(ket.hh, 1.0, 0.0)
    end
    eta_ph = zeros(Float64, size(Γ))
    Γ_ph = zeros(Float64, size(Γ))
    dummy1 = zeros(Float64, size(Γ))
    dummy2 = zeros(Float64, size(Γ))
    dΓ_ph = zeros(Float64, size(Γ))
    return Mat222ph(Mat_Nab, Mat_Npp, Mat_Nhh, eta_ph, Γ_ph, dummy1, dummy2, dΓ_ph)
end

function _main_IMSRG(HNO, holes, particles, gval, Nocc, to, debug_mode;
                    smax = 15.0, tol_eta=1.e-6, ds=1.e-2)
    F = HNO.f; dim1b = size(F,1)
    # Define 2b-kets for the model space
    idx2b, kets = generate_2b_kets(dim1b, holes, particles)
    if debug_mode > 0
        println("# of 2b kets = ", size(kets, 1))
    end

    # Define the Hamiltonian operator
    Hamiltonian = set_zero_Operator(dim1b, kets)
    Hamiltonian.zerobody .= HNO.E
    set_HNO!(Hamiltonian, F, gval, kets)
    if debug_mode > 2
        show_matrix( "F(H1b)", Hamiltonian.onebody)
    end
    Γ = Hamiltonian.twobody
    dim_fvec = length(Hamiltonian.zerobody) + length(Hamiltonian.onebody) + length(Γ) 
    ys = zeros(Float64, dim_fvec)
    flatten_operator!(Hamiltonian, ys)

    # Define the eta, generator of IMSRG-flow 
    eta_update_func = eta_white_atan!
    dim1b = size(F, 1);   eta_1b = zeros(Float64, dim1b, dim1b)
    dim2b = size(kets, 1);eta_2b = zeros(Float64, dim2b, dim2b)
    eta_update_func(eta_1b, eta_2b, F, Γ, idx2b, kets, holes, particles)
    idx_kets_ph = prepare_idx_kets_ph(kets)
    kets_ph = @view kets[idx_kets_ph]
    coppy_matrices = making_matrices_for_222ph(Γ, kets, holes)
    eval_2b_ph_trans!(eta_2b, coppy_matrices.eta_ph, Γ, coppy_matrices.Γ_ph, idx2b, kets)

    # Prepare derivatives
    dH = deepcopy(Hamiltonian)
    holes_2b, particles_2b = get_2b_holeparticle_idx(kets, holes, particles)   
    dyds = zeros(Float64, dim_fvec)
    flatten_operator!(dH, dyds)

    # Start the IMSRG flow
    s = 0.0
    if debug_mode > 0
        println("Dim. of flattened operator = ", length(dyds))
        println(@sprintf("%7s  %12s  %12s  %12s  %12s", "s", "Energy", "MBPT2", "||eta_1b||", "||eta_2b||"))
        show_status(s, Hamiltonian, NaN, eta_1b, eta_2b)
    end
    @timeit to "flow" while s < smax && norm(eta_1b) + norm(eta_2b) > tol_eta
        s += ds
        # Evaluate derivatives R.H.S. of the flow equation => updating flatten vecotor dyds
        @timeit to "eval" begin
            eval_derivatives!(dH, eta_1b, eta_2b, F, Γ, idx2b, coppy_matrices, kets, kets_ph,
                             holes, particles, holes_2b, particles_2b, to, debug_mode)
            flatten_operator!(dH, dyds) 
        end

        # y(s+ds) = y(s) + dy/ds * ds + ... by a solver (Euler method for now)
        update_ys!(ys, dyds, ds) 

        # H(s+ds) <= y(s+ds)
        update_Op_from_fvec!(Hamiltonian, ys) # update Hamiltonian from ys

        # Update eta(s+ds) by the new Hamiltonian
        F = Hamiltonian.onebody; Γ = Hamiltonian.twobody
        eta_update_func(eta_1b, eta_2b, F, Γ, idx2b, kets, holes, particles)    
        
        # developing
        eval_2b_ph_trans!(eta_2b, coppy_matrices.eta_ph, Γ, coppy_matrices.Γ_ph, idx2b, kets)

        PT2 = PT2_from_HNO(Hamiltonian, holes, particles, idx2b)
        if debug_mode > 0
            show_status(s, Hamiltonian, PT2, eta_1b, eta_2b)
        end
        if abs(PT2) < tol_eta
            break
        end

    end    
    return Hamiltonian.zerobody[1]
end
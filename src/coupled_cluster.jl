struct chi_intermediate
    χ_hh::Array{Float64, 2}
    χ_pp::Array{Float64, 2}
    χ_hhhh::Array{Float64, 2}
    χ_pppp::Array{Float64, 2}
    χ_hphp::Array{Float64, 2}
end

"""
    hash_key4(i, j, a, b)::UInt64

Function to encode `i,j,a,b` into a single UInt64.
From the lowest 8 bits, `b`, `a`, `j`, `i` are stored.
If you want to work with larger model space, requiring e.g. single particle states > 255, you need to modify this function.
"""
function hash_key4(i, j, a, b)::UInt64
    return (UInt64(i) << 24) + (UInt64(j) << 16) + (UInt64(a) << 8) + UInt64(b)
end

"""
    unhash_key4(key::UInt64)::Tuple{Int, Int, Int, Int}

Function to decode `i,j,a,b` from a single UInt64.
"""
function unhash_key4(key::UInt64)
    i = Int(key >> 24)
    j = Int((key >> 16) & 0xFF)
    a = Int((key >> 8) & 0xFF)
    b = Int(key & 0xFF)
    return i, j, a, b
end

"""
    CCD(F, gval, Nocc, to, debug_mode; itnum=100, tol=1e-9, mix_ratio=0.25)

Main function to calculate the ground state energy of the pairing Hamiltonian using the coupled cluster method.
Since the Hamiltonian is the pairing model, we don't need to consider 1p1h excitations, CCD instead of CCSD, etc.

# Arguments
- `F::Array{Float64, 2}`: Fock matrix from the Hartree-Fock calculation
- `gval::Float64`: strength of the pairing interaction
- `Nocc::Int`: number of occupied states
- `to::TimerOutput`: timer object to measure the elapsed time
- `debug_mode::Int`: specifying the debug mode

# Optional arguments
- `itnum::Int(100)`: maximum number of iterations
- `tol::Float64(1e-9)`: convergence criterion
- `mix_ratio::Float64(0.25)`: mixing ratio for the update of the t2 amplitudes

# References
- [Shavitt, I., & Bartlett, R. J. (2009). Many-body methods in chemistry and physics: MBPT and coupled-cluster theory. Cambridge University Press.](https://www.cambridge.org/jp/academic/subjects/chemistry/physical-chemistry/many-body-methods-chemistry-and-physics-mbpt-and-coupled-cluster-theory)
- [Justin G. Lietz, Samuel Novario, Gustav R. Jansen, Gaute Hagen & Morten Hjorth-Jensen, Computational Nuclear Physics and Post Hartree-Fock Methods](https://link.springer.com/chapter/10.1007/978-3-319-53336-0_8), part of the book series: [Lecture Notes in Physics ((LNP,volume 936))](https://link.springer.com/book/10.1007/978-3-319-53336-0)
"""
function CCD(F, gval, Nocc, to, debug_mode; itnum=100, tol=1e-9, mix_ratio=0.25)
    dim = size(F, 1)
    dim_p = dim - Nocc
    dim_h = Nocc

    dim_hh, dim_hp, dim_pp = eval_2bket_dim(dim_h, dim_p)
    t2_pphh = zeros(Float64, dim_pp, dim_hh)
    t2_pphh_new =  zeros(Float64, dim_pp, dim_hh)
    Hbar_pphh = zeros(Float64, dim_pp, dim_hh)

    # Initialize CC coeffs by MBPT2
    for a = 1:2:dim_p
        b = a + 1 
        idx_pp = get_vecidx(a, b, dim_p)
        for i = 1:2:Nocc
            j = i+1
            idx_hh = get_vecidx(i, j, dim_h)
            t2b = -gval / (F[i,i] + F[j,j] - F[a+dim_h,a+dim_h] - F[b+dim_h,b+dim_h])
            t2_pphh[idx_pp, idx_hh] = t2b 
        end
    end
    E_CC = eval_ECC(F, t2_pphh, Nocc, gval)

    χ = construct_chi(F,t2_pphh, Nocc, gval, debug_mode)
    for it = 1:itnum
        update_χ!(χ, F, t2_pphh, Nocc, gval)
        eval_Hbar_pphh!(Hbar_pphh, t2_pphh, t2_pphh_new, F, χ, Nocc, gval)
        E = eval_ECC(F, t2_pphh_new, Nocc, gval)
        if debug_mode > 0
            println("It $it, E = $E norm t2 $(2*norm(t2_pphh)) new $(2*norm(t2_pphh_new))")
        end
        t2_pphh .= mix_ratio * t2_pphh_new .+ (1-mix_ratio) * t2_pphh
        if abs(E-E_CC) < tol
            break
        end 
        E_CC = E
    end
    return E_CC
end

function eval_2bket_dim(Nh, Np)
    dim_hh = div(Nh*(Nh-1), 2)
    dim_hp = Np * Nh
    dim_pp = div(Np*(Np-1), 2)
    return dim_hh, dim_hp, dim_pp
end

function eval_χ_hh!(χ_hh, F, t2_pphh, Nocc, gval)
    dim_h = size(χ_hh, 1)
    dim_p = size(F, 1) - Nocc
    for k = 1:2:Nocc
        j = k
        l = k + 1
        t2b = F[k, j]
        for c = 1:2:dim_p
            d = c + 1
            idx_hh = get_vecidx(j, l, dim_h)
            idx_pp = get_vecidx(c, d, dim_p)
            t2b += t2_pphh[idx_pp, idx_hh] * (-gval)
        end
        χ_hh[k, j] = t2b
        χ_hh[k+1, j+1] = t2b
    end
    return nothing
end

function eval_χ_pp!(χ_pp, F, t2_pphh, Nocc, gval)
    dim_h = Nocc
    dim_p = size(F,1) - Nocc 
    for b = 1:2:dim_p
        c = b 
        d = c + 1
        tmp = F[b+dim_h, b+dim_h] 
        for k = 1:2:dim_h
            l = k + 1
            idx_hh = get_vecidx(k, l, dim_h)
            idx_pp = get_vecidx(b, d, dim_p)
            tmp += - t2_pphh[idx_pp, idx_hh] * (-gval)
        end
        χ_pp[b, c] = tmp
        χ_pp[b+1, c+1] = tmp
    end
    return nothing
end

"""
To get continuous idx for pp or hh from `i` and `j`.
"""
function get_vecidx(i, j, N)
    @assert i != j "i == j should not happen in get_vecidx"
    if i > j
        i, j = j, i
    end
    return N*(i-1) - div(i*(i-1), 2) + j - i
end

"""
To get `i` and `j` for pp or hh from a vectorized index.
Note that `i` must be less than `j`.
"""
function vecidx_to_idx(idx, N)
    if idx <= N-1
        return 1, idx+1
    end
    i = 1 
    ksum = N - 1
    while ksum < idx
        i += 1
        ksum += N - i
    end
    j = idx + i - ( ksum  - (N-i) )
    return i, j
end

"""
```math
\\langle kl | \\chi | ij \\rangle  = \\langle kl | v | ij \\rangle + \\frac{1}{2} \\sum_{cd} \\langle kl | v | cd \\rangle \\langle cd | t | ij \\rangle
```
"""
function eval_χ_hhhh!(χ_hhhh, t2_pphh, dim, Nocc, gval) # 8.48
    dim_p = dim - Nocc
    dim_h = Nocc
    dim_hh, dim_hp, dim_pp = eval_2bket_dim(dim_h, dim_p)
    for idx_bra = 1:dim_hh
        k, l = vecidx_to_idx(idx_bra, dim_h)
        is_pair = (k%2==1 && l == k+1)
        if !is_pair
            continue
        end
        for idx_ket = 1:dim_hh
            i, j = vecidx_to_idx(idx_ket, dim_h)
            is_pair = (i%2==1 && j == i+1)
            if !is_pair
                continue
            end            
            t2b = -gval 
            for c = 1:2:dim_p
                d = c + 1
                idx_cd = get_vecidx(c, d, dim_p)
                t2b += t2_pphh[idx_cd, idx_ket] * (-gval)
            end
            χ_hhhh[idx_bra, idx_ket] = t2b
        end
    end
    return nothing
end

function eval_χ_pppp!(χ_pppp, t2_pphh, dim, Nocc, gval) # 8.50
    dim_p = dim - Nocc
    dim_h = Nocc
    dim_hh, dim_hp, dim_pp = eval_2bket_dim(dim_h,dim_p)
    for idx_bra = 1:dim_pp
        a,b = vecidx_to_idx(idx_bra, dim_p)
        is_pair = (a%2==1 && b == a+1)
        if !is_pair
            continue
        end
        for idx_ket = idx_bra:dim_pp
            c, d = vecidx_to_idx(idx_ket, dim_p)
            is_pair = (c%2==1 && d == c+1)
            if !is_pair
                continue
            end            
            χ_pppp[idx_bra, idx_ket] = χ_pppp[idx_ket, idx_bra] = -gval
        end
    end
    return nothing
end

function eval_χ_hphp!(χ_hphp, t2_pphh, dim, Nocc, gval) # 8.49
    dim_p = dim - Nocc
    dim_h = Nocc
    for k = 1:dim_h
        j = k 
        l = ifelse(k%2==1, k+1, k-1)
        for c = 1:dim_p
            b = c
            d = ifelse(c%2==1, c+1, c-1)
            idx_bra = (k-1)*dim_p + b
            idx_ket = (j-1)*dim_p + c          
            idx_hh = get_vecidx(k, l, dim_h)
            idx_pp = get_vecidx(c, d, dim_p)
            χ_hphp[idx_bra,idx_ket] = 0.5 * t2_pphh[idx_pp,idx_hh] * (-gval)
        end
    end
    return nothing
end

"""
`i,j,k` are hole indices, `a,b,c` are particle indices.
"""
function construct_chi(F, t2_pphh, Nocc, gval, debug_mode=0)
    dim = size(F, 1)
    dim_p = dim - Nocc
    dim_h = Nocc
    χ_hh = zeros(Float64, dim_h, dim_h)
    χ_pp = zeros(Float64, dim_p, dim_p)
    eval_χ_hh!(χ_hh, F, t2_pphh, Nocc, gval)
    eval_χ_pp!(χ_pp, F, t2_pphh, Nocc, gval)
    if debug_mode > 1
        show_matrix("χpp", χ_pp)
    end

    dim_hh, dim_hp, dim_pp = eval_2bket_dim(dim_h, dim_p)
    χ_hhhh = zeros(Float64, dim_hh, dim_hh)
    χ_pppp = zeros(Float64, dim_pp, dim_pp)
    χ_hphp = zeros(Float64, dim_hp, dim_hp)

    # evaluate χ_hhhh
    eval_χ_hhhh!(χ_hhhh, t2_pphh, dim, Nocc, gval)
    if debug_mode > 1
        show_matrix("χhhhh", χ_hhhh)
    end

    # evaluate χ_pppp
    eval_χ_pppp!(χ_pppp, t2_pphh, dim, Nocc, gval)
    if debug_mode > 1
        show_matrix("χpppp", χ_pppp)
    end

    # evaluate χ_hphp
    eval_χ_hphp!(χ_hphp, t2_pphh, dim, Nocc, gval)
    if debug_mode > 1
        show_matrix("χhphp", χ_hphp)
    end

    return chi_intermediate(χ_hh, χ_pp, χ_hhhh, χ_pppp, χ_hphp)
end

function update_χ!(χ, F, t2_pphh, Nocc, gval)
    dim = size(F, 1)
    eval_χ_hh!(χ.χ_hh, F, t2_pphh, Nocc, gval)
    eval_χ_pp!(χ.χ_pp, F, t2_pphh, Nocc, gval)
    eval_χ_hhhh!(χ.χ_hhhh, t2_pphh, dim, Nocc, gval)
    eval_χ_pppp!(χ.χ_pppp, t2_pphh, dim, Nocc, gval)
    eval_χ_hphp!(χ.χ_hphp, t2_pphh, dim, Nocc, gval)
    return nothing
end

function eval_Hbar_pphh!(Hbar_pphh, t2_pphh, t2_pphh_new, F, χ::chi_intermediate, Nocc, gval)
    dim = size(F, 1)
    dim_p = dim - Nocc
    dim_h = Nocc
    Hbar_pphh .= 0.0
    t2_pphh_new .=0.0
    for i = 1:2:dim_h
        j = i + 1
        idx_ij = get_vecidx(i, j, dim_h)
        for a = 1:2:dim_p
            b = a + 1
            idx_ab = get_vecidx(a, b, dim_p)
            tmp =  -gval 

            ## chi_p 
            χp = χ.χ_pp
            c = a + 1
            idx_ac = get_vecidx(a, c, dim_p)
            tmp += χp[b, c] * t2_pphh[idx_ac, idx_ij] # 2 come from Pab

            # a <-> b, then c = var{b} = a 
            c = a
            idx_bc = get_vecidx(b, c, dim_p)
            tmp += χp[a, c] * t2_pphh[idx_bc, idx_ij] # sign cancels due to the b > c 

            ## chi_h
            χh = χ.χ_hh
            k = i + 1
            idx_ik = get_vecidx(i, k, dim_h)
            tmp += - χh[k, j] * t2_pphh[idx_ab, idx_ik]

            k = i # j-1 Pij
            idx_jk = get_vecidx(j, k, dim_h)
            tmp += χh[k, i] * t2_pphh[idx_ab, idx_jk] * (-1.0) # sign due to the j > k
            
            #chi_hphp            
            χhphp = χ.χ_hphp
            for k = 1:dim_h                
                for c = 1:dim_p
                    idx_ka = (k-1)*dim_p + a
                    idx_kb = (k-1)*dim_p + b
                    idx_ci = (i-1)*dim_p + c
                    idx_cj = (j-1)*dim_p + c
                    if a != c 
                        idx_ac = get_vecidx(a, c, dim_p)
                        if i != k                            
                            tmp += χhphp[idx_kb, idx_cj] * t2_pphh[idx_ac,idx_ik]
                        end
                        # i <-> j
                        if j != k
                            idx_jk = get_vecidx(j, k, dim_h)
                            tmp += - χhphp[idx_kb, idx_ci] * t2_pphh[idx_ac,idx_jk] * ifelse(j>k,-1,1) * ifelse(a>c,-1,1)
                        end
                    end
                    # a <-> b
                    if b != c
                        idx_bc = get_vecidx(b, c, dim_p)
                        if i != k
                            idx_ik = get_vecidx(i, k, dim_h)
                            tmp += - χhphp[idx_ka, idx_cj] * t2_pphh[idx_bc,idx_ik] * ifelse(b>c,-1,1) * ifelse(i>k,-1,1)
                        end                       
                    end
                    # i <-> j, a <-> b
                    if j != k && b != c
                        idx_jk = get_vecidx(j, k, dim_h)
                        idx_bc = get_vecidx(b, c, dim_p)
                        idx_ka = (k-1)*dim_p + a
                        idx_ci = (i-1)*dim_p + c
                        tmp += χhphp[idx_ka, idx_ci] * t2_pphh[idx_bc,idx_jk] * ifelse(j>k,-1,1) * ifelse(b>c,-1,1)
                    end
                end
            end            
            Hbar_pphh[idx_ab, idx_ij] = tmp 
        end
    end
    BLAS.gemm!('N', 'N', 1.0, χ.χ_pppp, t2_pphh, 1.0, Hbar_pphh)
    BLAS.gemm!('N', 'N', 1.0, t2_pphh, χ.χ_hhhh, 1.0, Hbar_pphh)
  
    for i = 1:2:dim_h
        j = i + 1
        idx_ij = get_vecidx(i, j, dim_h)
        for a = 1:2:dim_p
            b = a + 1
            idx_ab = get_vecidx(a, b, dim_p) 
            t2_pphh_new[idx_ab, idx_ij] = t2_pphh[idx_ab, idx_ij] + Hbar_pphh[idx_ab, idx_ij] / (F[i,i] + F[j,j] - F[a+dim_h,a+dim_h] - F[b+dim_h,b+dim_h])
        end
    end
    return nothing
end

function eval_ECC(F, t2_pphh, Nocc, gval)
    dim = size(F, 1)
    E2b = 0.0
    for i = 1:2:Nocc
        j = i + 1 
        idx_hh = get_vecidx(i, j, Nocc)
        for a = 1:2:dim-Nocc
            b = a + 1
            idx_pp = get_vecidx(a, b, dim-Nocc)
            E2b += (-gval) * t2_pphh[idx_pp, idx_hh]
        end
    end
    return E2b
end
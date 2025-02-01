struct Hamiltonian_NormalOrdered
    E::Vector{Float64}
    f::Array{Float64,2}
    Gamma::Vector{Float64}
end

function eval_EHF(h1b, h2b, gval, rho, Nocc)
    Nq = size(h1b, 1)
    e1b = e2b = 0.0
    for α = 1:Nocc
        γ = ifelse(α%2==0, α-1, α+1)
        for β = 1:Nocc
            e1b += rho[α, β] * h1b[α, β]
            δ = ifelse(β%2==0, β-1, β+1)
            e2b += 0.5 * (-gval) * rho[α, β] * rho[γ, δ] 
        end        
    end
    EHF = e1b + e2b
    return e1b, e2b, EHF
end

function eval_Fock!(F, rho, h1b, h2b, Nocc, gval)
    Nq = size(h1b, 1)
    F .= 0.0
    for α = 1:Nq
        F[α, α] = h1b[α, α]
    end
    for α = 1:Nq
        γ = ifelse(α%2==0, α-1, α+1)
        if γ > Nocc; continue; end
        for β = α:Nq
            δ = ifelse(β%2==0, β-1, β+1)
            if δ > Nocc; continue; end
            F[α, β] += - gval * rho[γ, δ]
            F[β, α] = F[α, β]
        end
    end
    return nothing
end

function eval_rho!(rho, U, Nocc)
    tU = @view U[:,1:Nocc]
    BLAS.gemm!('N', 'T', 1.0, tU, tU, 0.0, rho)
    return nothing
end

function check_rho(rho, Nocc)
    @assert tr(rho) ≈ Nocc
    return nothing
end

"""
Function to return the indices of the holes, which are below the Fermi level, and the particles, which are above the Fermi level.
"""
function define_holes_particles(F, Nocc)
    idxs = sortperm(diag(F))
    holes = idxs[1:Nocc]
    particles = idxs[Nocc+1:end]
    return holes, particles
end

"""
To get normal ordered Hamiltonian upto two-body level.
"""
function normal_ordering(h1b, h2b, holes, particles, gval)
    Nq = size(h1b, 1)

    # 0-body 
    E0 = 0.0
    for i in holes
        E0 += h1b[i,i]
    end
    idx = 0
    for i = 1:Nq
        for j =i+1:Nq
            idx += 1
            if i in particles; continue; end
            if j in particles; continue; end
            if (i % 2 == 1 && j == i + 1) || (i % 2 == 0 && j == i - 1)
                E0 += h2b[idx]
            end
        end
    end
    E = [E0]

    # 1-body <ih|f|jh> => δ_{ij} (-g)
    f = copy(h1b)
    for h in holes
        hbar = ifelse(h%2==0, h-1, h+1)
        i = j = hbar
        f[i,j] += (-gval)
    end

    # 2-body
    Gamma = copy(h2b)
    return Hamiltonian_NormalOrdered(E,f,Gamma)
end

"""
Compute the energies of the Hamiltonian by the Hartree-Fock method and perturbation theory.
"""
function get_Egs_HF(Nocc, h1b, h2b, gval, to, debug_mode=0; return_only_E0=true, itnum_max=30, tol=1e-9)
    Nq = size(h1b, 1)
    F = zeros(Float64, Nq, Nq)
    rho = zeros(Float64, Nq, Nq)
    U = zeros(Float64, Nq, Nq)
    for i in 1:Nocc # sample(1:Nq, Nocc, replace=false)
        U[i,i] = 1.0
    end    

    Eprev = EHF = 1.e5
    @timeit to "HF" for i = 1:itnum_max
        eval_rho!(rho, U, Nocc)
        eval_Fock!(F, rho, h1b, h2b, Nocc, gval)
        @timeit to "solver: eigen Fock." evals, evecs = eigen(F)
        U .= evecs'
        e1b, e2b, EHF = eval_EHF(h1b, h2b, gval, rho, Nocc)
        if debug_mode > 0
            println("it = $i \t EHF $EHF = ($e1b + $e2b)")
        end
        if abs(EHF - Eprev) < tol
            Eprev = EHF
            break
        end
        Eprev = EHF
    end
    EPT2 = PT2(F, Nocc, gval)
    EPT3 = PT3(F, Nocc, gval)

    holes, particles = define_holes_particles(F, Nocc)
    HNO = normal_ordering(h1b, h2b, holes, particles, gval)
    return EHF, EPT2, EPT3, HNO, holes, particles
end
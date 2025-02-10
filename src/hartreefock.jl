struct Hamiltonian_NormalOrdered
    E::Vector{Float64}
    f::Array{Float64,2}
    Gamma::Vector{Float64}
end

"""
    eval_EHF(h1b, gval, rho, Nocc)

Function to evaluate the energy of the Hamiltonian in the Hartree-Fock method.

```math
\\begin{align}
E &= E_{1b} + E_{2b} \\nonumber \\\\
E_{1b} & = \\sum_{i \\leq F} \\sum_{\\alpha\\beta} C^T_{i \\alpha} C_{i \\beta} \\langle \\alpha | h | \\beta \\rangle =
\\sum_{\\alpha \\beta} \\rho_{\\alpha \\beta} \\langle \\alpha | h | \\beta \\rangle\\nonumber  \\\\
E_{2b} &= \\frac{1}{2}\\sum_{\\alpha\\beta\\gamma\\delta} \\rho_{\\alpha \\beta} \\rho_{\\gamma \\delta} \\langle \\alpha \\gamma | v | \\beta \\delta \\rangle \\nonumber 
\\end{align}
```
"""
function eval_EHF(h1b, gval, rho, Nocc)
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

"""
    eval_Fock!(F, rho, h1b, Nocc, gval)

Function to evaluate the Fock matrix destructively.

```math
F_{\\alpha\\beta} = \\langle \\alpha | h | \\beta \\rangle + \\sum_{\\gamma \\delta} \\rho_{\\gamma \\delta} \\langle \\alpha \\gamma | v | \\beta \\delta \\rangle
```
"""
function eval_Fock!(F, rho, h1b, Nocc, gval)
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

"""
    eval_rho!(rho, U, Nocc)

Update the density matrix with the coefficient matrix.

```math
\\rho_{\\alpha\\beta} = \\sum_{i \\leq F} C^T_{i \\alpha} C_{i \\beta}
```
"""
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
    define_holes_particles(F, Nocc)

Function to return the indices of the holes, which are below the Fermi level, and the particles, which are above the Fermi level.
"""
function define_holes_particles(F, Nocc)
    idxs = sortperm(diag(F))
    holes = idxs[1:Nocc]
    particles = idxs[Nocc+1:end]
    return holes, particles
end

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
    HF(Nocc, h1b, h2b, gval, to, debug_mode=0; itnum_max=100, tol=1e-9)

Main function to compute the energies of the Hamiltonian by the Hartree-Fock method and perturbation theory.

# Arguments
- `Nocc::Int64`: Number of occupied orbitals.
- `h1b::Matrix{Float64,2}`: One-body Hamiltonian.
- `h2b::Vector{Float64,1}`: Two-body Hamiltonian flattened. This will be used only for the normal ordering.
- `gval::Float64`: Interaction strength.
- `to::TimeOutput`: TimeOutput object to measure the time.

# Optional arguments
- `debug_mode::Int64(0)`: If 1, print the energies at each iteration.
- `itnum_max::Int64(100)`: Maximum number of iterations.
- `tol::Float64(1e-9)`: Tolerance to stop the iteration.
"""
function HF(Nocc, h1b, h2b, gval, to, debug_mode=0; itnum_max=100, tol=1e-9)
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
        eval_Fock!(F, rho, h1b, Nocc, gval)
        @timeit to "solver: eigen Fock." evals, evecs = eigen(F)
        U .= evecs'
        e1b, e2b, EHF = eval_EHF(h1b, gval, rho, Nocc)
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

    # Developing functions to eval. natural orbitals (not used for now)
    eval_OnebodyDensityMatrix(h1b, gval, holes, particles)

    return EHF, EPT2, EPT3, HNO, holes, particles
end

"""

```math
\\begin{align}
\\rho_{pq} & = \\bra{\\Psi} c^\\dagger_p c_q \\ket{\\Psi}
\\ket{\\Psi} & \\approx \\ket{\\Psi^{(0)}} + \\ket{\\Psi^{(1)}} + \\ket{\\Psi^{(2)}}
\\end{align}
```

```math
\\rho \\approx \\rho^{(00)} + \\rho^{(11)} + \\rho^{(20)} + \\rho^{(02)}
```

where ``\\rho^{(00)}`` is the HF density matrix and the remaining terms are the correlation terms:

```math
\\begin{align}
\\rho^{(02)}_{pq} &= \\rho^{(20)}_{qp} = \\bra{\\Psi^{(0)}} c^\\dagger_p c_q \\ket{\\Psi^{(2)}} \\\\
\\rho^{(11)}_{pq} &=\\bra{\\Psi^{(1)}} c^\\dagger_p c_q \\ket{\\Psi^{(1)}} 
\\end{align}
```

These can be further decomposed into the following terms:

```math
\\begin{align}
\\rho^{(02)} &= D^{(A)} + D^{(B)} \\\\
\\rho^{(11)} &= D^{(C)} + D^{(D)} \\\\
D^{(A)}_{i'a'} &= \\frac{1}{2} \\sum_{abi} \\frac{H_{i'iab}H_{aba'i}}{(\\epsilon_{i'}-\\epsilon_{a'})(\\epsilon_{i'}+\\epsilon_{i}-\\epsilon_{a}-\\epsilon_{b})} \\\\
D^{(B)}_{i'a'} &= -\\frac{1}{2} \\sum_{aij} \\frac{H_{i'aij}H_{ija'a}}{(\\epsilon_{i'}-\\epsilon_{a'})(\\epsilon_{i}+\\epsilon_{j}-\\epsilon_{a}-\\epsilon_{a'})} \\\\
D^{(C)}_{i'j'} &= -\\frac{1}{2} \\sum_{abi} \\frac{H_{i'iab}H_{abj'i}}{(\\epsilon_{i'}+\\epsilon_{i}-\\epsilon_{a}-\\epsilon_{b})(\\epsilon_{j'}+\\epsilon_{i}-\\epsilon_{a}-\\epsilon_{b})} \\\\
D^{(D)}_{a'b'} & = \\frac{1}{2} \\sum_{aij} \\frac{H_{a'aij}H_{ijb'a}}{(\\epsilon_{i}+\\epsilon_{j}-\\epsilon_{a}-\\epsilon_{a'})(\\epsilon_{i}+\\epsilon_{j}-\\epsilon_{a}-\\epsilon_{b'})}
\\end{align}
```

The first two terms are exactly zero for global pairing Hamiltonian.
"""
function eval_OnebodyDensityMatrix(h1b, gval, holes, particles)
    D_hh = zeros(Float64, length(holes), length(holes))
    D_pp = zeros(Float64, length(particles), length(particles))
    
    for ip = 1:length(holes)
        jp = ip # through Γ_{i'iab}Γ_{abj'i}, i = bar(i') = bar(j')
        i = ip + ifelse(ip%2==0, -1, 1)
        tmp = 0.0
        for a in particles
            b = a + ifelse(a%2==0, -1, 1) # b must be bar(a)
            nume = (-gval)^2
            deno = h1b[i,i] + h1b[ip,ip] - h1b[a,a] - h1b[a,a]
            deno = deno * (h1b[jp,jp] + h1b[i,i] - h1b[a,a] - h1b[b,b])
            tmp += nume/deno
        end
        D_hh[ip, jp] = -0.5 * tmp + 1.0 # 1.0 is the HF density matrix
    end

    for (idx_ap,ap) in enumerate(particles)
        bp = ap # through Γ_{a'aij}Γ_{ijb'a}, a = bar(a') = bar(b')
        idx_bp = idx_ap
        a = ap + ifelse(ap%2==0, -1, 1)
        tmp = 0.0
        for i in holes
            j = i + ifelse(i%2==0, -1, 1) # j must be bar(i)
            nume = (-gval)^2
            deno = h1b[i,i] + h1b[j,j] - h1b[a,a] - h1b[ap,ap]
            deno = deno * (h1b[i,i] + h1b[j,j] - h1b[a,a] - h1b[bp,bp])
            tmp += nume/deno
        end
        D_pp[idx_ap, idx_bp] = 0.5 * tmp 
    end
    @assert tr(D_hh) + tr(D_pp) ≈ length(holes)
    return D_hh, D_pp
end
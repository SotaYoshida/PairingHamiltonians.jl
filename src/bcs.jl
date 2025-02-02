"""
    uv_from_Lam_Delta!(Delta, lambda, epsilon_p::Vector{Float64}, us::Vector{Float64}, vs::Vector{Float64})

Destructively update the Bogoliubov coefficients `us` and `vs` from the gap parameter `Delta` and the chemical potential `lambda`:

```math
u^2_a = \\frac{1}{2} \\left( 1 + \\frac{\\epsilon'_a - \\lambda}{\\sqrt{(\\epsilon'_a - \\lambda)^2 + \\Delta^2}} \\right)  \\\\
v^2_a = \\frac{1}{2} \\left( 1 - \\frac{\\epsilon'_a - \\lambda}{\\sqrt{(\\epsilon'_a - \\lambda)^2 + \\Delta^2}} \\right) 
```

# Arguments
- `Delta::Float64`: gap parameter
- `lambda::Float64`: chemical potential
- `epsilon_p::Vector{Float64}`: `epsilon_p[i] = epsilon[i] - g v[i]^2`
- `us::Vector{Float64}`: `u` coefficients of Bogoliubov transformation
- `vs::Vector{Float64}`: `v` coefficients of Bogoliubov transformation
"""
function uv_from_Lam_Delta!(Delta, lambda, epsilon_p::Vector{Float64}, us::Vector{Float64}, vs::Vector{Float64})
    for i = 1:length(epsilon_p)
        ep = epsilon_p[i]
        us[i] = sqrt( 0.5 * ( 1.0 + (ep - lambda) / denom(ep, lambda, Delta) ) )
        vs[i] = sqrt( 0.5 * ( 1.0 - (ep - lambda) / denom(ep, lambda, Delta) ) )
        @assert us[i]^2 + vs[i]^2 ≈ 1.0 "u2 + v2 (i=$i)= $(us[i]^2 + vs[i]^2)"
    end
    return nothing
end

"""
    eval_eps_p!(epsilon_p, epsilon, vs, gval)

Calculate the `epsilon_p` from `epsilon` and `vs`.

# Arguments
- `epsilon_p::Vector{Float64}`: `epsilon_p[i] = epsilon[i] - g v[i]^2`
- `epsilon::Vector{Float64}`: `epsilon[i]`
- `vs::Vector{Float64}`: `v[i]`
"""
function eval_eps_p!(epsilon_p, epsilon, vs, gval)
    for i in 1:length(epsilon)
        epsilon_p[i] = epsilon[i] - gval * vs[i]^2
    end
    return nothing
end

"""
    eval_Delta(vs, us, gval)

Calculate the gap parameter `Delta`: `` \\Delta \\equiv \\frac{g}{2} \\sum_{\\alpha} u_a v_a ``

# Arguments
- `vs::Vector{Float64}`: `v[i]`
- `us::Vector{Float64}`: `u[i]`
- `gval::Float64`: `g`
"""
function eval_Delta(vs, us, gval)::Float64
    delta = 0.5 * gval  * dot(vs, us)
    return delta
end

function denom(epsilon, lambda, Delta)
    return sqrt( (epsilon - lambda)^2 + Delta^2 )
end

"""
    reeval_Delta(lambda, Delta, epsilon_p, gval, debug_mode, alpha=0.5)::Float64

Update the gap parameter `Delta` iteratively with the following equation:

```math
\\Delta^\\mathrm{new} := \\alpha \\Delta^\\mathrm{new} + (1-\\alpha) \\frac{g}{4} \\sum_{\\alpha} \\frac{\\Delta}{\\sqrt{(\\epsilon'_a - \\lambda)^2 + \\Delta^2}}
```

# Arguments
- `lambda::Float64`: chemical potential
- `Delta::Float64`: current gap parameter
- `epsilon_p::Vector{Float64}`: `epsilon_p[i] = epsilon[i] - g v[i]^2`
- `gval::Float64`: `g`
- `debug_mode::Int`: specify the debug mode
- `alpha::Float64`: mixing parameter to update `Delta`, which may be useful to stabilize the calculation
"""
function reeval_Delta(lambda, Delta, epsilon_p, gval, debug_mode, alpha=0.5)::Float64
    ret = 0.0
    for i in eachindex(epsilon_p)
        ret += 1.0 / denom(epsilon_p[i], lambda, Delta)
    end
    ret = alpha * Delta + (1-alpha) * Delta * ret * gval/ 4 
    if debug_mode > 1
        println(@sprintf("updated Delta %7.4f (<=%7.4f)", ret, Delta))
    end
    return ret
end

"""
    Neval(epsilon_p, lambda, Delta)::Float64

Calculate the expectation value of the particle number from the chemical potential `lambda` and the gap parameter `Delta` via the following equation:

```math
N  = \\sum_{\\alpha} \\frac{1}{2} \\left( 1 - \\frac{\\epsilon'_a - \\lambda}{\\sqrt{(\\epsilon'_a - \\lambda)^2 + \\Delta^2}} \\right)
```
"""
function Neval(epsilon_p, lambda, Delta)::Float64
    tot = 0.0
    for i = 1:length(epsilon_p)
        tot += 0.5 
        tot += -0.5 * (epsilon_p[i] - lambda) / denom(epsilon_p[i], lambda, Delta)
    end
    return tot
end

"""
    update_lambda(Nocc, lambda, Delta, epsilon_p, gval, debug_mode, alpha=0.5)

Update the chemical potential `lambda` iteratively with the following equation:

```math
\\lambda^\\mathrm{new} := \\frac{g}{2} \\left( N_\\mathrm{occ.} - \\frac{N_\\mathrm{orb.}}{2} + \\frac{1}{2} \\sum_{\\alpha} \\frac{\\epsilon'_a }{\\sqrt{(\\epsilon'_a - \\lambda)^2 + \\Delta^2}} \\right)
```

# Arguments
- `Nocc::Int`: number of occupied states
- `lambda::Float64`: chemical potential
- `Delta::Float64`: gap parameter
- `epsilon_p::Vector{Float64}`: `epsilon_p[i] = epsilon[i] - g v[i]^2`
- `gval::Float64`: `g`
- `debug_mode::Int`: specify the debug mode
- `alpha::Float64`: mixing parameter to update `lambda`
"""
function update_lambda(Nocc, lambda, Delta, epsilon_p, gval, debug_mode, alpha=0.5)
    Norb = length(epsilon_p)
    lambda_new = denosum = 0.0    
    for i in eachindex(epsilon_p)
        deno = denom(epsilon_p[i], lambda, Delta)
        denosum += 1.0 / deno
        lambda_new += epsilon_p[i] / deno  
    end
    lambda_new = alpha * lambda + (1-alpha) * gval * (lambda_new/2 + Nocc - Norb/2) / 2
    if debug_mode > 0
        println("lambda_new = $lambda_new |1-g/4*denosum| = $(@sprintf("%9.1e", abs(1.0-gval*denosum/4)))")
    end
    return lambda_new
end

"""
    E_BCS_from_uv(epsilon, epsilon_p, vs, Delta, lambda, gval, debug_mode)

Calculate the ground state energy of the BCS Hamiltonian from the Bogoliubov coefficients `vs` and `us`:

```math
E_\\mathrm{BCS} = \\sum_{\\alpha} \\epsilon'_a  v^2_a - \\frac{\\Delta^2}{g}
```
"""
function E_BCS_from_uv(epsilon, epsilon_p, vs, Delta, lambda, gval, debug_mode)
    E_BCS = term_0 = term_1 = 0.0
    for i in 1:length(epsilon)
        term_0 += epsilon[i] * vs[i]^2
        term_1 += - gval *vs[i]^4
        E_BCS += (   epsilon[i] -  0.5 * gval * vs[i]^2 )* vs[i]^2
    end
    term_2 = - Delta^2 / gval
    E_BCS += term_2 
    if debug_mode > 0 
        println(@sprintf("E = %7.4f (Σe_i v2i %7.4f, 1b %7.4f, 2b %7.4f) ", E_BCS, term_0, term_1, term_2))
        println(@sprintf("<N> sum(v^2) %5.2f sigma %5.2f ", sum(vs.^2),  Neval(epsilon_p, lambda, Delta)))
    end
    return E_BCS
end

"""
    _main_BCS(Norb, Nocc, epsilon, gval, to, debug_mode; max_iter=1000)

Main function to calculate the ground state energy of the BCS Hamiltonian.

# Arguments
- `Norb::Int`: number of orbitals
- `Nocc::Int`: number of occupied states
- `epsilon::Vector{Float64}`: single-particle energies
- `gval::Float64`: interaction strength
- `to::TimerOutput`: timer object to measure the elapsed time
- `debug_mode::Int`: specifying the debug mode

# Optional arguments
- `max_iter::Int(1000)`: maximum number of iterations
"""
function _main_BCS(Norb, Nocc, epsilon, gval, to, debug_mode;
                     max_iter=1000)
    Random.seed!(123)
    # initialize Bogoliubov coefficients
    thetas = [ rand() * pi/2 for _ in 1:Norb ] 

    us = [ cos(theta) for theta in thetas ]
    vs = [ sin(theta) for theta in thetas ]
    lambda = ( epsilon[Nocc]  + epsilon[Nocc+1]) / 2.0
    Delta = max(1.e-6, eval_Delta(vs, us, gval))
    if debug_mode > 0
        show_vector("initial v2s", vs.^2, 2)           
        println("initial Delta = $Delta")
    end
    E_BCS = 0.0

    # It is convenient to define the array to store epsilon' = epsilon - g v^2
    epsilon_p = zeros(Float64, length(epsilon))
    eval_eps_p!(epsilon_p, epsilon, vs, gval)
    E_BCS_prev = E_BCS = E_BCS_from_uv(epsilon, epsilon_p, vs, Delta, lambda, gval, debug_mode)

    # Iteration for self-consistent calculation
    for iter = 1:max_iter
        if debug_mode > 0
            println(@sprintf("------------ iter = %4i ------------", iter))
        end
        # calculate Delta
        Delta = reeval_Delta(lambda, Delta, epsilon_p, gval, debug_mode)
        uv_from_Lam_Delta!(Delta, lambda, epsilon_p, us, vs)
        eval_eps_p!(epsilon_p, epsilon, vs, gval)

        # calculate lambda
        lambda = update_lambda(Nocc, lambda, Delta, epsilon_p, gval, debug_mode)
        uv_from_Lam_Delta!(Delta, lambda, epsilon_p, us, vs)
        eval_eps_p!(epsilon_p, epsilon, vs, gval)
        
        # calculate E_BCS
        E_BCS = E_BCS_from_uv(epsilon, epsilon_p, vs, Delta, lambda, gval, debug_mode)

        # check convergence
        if abs(E_BCS - E_BCS_prev) < 1.e-6
            break
        end
        E_BCS_prev = E_BCS
        if debug_mode > 0
            println("")
        end
    end

    if debug_mode > 0
        println("-------------------------------------------")
        println("fermi surface $lambda")
        show_vector("v2s", vs.^2, 2)           
        show_vector("epslion", epsilon, 2)
        show_vector("epsilon_p", epsilon_p, 2)
    end

    # Check the number of particles
    N = sum(vs.^2)
    if abs(N - Nocc) > 1.e-4
        println("N_BCS $N E $E_BCS")
        E_BCS = NaN
    end
    return E_BCS 
end
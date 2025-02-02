"""
    PT2(F, Nocc, gval)

Function to calculate the second order perturbation theory contribution to the energy.

```math
\\frac{1}{4} \\sum^{\\leq \\epsilon_F}_{ij}\\sum^{\\geq \\epsilon_F}_{pq} \\frac{ \\langle ij | V | pq \\rangle \\langle pq | V | ij \\rangle
}{\\epsilon_i + \\epsilon_j - \\epsilon_p - \\epsilon_q} 
```
For the pairing Hamiltonian, `j` should be the pair of `i` and `p` should be the pair of `q`, and the numerator will be `g^2`.
"""
function PT2(F, Nocc, gval)
    Nq = size(F, 1)
    EPT2 = 0.0
    for i = 1:Nocc #hole
        for p = Nocc+1:Nq #particle 
            denominator = 8* (F[i,i] - F[p, p])
            EPT2 += gval^2 /denominator
        end
    end
    return EPT2
end

"""
    PT3(F, Nocc, gval)

Function to calculate the third order perturbation theory contribution to the energy.

```math
\\begin{align}
\\Delta E^{(3)}_{pp} &
= \\frac{1}{8} \\sum_{ij \\leq F} \\sum_{abcd > F} \\frac{\\bra{ij} \\hat{V} \\ket{ab} \\bra{ab} \\hat{V} \\ket{cd} \\bra{cd} \\hat{V} \\ket{ij}}{
\\epsilon^{ab}_{ij} \\epsilon^{cd}_{ij}}  \\\\
\\Delta E^{(3)}_{hh} &
= \\frac{1}{8} \\sum_{ijkl \\leq F} \\sum_{ab > F} \\frac{\\bra{ij} \\hat{V} \\ket{ab} \\bra{kl} \\hat{V} \\ket{kl} \\bra{kl} \\hat{V} \\ket{ij}}{
    \\epsilon^{ab}_{ij} \\epsilon^{ab}_{kl}} \\\\
\\Delta E^{(3)}_{ph} &
= - \\sum_{ijk \\leq F} \\sum_{abc > F} \\frac{\\bra{ij} \\hat{V} \\ket{ab} \\bra{kb} \\hat{V} \\ket{ic} \\bra{ac} \\hat{V} \\ket{kj}}{
    \\epsilon^{ab}_{ij} \\epsilon^{ac}_{kj}}
\\end{align}
```

For the global pairing Hamiltonian, the PT3 is much simpler than the general case as in the PT2,
since the numerator will be proportional to `g^3` and the `ph` term is zero.
"""
function PT3(F, Nocc, gval; verbose=false)
    Nq = size(F, 1)
    Epp = Ehh = Eph = 0.0
    nume = (-gval)^3

    for i = 1:Nocc
        # Epp
        for p = Nocc+1:Nq
            for q = Nocc+1:Nq
                deno = 32 * (F[i,i] - F[p,p]) * (F[i,i] - F[q,q])
                tmp = nume /deno
                Epp += nume/deno
            end
        end
        # Ehh
        for j = 1:Nocc
            for p = Nocc+1:Nq
                deno = 32 * (F[j,j] - F[p,p]) * (F[i,i] - F[p,p])
                Ehh += nume/deno
            end
        end
    end
    EPT3 = Epp + Ehh + Eph
    return EPT3
end
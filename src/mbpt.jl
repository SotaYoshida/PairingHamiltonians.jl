"""
```math
\\frac{1}{4} \\sum^{\\leq \\epsilon_F}_{ij}\\sum^{\\geq \\epsilon_F}_{pq} \\frac{g^2}{\\epsilon_i + \\epsilon_j - \\epsilon_p - \\epsilon_q} 
```
For the pairing Hamiltonian, `j` should be the pair of `i` and `p` should be the pair of `q`.
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
For the simple pairing Hamiltonian, the PT3 contribution is again much simpler than the general case.
Eph is to be zero.
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
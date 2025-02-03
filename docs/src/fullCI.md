# Full-CI

Functions to solve the pairing Hamiltonian with the Full-CI method.
One can specify the `solver` for main API, i.e. `Pairing_Hamiltonian` function in `src/hamiltonian.jl`, to be `"Full-CI"` and `"Full-CI(2-fold)"`.

The latter method is more efficient than the former, but it is only applicable to pairing Hamiltonian with 2-fold degenerate single-particle states.
The former method is more general and applicable to e.g. odd-number of particles.

Naively, the dimension of the Hamiltonian matrix scales ${}_{N_\mathrm{orb}} C_{N_\mathrm{occ}}$, where $N_\mathrm{orb}$ is the number of single-particle states and $N_\mathrm{occ}$ is the number of particles.
For half-filled case, current limitation in a laptop may be around $N_\mathrm{orb} \approx 30$.

```@autodocs
Modules = [PairingHamiltonian]
Pages = ["fullCI.jl"]
``` 


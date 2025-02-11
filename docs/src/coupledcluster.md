# Coupled Cluster

Functions for performing coupled cluster calculations:

```math
\begin{align*}
\ket{\Psi_\mathrm{CC}} & = e^{\hat{T}} \ket{\Phi_0}\\
\hat{T} & = \hat{T}_1 + \hat{T}_2 + \hat{T}_3 + \cdots + \hat{T}_A \\
\hat{T}_n & = \left( \frac{1}{n!} \right)^2
\sum^{a_1 a_2 \ldots a_n}_{i_1 i_2 \ldots i_n} 
t^{a_1 a_2 \ldots a_n}_{i_1 i_2 \ldots i_n}
a^\dagger_{a_1} a^\dagger_{a_2} \ldots a^\dagger_{a_n}
a_{i_n} a_{i_{n-1}} \ldots a_{i_1}
\end{align*}
```

Since the current Hamiltonian is a pairing Hamiltonian, we will only consider the CCD method, $\hat{T} = \hat{T}_2$.

```@autodocs
Modules = [PairingHamiltonians]
Pages = ["coupled_cluster.jl"]
``` 


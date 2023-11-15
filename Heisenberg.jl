## 11/06/2023
# Obtain the ground-state of the Heisenberg model and save the wavefuction 
# Use the ground-state wavefunction as the input for the disentangler

using Pkg
using ITensors
using HDF5

let
    N = 100
    sites = siteinds("S=1/2", N)

    os = OpSum()
    for index = 1 : N - 1
        os += "Sz", index, "Sz", index + 1
        os += 1/2, "S+", index, "S-", index + 1
        os += 1/2, "S-", index, "S+", index + 1
    end

    Hamiltonian = MPO(os, sites)
    ψ₀ = randomMPS(sites, 10)
    nsweeps = 10
    maximum_bond_dimension = [10, 20, 100, 100, 200, 200, 500, 500, 800, 800]
    cutoff = [1E-8]

    energy, ψ = dmrg(Hamiltonian,ψ₀; nsweeps, maximum_bond_dimension, cutoff)
    @show energy
    
    h5open("Heisenberg_wave_function.h5", "w") do file
        write(file, "psi", ψ)
    end
    
    return
end

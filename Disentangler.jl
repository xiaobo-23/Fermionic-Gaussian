## 11/15/2023
# Implement the disentangler to perform basis rotation

using Pkg
using ITensors
using HDF5
using Random 
using NLopt

include("generators.jl")
include("obtain_probability.jl")
include("entanglement_entropy.jl")

let 
    # # Load the wavefunction from the file
    # ψ = h5open("Heisenberg_wave_function.h5", "r") do file
    #     read(file, "psi")
    # end
    # println(typeof(ψ))

    # Obtain the ground state of the Heisenberg model
    #************************************************************************************
    N = 5
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
    maximum_bond_dimension = [10, 20, 100, 100, 200, 200, 500, 500, 800, 1000]
    cutoff = [1E-8]

    energy, ψ = dmrg(Hamiltonian,ψ₀; nsweeps, maximum_bond_dimension, cutoff)
    @show energy
    #************************************************************************************

    # Define the parameters to set up the disentangler 
    subsystem_size = 5           # System size to compute the correlation matrix 
    number_of_generators = 15    # This number is fixed by the number of SU(4) generators

    # The parameters to be optimized
    # Initialize the parameters to be one
    alpha = ones(subsystem_size - 1, number_of_generators)    
    
    # Random initialization of the parameters
    # random_seed=0
    # Random.seed!(random_seed)
    # upper_bound = 5.0
    # lower_bound = -5.0
    # alpha = rand(subsystem_size - 1, number_of_generators) * (upper_bound - lower_bound) .+ lower_bound     # The parameters to be optimized
    # @show alpha
    # @show typeof(alpha)

    # Define the generators of SU(4) group
    lambda_set = [lambda₁, lambda₂, lambda₃, lambda₄, lambda₅, lambda₆, lambda₇, lambda₈, 
    lambda₉, lambda₁₀, lambda₁₁, lambda₁₂, lambda₁₃, lambda₁₄, lambda₁₅]

    # Define the observables to be measured
    SvN = zeros(Float64, 2, N - 1)  # N-1 is the number of bonds 
    Sz  = zeros(Float64, 2, N)      # N is the number of sites

    # Build up the gates for the disentangler
    gates = ITensor[]
    for index = subsystem_size - 1 : -1 : 1
        @show index
        tmp_matrix = zeros(4, 4)

        for generator_index in 1 : number_of_generators
            @show alpha[index, generator_index] * lambda_set[generator_index]
            tmp_matrix += alpha[index, generator_index] * lambda_set[generator_index]
        end 
        # @show tmp_matrix
        s1 = sites[index]
        s2 = sites[index + 1]
        tmp_gate = itensor(exp(im * tmp_matrix), s1', s2', s1, s2)

        # @show typeof(tmp_gate)
        push!(gates, tmp_gate)
    end

    # Apply the generalized gates to the MPS
    ψ_copy = deepcopy(ψ)
    SvN[1, :] = entanglement_entropy(ψ_copy, N)
    Sz[1, :] = expect(ψ_copy, "Sz"; sites = 1 : N)
    
    ψ_copy = apply(gates, ψ_copy; cutoff = 1e-8)
    Sz[2, :] = expect(ψ_copy, "Sz"; sites = 1 : N)
    SvN[2, :] = entanglement_entropy(ψ_copy, N)
    @show Sz[1, :]
    @show Sz[2, :]

    tmp_probability = project_probability(ψ_copy, 1, "Sz")
    cost_function = tmp_probability[1, 2]
    @show tmp_probability[1, 1], tmp_probability[1, 2], tmp_probability[1, 1] + tmp_probability[1, 2] == 1
    @show SvN[1, :]
    @show SvN[2, :]


    ## 12/06/2023
    # Implement the optimization process for the first gate in a reverse order 
    psi_bra = deepcopy(ψ)
    psi_ket = deepcopy(ψ)
    
    # Apply all the gates in the original order
    @show siteinds(psi_ket)
    psi_ket = apply(gates, psi_ket; cutoff = 1e-8)
    @show siteinds(psi_ket)

    # Apply all the gates except for the gate that is being optimized
    reverse_gates = reverse(gates)
    psi_ket = apply(reverse_gates[1 : subsystem_size - 2], psi_ket; cutoff = 1e-8)
    @show length(reverse_gates)
    
    # @show gates, reverse_gates
    # @show siteinds(psi_ket)
    # @show siteinds(psi_bra)
    @show length(psi_bra), length(psi_ket)
    @show inner(psi_bra, psi_ket)

    # h5open("Data/disentangler_N100_random$(random_seed).h5", "w") do file
    #     write(file, "psi", ψ_copy)
    #     write(file, "SvN", SvN)
    #     write(file, "Sz", Sz)
    # end

    return
end
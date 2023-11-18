## 11/15/2023
# Implement the disentangler to perform basis rotation

using Pkg
using ITensors
using HDF5
using NLopt

include("generators.jl")
include("obtain_probability.jl")

let 
    # # Load the wavefunction from the file
    # ψ = h5open("Heisenberg_wave_function.h5", "r") do file
    #     read(file, "psi")
    # end

    # println(typeof(ψ))

    # Obtain the ground state of the Heisenberg model
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
    maximum_bond_dimension = [10, 20, 100, 100, 200, 200, 500, 500, 800, 1000]
    cutoff = [1E-8]

    energy, ψ = dmrg(Hamiltonian,ψ₀; nsweeps, maximum_bond_dimension, cutoff)
    @show energy

    # Define the parameters to set up the disentangler 
    subsystem_size = 5          # System size to compute the correlation matrix 
    number_of_generators = 15   # This number is fixed by the number of SU(4) generators
    # alpha = ones(subsystem_size - 1, number_of_generators)      # The parameters to be optimized
    lambda_set = [lambda₁, lambda₂, lambda₃, lambda₄, lambda₅, lambda₆, lambda₇, lambda₈, 
    lambda₉, lambda₁₀, lambda₁₁, lambda₁₂, lambda₁₃, lambda₁₄, lambda₁₅]
    # @show typeof(alpha)

    # # Build up the gates for the disentangler
    # gates = ITensor[]
    # for index = subsystem_size - 1 : -1 : 1
    #     @show index
    #     tmp_matrix = zeros(4, 4)

    #     for generator_index in 1 : number_of_generators
    #         @show alpha[index, generator_index] * lambda_set[generator_index]
    #         tmp_matrix += alpha[index, generator_index] * lambda_set[generator_index]
    #     end 
    #     @show tmp_matrix
    #     s1 = sites[index]
    #     s2 = sites[index + 1]
    #     tmp_gate = itensor(exp(im * tmp_matrix), s1', s2', s1, s2)

    #     @show typeof(tmp_gate)
    #     push!(gates, tmp_gate)
    # end

    # # Apply the list of the gates to the wavefunction
    # ψ_copy = deepcopy(ψ)
    # ψ_copy = apply(gates, ψ_copy; cutoff = 1e-8)
    # # @show expect(ψ, "Sz"; sites = 1 : N)
    # # @show expect(ψ_copy, "Sz"; sites = 1 : N)

    # tmp_probability = sample(ψ_copy, 1, "Sz")
    # cost_function = tmp_probability[1, 2]

    alpha = ones(Float64, (subsystem_size - 1) * number_of_generators)
    @show typeof(alpha)
    function disentangler(tmp_alpha :: Vector)
        # Build up the gates for the disentangler
        gates = ITensor[]
        for index = subsystem_size - 1 : -1 : 1
            @show index
            tmp_matrix = zeros(4, 4)

            for generator_index in 1 : number_of_generators
                @show tmp_alpha[(index - 1) * number_of_generators + generator_index] * lambda_set[generator_index]
                tmp_matrix += tmp_alpha[(index - 1) * number_of_generators + generator_index] * lambda_set[generator_index]
            end 
            @show tmp_matrix
            s1 = sites[index]
            s2 = sites[index + 1]
            tmp_gate = itensor(exp(im * tmp_matrix), s1', s2', s1, s2)

            @show typeof(tmp_gate)
            push!(gates, tmp_gate)
        end

        # Apply the list of the gates to the wavefunction
        ψ_copy = apply(gates, ψ_copy; cutoff = 1e-8)
        # @show expect(ψ, "Sz"; sites = 1 : N)
        # @show expect(ψ_copy, "Sz"; sites = 1 : N)

        tmp_probability = sample(ψ_copy, 1, "Sz")
        # cost_function = tmp_probability[1, 1]
        @show tmp_probability[1, 2]
        return tmp_probability[1, 2]
    end

    ψ_copy = deepcopy(ψ)
    opt = Opt(:LN_BOBYQA, (subsystem_size - 1) * number_of_generators)
    opt.min_objective = disentangler
    (minf, minx, ret) = optimize(opt, alpha)
    numevals = opt.numevals
    println("got $minf at $minx after $numevals iterations (retunrned $ret))")


    # tmp_probability = disentangler(alpha) 
    # @show tmp_probability[1, 1], tmp_probability[1, 2], tmp_probability[1, 1] + tmp_probability[1, 2] == 1
    return
end
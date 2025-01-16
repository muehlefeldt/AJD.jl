module AJD 
using LinearAlgebra
using BenchmarkTools

# Import different algorithms.
include("jdiag_algorithms/jdiag_cardoso.jl")
include("jdiag_algorithms/jdiag_gabrieldernbach.jl")
include("jdiag_algorithms/jdiag_edourdpineau.jl")
include("jdiag_algorithms/FFDiag.jl")

# Utility functions, plotting functions and global constanst imported.
include("utils.jl")
include("plotting.jl")
include("global_constants.jl")



"""
    diagonalize(
        A::Vector{<:AbstractMatrix{<:Number}};
        algorithm::String = "jdiag_gabrieldernbach",
        max_iter::Int = 1000,
        threshold::AbstractFloat = eps()
        )

Calculate joint diagonalization of multiple input matrices using the requested algorithms.

Main function of the AJD package.
Implemented algorithms at this point in time are limited to the [JDiag algorithm](https://doi.org/10.1137/S0895479893259546)  in different versions.
Input of matrices to be diagonalized need to be a vector of matrices.
The matrices can be of types Float64 or Complex. Limitations of the different implementations apply.

Supported algorithms are `jdiag_gabrieldernbach`, `jdiag_cardoso` and `jdiag_edourdpineau`.
See the Getting Started Guide for information on the algorithms.
"""
function diagonalize(
    A::Vector{<:AbstractMatrix{<:Number}};
    algorithm::String = "jdiag_gabrieldernbach",
    max_iter::Int = 1000,
    threshold::AbstractFloat = eps(),
    plot_matrix::Bool = false,
    plot_convergence::Bool = false
    )

    if !check_input(A)
        throw(ArgumentError("Invalid input."))
    end

    if algorithm in ["jdiag", "jdiag_gabrieldernbach"]
        F, B, error_array = jdiag_gabrieldernbach!(A, max_iter = max_iter, threshold = threshold, plot_convergence = plot_convergence)

    elseif algorithm == "jdiag_edourdpineau"
        F, B, error_array = jdiag_edourdpineau(A, iter = max_iter)

    elseif algorithm == "jdiag_cardoso"
        if typeof(A) <: AbstractArray{<:AbstractArray{<:Real}} 
            F, B, error_array = jdiag_cardoso(A, threshold, plot_convergence = plot_convergence)
        else
            throw(ArgumentError("Not supported for set of Matrices containing imaginary values!"))
        end

    elseif algorithm in ["FFD", "ffd", "ffdiag"]
        F, B, error_array = FFD!(
            A,
            threshold=threshold,
            max_iter=max_iter,
            plot_convergence=plot_convergence
        )

    else
        # If no vaild algorithm selected, throw an error.
        throw(ArgumentError("No valid algorithm selected."))
    end

    # Plotting output if so selected by the user.
    if plot_matrix
        # Illustrate Filter and diagonlised matrices.
        display(plot_matrix_heatmap(F, B))
    end

    if plot_convergence
        # Show convergence of the error.
        display(plot_convergence_lineplot(error_array, algorithm))
    end

    return create_linear_filter(F)
end

"""
    ajd_benchmark(n_dims::Int, n_matrices::Int)

Run benchmark of implemented algorithms with random inputs.
Prints basic overview of median execution times.
Returns BenchmarkGroup containing detailed results.
"""
function ajd_benchmark(n_dims::Int, n_matrices::Int)
    # Define a parent BenchmarkGroup to contain our suite
    suite = BenchmarkGroup()
    for name in ["jdiag_gabrieldernbach", "jdiag_cardoso", "jdiag_edourdpineau"]
        suite[name] = BenchmarkGroup(["jdiag"])
        # Set the function to be benchmarked.
        suite[name]["real"] = begin
            @benchmarkable diagonalize(data, algorithm=$name) setup=(data=AJD.random_normal_commuting_matrices($n_dims, $n_matrices)) 
        end
    end

    name = "ffdiag"
    suite[name] = BenchmarkGroup([name])
    suite[name]["real"] = begin
        @benchmarkable diagonalize(data, algorithm=$name) setup=(data=AJD.random_normal_commuting_matrices($n_dims, $n_matrices)) 
    end

    # Run the actual benchmark.
    tune!(suite)
    results = run(suite, verbose = true)

    # Print basic overview.
    for name in ["jdiag_gabrieldernbach", "jdiag_cardoso", "jdiag_edourdpineau"]
        print(name)
        print(median(results[name]["real"]))
    end

    # Return BenchmarkGroup for further evaluation.
    return results
end

export diagonalize, ajd_benchmark

end



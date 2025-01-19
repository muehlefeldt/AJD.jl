module AJD 
using LinearAlgebra
using BenchmarkTools
using Plots: Plot

# Import different algorithms.
include("jdiag_algorithms/jdiag_cardoso.jl")
include("jdiag_algorithms/jdiag_gabrieldernbach.jl")
include("jdiag_algorithms/jdiag_edourdpineau.jl")
include("jdiag_algorithms/FFDiag.jl")

# Utility functions, plotting functions and global constanst imported.
include("utils.jl")
include("utils_test_data.jl")
include("plotting.jl")
include("global_constants.jl")



"""
    diagonalize(
        A::Vector{<:AbstractMatrix{<:Number}};
        algorithm::String = "jdiag_gabrieldernbach",
        max_iter::Int = 1000,
        threshold::AbstractFloat = eps(),
        plot_matrix::Bool = false,
        plot_convergence::Bool = false
        )

Calculate joint diagonalization of multiple input matrices ``M_k``.

Main function of the AJD package.
Implemented algorithms are [JDiag](https://doi.org/10.1137/S0895479893259546) and FFDiag.
Input of matrices ``M_k`` need to be a vector of matrices.
The matrices can be of types Float64 or Complex.

Supported algorithms are `jdiag_gabrieldernbach`, `jdiag_cardoso` and `jdiag_edourdpineau`.
See the Getting Started Guide for information on the algorithms. Test

``M_k``.

"""
function diagonalize(
    A::Vector{<:AbstractMatrix{<:Number}};
    algorithm::String = "jdiag_gabrieldernbach",
    max_iter::Int = 1000,
    threshold::AbstractFloat = eps(),
    )::LinearFilter

    F, _, _ = get_diagonalization(A, algorithm=algorithm, max_iter=max_iter, threshold=threshold, only_plot=:no_plot)
    return create_linear_filter(F) 
end

function diagonalize(
    A::Vector{<:AbstractMatrix{<:Number}},
    only_plot::Symbol;
    algorithm::String = "jdiag_gabrieldernbach",
    max_iter::Int = 1000,
    threshold::AbstractFloat = eps(),
    )::Plot

    if only_plot == :plot
        F, B, error_array = get_diagonalization(
            A, 
            algorithm=algorithm,
            max_iter=max_iter, 
            threshold=threshold,
            only_plot=only_plot
        )
        p = get_plot(F, B, error_array, algorithm)
    else
        throw(ArgumentError("Please use symbol ony_plot=:plot to generate plots."))
    end
    return p
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



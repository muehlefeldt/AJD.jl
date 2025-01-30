module AJD
using LinearAlgebra: eigen, norm, Symmetric, Hermitian, I, qr, dot, diag
using BenchmarkTools
using ProgressMeter

# Import different algorithms.
include("jdiag_algorithms/jdiag_cardoso.jl")
include("jdiag_algorithms/jdiag_gabrieldernbach.jl")
include("jdiag_algorithms/jdiag_edourdpineau.jl")
include("FFDiag.jl")

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

Supported algorithms are `jdiag_gabrieldernbach`, `jdiag_cardoso`, `jdiag_edourdpineau` and `ffdiag`.
See the Getting Started Guide for information on the algorithms.

"""
function diagonalize(
    A::Vector{<:AbstractMatrix{<:Number}};
    algorithm::String = "jdiag_gabrieldernbach",
    max_iter::Int = 1000,
    threshold::AbstractFloat = eps(),
)
    check_input(A, max_iter, threshold)

    F, _, _, n_iter = get_diagonalization(
        A,
        algorithm = algorithm,
        max_iter = max_iter,
        threshold = threshold,
        only_plot = :no_plot,
    )
    
    if n_iter >= max_iter
        @warn "Max iteration was reached. Consider increasing max_iter: diagonalize(M, max_iter=...)."
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

    for name in ["jade", "ffdiag"]
        suite[name] = BenchmarkGroup([name])
        suite[name]["exact_diag"] = begin
            @benchmarkable diagonalize(data, algorithm = $name) setup = (
                data = AJD.get_test_data(
                    :exact_diag,
                    n_dims = $n_dims,
                    n_matrices = $n_matrices,
                )
            )
        end
        suite[name]["approx_diag_large"] = begin
            @benchmarkable diagonalize(data, algorithm = $name) setup = (
                data = AJD.get_test_data(
                    :approx_diag_large,
                    n_dims = $n_dims,
                    n_matrices = $n_matrices,
                )
            )
        end
        suite[name]["random"] = begin
            @benchmarkable diagonalize(data, algorithm = $name) setup = (
                data = AJD.get_test_data(
                    :random_noise,
                    n_dims = $n_dims,
                    n_matrices = $n_matrices,
                )
            )
        end
    end

    # Run the actual benchmark.
    tune!(suite)
    results = run(suite, verbose = true)

    # Return BenchmarkGroup for further evaluation.
    return results
end

export diagonalize, ajd_benchmark

end

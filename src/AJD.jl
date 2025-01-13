module AJD 
using LinearAlgebra
using BenchmarkTools

# Import different JDiag algorithms.
include("jdiag_algorithms/jdiag_cardoso.jl")
include("jdiag_algorithms/jdiag_gabrieldernbach.jl")
include("jdiag_algorithms/jdiag_edourdpineau.jl")
include("jdiag_algorithms/FFDiag.jl")
# Utility function import.
include("utils.jl")

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
    threshold::AbstractFloat = eps()
    )

    if !check_input(A)
        throw(ArgumentError("Invalid input."))
    end

    if algorithm in ["jdiag", "jdiag_gabrieldernbach"]
        A,F = jdiag_gabrieldernbach!(A, max_iter = max_iter, threshold = threshold)
        @info "Jdiag", A
        return AJD.create_linear_filter(F)
    end

    if algorithm == "jdiag_cardoso"
        if typeof(A) <: AbstractArray{<:AbstractArray{<:Real}} 
            _,F ,_ = jdiag_cardoso(hcat(A...), threshold)
            return AJD.create_linear_filter(F)
        else
            throw(ArgumentError("Not supported for set of Matrices containing imaginary values!"))
        end
    end

    if algorithm == "jdiag_edourdpineau"
        F, _, _ = jdiag_edourdpineau(A)
        return AJD.create_linear_filter(F)
    end

    if algorithm == "FFD"
        _,F = FFD!(copy(A))
        return AJD.create_linear_filter(Matrix(F'))
    end
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

# Only export diagonalize().
# All provided functionality is available through the function.
export diagonalize, ajd_benchmark

end



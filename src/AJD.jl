module AJD 
using LinearAlgebra
using BenchmarkTools
using Plots

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
    threshold::AbstractFloat = eps(),
    plot_matrix::Bool = false,
    plot_convergence::Bool = false
    )

    if !check_input(A)
        throw(ArgumentError("Invalid input."))
    end

    if algorithm in ["jdiag", "jdiag_gabrieldernbach"]
        F, B, error_array = jdiag_gabrieldernbach!(A, max_iter = max_iter, threshold = threshold, plot_convergence = plot_convergence)
        #plotting(mean(B, dims=3)[:, :, 1])
        #plotting(Matrix(F))
        #@info A
        if plot_matrix
            plot_matrix_heatmap(F, B)
        end
        if plot_convergence
            plot_convergence_lineplot(error_array, algorithm)
        end
        return error_array, AJD.create_linear_filter(F)
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
        #plotting(F)
        return AJD.create_linear_filter(F)
    end

    if algorithm in ["FFD", "ffd", "ffdiag"]
        _,F = FFD!(copy(A))
        return AJD.create_linear_filter(Matrix(F'))
    end

    return throw(ArgumentError("No valid algorithm selected from available"))
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

"""
    plot_matrix_heatmap(filter::AbstractMatrix, diag_matrices)

Plot a heatmap of the calculated filter matrix and the mean of the diagonlized matrices.
"""
function plot_matrix_heatmap(filter::AbstractMatrix, diag_matrices) 
    theme(:dark)
    # Subplot of the filter matrix.
    filter_plot = heatmap(filter, yflip=true, title="Filter Matrix")

    # Subplot of the mean of all the diagonalised matrices.
    mean_diag_plot = heatmap(mean(diag_matrices, dims=3)[:, :, 1], yflip=true, title="Mean Diagonalized Matrices")

    # Combine subplots and 
    combined_plot = plot(filter_plot, mean_diag_plot, layout=(2, 1), legend=false, size=(400, 800))
    display(combined_plot)
end

"""
    plot_convergence_lineplot(error_array::AbstractArray, name::String)

Plot the convergence error as recorded during the algorithm execution.
"""
function plot_convergence_lineplot(error_array::AbstractArray, name::String)
    theme(:dark)
    line_plot = plot(error_array, w=3, title="Error Convergence", label=name)
    display(line_plot)
end

export diagonalize, ajd_benchmark

end



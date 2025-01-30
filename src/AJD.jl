module AJD
using LinearAlgebra: eigen, norm, Symmetric, Hermitian, I, qr, dot, diag
using ProgressMeter

# Import different algorithms.
include("jdiag_algorithms/jdiag_cardoso.jl")
include("jdiag_algorithms/jdiag_gabrieldernbach.jl")
include("jdiag_algorithms/jdiag_edourdpineau.jl")
include("FFDiag.jl")

# Utility functions, plotting functions and global constanst imported.
include("utils.jl")
include("utils_test_data.jl")
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

export diagonalize

end

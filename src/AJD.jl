module AJD
using LinearAlgebra: eigen, norm, Symmetric, Hermitian, I, qr, dot, diag
using ProgressMeter

abstract type AbstractDiagonalization end

# Algorithm types
struct JDiagGabrielDernbach <: AbstractDiagonalization end
struct JDiagEdourdPineau <: AbstractDiagonalization end
struct JDiagCardoso <: AbstractDiagonalization end
struct FFDiag <: AbstractDiagonalization end

# Traits
supportscomplex(::AbstractDiagonalization) = false
supportscomplex(::JDiagGabrielDernbach) = true
supportscomplex(::JDiagEdourdPineau) = true

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
(1)
    diagonalize(
        M::Vector{<:AbstractMatrix{<:Number}};
        algorithm::AbstractDiagonalization = JDiagEdourdPineau(),
        max_iter::Int = 1000,
        threshold::AbstractFloat = eps(),
    )

(2)
    diagonalize(
        M::Vector{<:AbstractMatrix{<:Number}},
        only_plot::Symbol;
        algorithm::AbstractDiagonalization = JDiagGabrielDernbach(),
        max_iter::Int = 1000,
        threshold::AbstractFloat = eps(),
    )

(3)
    diagonalize(
        benchmark::Symbol,
        n_dims::Int,
        n_matrices::Int,
    )

Calculate joint diagonalization of multiple input matrices ``M_k``.

Main function of the AJD package.
Implemented algorithms are JDiag and FFDiag.
Input of matrices ``M_k`` need to be a vector of matrices.
The matrices can be of types Float64 or Complex.

See the [Getting Started Guide](https://muehlefeldt.github.io/AJD.jl/dev/getting-started/) for information on the algorithms.

# Dispatch (1)
## Inputs
* `M`: Vector of matrices (requiered).
* `algorithm`: Selected algorithm from `JDiagGabrielDernbach()`, `JDiagEdourdPineau()`, `JDiagCardoso()` or `FFDiag()`.
* `max_iter`: Maximum iteration step as integer.
* `threshold`: Desired threshold minimizing the off-diagonal elements.
## Output
* Return LinearFilter object. Filter `fil` contains filter matrix `fil.F` and the inverse `fil.iF.`


# Dispatch (2)
## Inputs
* Additional symbol used to generate overview plot of the result. `Use diagonalize(M, :plot)`.
## Output
* Overview plot. Shows heatmap of the filter matrix, heatmap of the mean of the diagonalized matrices and the vonvergence behaviour of the algorithm.

# Dispatch (3)
## Inputs
* Symbol `:benchmark` used as `diagonalize(:benchmark, 10, 10)`.
* Automatic benchmark is run comparing JDiag and FFDiag algorithms. Using `n_dims` ``\\times`` `n_dims` matrices of count `n_matrices.` 
## Output
* BenchmarkGroup of package `BenchmarkTools` comparing the algorithms using diffrent types of test data.
"""
function diagonalize(
    M::Vector{<:AbstractMatrix{<:Number}};
    algorithm::AbstractDiagonalization = JDiagEdourdPineau(),
    max_iter::Int = 1000,
    threshold::AbstractFloat = eps(),
)
    check_input(M, max_iter, threshold)

    #convert integers to float in case 
    #input is of type Int

    if typeof(M) <: AbstractArray{<:AbstractArray{<:Int}}
        M = float.(M)
    end
    # Check complex support
    if !supportscomplex(algorithm) && any(x -> eltype(x) <: Complex, M)
        throw(ArgumentError("Selected algorithm doesn't support complex matrices"))
    end

    F, _, _, n_iter = get_diagonalization(
        M,
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


export AbstractDiagonalization,
    JDiagGabrielDernbach, JDiagEdourdPineau, JDiagCardoso, FFDiag
export diagonalize, get_diagonalization, supportscomplex
export ALL_ALGORITHMS, COMPLEX_ALGORITHMS

end

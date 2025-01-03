module AJD 
using LinearAlgebra

# Import different JDiag algorithms.
include("jdiag_algorithms/jdiag_cardoso.jl")
include("jdiag_algorithms/jdiag_gabrieldernbach.jl")
include("jdiag_algorithms/jdiag_edourdpineau.jl")

include("check_requirements.jl")
include("utils.jl")
"""
    diagonalize(
        A::Vector{<:AbstractMatrix{<:Union{Float64, ComplexF64}}};
        algorithm::String)

Diagonalize input matrices using the requested algorithm.

Main function of the AJD package.
Implemented algorithms at this point in time are limited to the [JDiag algorithm](https://doi.org/10.1137/S0895479893259546)  in different versions.
Input of matrices to be diagonalized need to be a vector of matrices.
The matrices can be of types Float64 or Complex. Limitations of the different implementations apply.

Supported algorithms are `jdiag_gabrieldernbach`, `jdiag_cardoso` and `jdiag_edourdpineau`.
See the Getting Started Guide for information on the algorithms.
"""

function diagonalize(
    A::Vector{<:AbstractMatrix{<:Number}};
    algorithm::String,
    max_iter::Int = 1000,
    threshold::AbstractFloat = eps()
    )
    if algorithm in ["jdiag", "jdiag_gabrieldernbach"]
        _, F = jdiag_gabrieldernbach!(A, max_iter = max_iter, threshold = threshold)
        return AJD.create_linear_filter(F)
    end

    if algorithm =="jdiag_cardoso"
        if typeof(A) <: AbstractArray{<:AbstractArray{<:Real}} 
            _, F ,_ = jdiag_cardoso(hcat(A...), 10e-8)
            return AJD.create_linear_filter(F)
        else
            throw(ArgumentError("Not supported for set of Matrices containing imaginary values!"))
        end
    end
    
    if algorithm == "jdiag_edourdpineau"
        F, _, _ = jdiag_edourdpineau(A)
        return AJD.create_linear_filter(F)
    end
    return error
end

export diagonalize

end



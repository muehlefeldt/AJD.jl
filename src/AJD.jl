module AJD 
using LinearAlgebra
include("jdiag_cardoso.jl")
include("jdiag_gabrieldernbach.jl")
include("jdiag_edourdpineau.jl")


function get_non_Diag_elements(A::AbstractMatrix)
    # Source: https://discourse.julialang.org/t/off-diagonal-elements-of-matrix/41169/4
    # Best method regarding compilation time found so far
        row, column = size(A)
        non_diag_elements_vector = [A[index_row, index_column] for index_row = 1:row, index_column = 1:column if index_row != index_column]
        return non_diag_elements_vector
end

function Jacobi_Rotation(G::Matrix)
    Eigenvalues, Eigenvector = eigen(G) #sorted by highest value last

    max_eigenvector = Eigenvector[:,end] #get the eigenvector of the corresponding highest eigenvalue
    #max_eigenvector = sign(max_eigenvector[1])*max_eigenvector #why is that? i don't know why i need to do that but the code says so?
    
    x = max_eigenvector[1]
    y = max_eigenvector[2]
    z = max_eigenvector[3]

    r = sqrt(x^2+y^2+z^2)

    c = sqrt((x+r)/2*r)

    s = (y - z*im)/(sqrt(2*r*(x+r)))
    R = [c conj(s); -s conj(c)]
    return R

end

"""
    diagonalize(A::Vector{<:AbstractMatrix{<:Union{Float64, ComplexF64}}}; algorithm::String)

Diagonalize input matrix using requested algorithm.

Main function of the AJD package.
Implemented algorithms at this point in time are limited to the algorithm in different versions.
Input of matrices to be diagonalized need to be a vector of matrices.
The matrices can be Float64 or complex. Limitations of the different algorithms apply.

Supported algorithms are "jdiag_gabrieldernbach", "jdiag_cardoso" and "jdiag_edourdpineau".
See the Getting Started Guide for information on the algorithms.
"""
function diagonalize(
    A::Vector{<:AbstractMatrix{<:Union{Float64, ComplexF64}}};
    algorithm::String
    )
    if algorithm in ["jdiag", "jdiag_gabrieldernbach"]
        return jdiag_gabrieldernbach(A)
    end
    if algorithm =="jdiag_cardoso"
        return jdiag_cardoso(hcat(A...), 10e-8)
    end
    if algorithm == "jdiag_edourdpineau"
        return jdiag_edourdpineau(A)
    end
    return error
end

export diagonalize


end



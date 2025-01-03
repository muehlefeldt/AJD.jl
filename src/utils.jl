using LinearAlgebra
using Diagonalizations

"""
    random_matrices(n::Int, m::Int)

Generate m random matrices of size n x n
"""
function random_matrices(n::Int, m::Int)
    return [rand(n,n) for i in 1:m]
end

"""
Generate m random symmetric matrices of size n x n
"""
function random_symmetric_matrices(n::Int, m::Int)
    return [Symmetric(rand(n,n)) for i in 1:m]
end

"""
Generate m random commuting matrices of size n x n
These will produce all real rotation matrices using the Jacobi method

M_i * M_j = M_j * M_i for all i,j
"""
function random_commuting_matrices(n::Int, m::Int)
    P = rand(n,n)
    return [P*Diagonal(rand(n))*inv(P) for _ in 1:m]
end

"""
    random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)
    
Generate m random normal commuting matrices of size n x n
These can be exactly diagonalized

M_i * M_j = M_j * M_i for all i,j
M_i*M_i' = M_i'*M_i for all i
"""
function random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)
    Q, _ = qr(rand(n,n))
    Q = Matrix(Q)
    if complex
        return [Q*Diagonal(rand(ComplexF64, n))*Q' for _ in 1:m]
    end
    return [Q*Diagonal(rand(n))*Q' for _ in 1:m]
end

"""
Create LinearFilter object as introduced by Diagonalizations.jl.
Output of AJD.jl follows convention of Diagonalizations.jl and produces a LinearFilter.
"""
function create_linear_filter(A::Matrix{T} where {T<:Number})
    args=("Approximate Joint Diagonalization", false)
    return LinearFilter(A, Matrix(A'), nothing, nothing, nothing, nothing, args...)
end
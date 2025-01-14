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
    return [P*Diagonal(rand(n))*inv(P) for _ in 1:m] #change inv to cholesky if possible?
end

"""
    random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)
    
Generate m random normal commuting matrices of size n x n
These can be exactly diagonalized

M_i * M_j = M_j * M_i for all i,j
M_i*M_i' = M_i'*M_i for all i
"""
function random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)
    # Just like the function below, this produces hermitian an symmetric matrices, they are
    # just not annotated as such
    # Q, _ = qr(rand(complex ? ComplexF64 : Float64, n,n))
    Q, _ = qr(rand(n,n))
    Q = Matrix(Q)
    if complex
        return [Q*Diagonal(rand(ComplexF64, n))*Q' for _ in 1:m]
    end
    return [Q*Diagonal(rand(n))*Q' for _ in 1:m]
end

function random_normal_commuting_symmetric_matrices(n::Int, m::Int; complex::Bool=false)
    Q, _ = qr(rand(complex ? ComplexF64 : Float64, n,n))
    Q = Matrix(Q)
    if complex
        return [Hermitian(Q*Diagonal(rand(n))*Q') for _ in 1:m]
    end
    return [Symmetric(Q*Diagonal(rand(n))*Q') for _ in 1:m]
end

function get_test_data_complex_real(n::Int, m::Int)
    Q, _ = qr(rand(n,n))
    Q = Matrix(Q)
    
    C = [Q * Diagonal(rand(ComplexF64, n)) * Q' for _ in 1:m]
    R = [Q * Diagonal(rand(n)) * Q' for _ in 1:m]
    return [C..., R...]
end

"""
Create LinearFilter object as introduced by Diagonalizations.jl.
Output of AJD.jl follows convention of Diagonalizations.jl and produces a LinearFilter.
"""
function create_linear_filter(A::Matrix{T} where {T<:Number})
    args=("Approximate Joint Diagonalization", false)
    return LinearFilter(A, Matrix(A'), nothing, nothing, nothing, nothing, args...)
end

"""
    frobenius_off_diag_norm(A::Array{<:Number,3})
Input
* A: Vector of matrices

Takes an array namely the Array of matrices A_k and gets the offdiagonal elements and applies the frobenius norm (∑ |a_{i,j}|^{2}). 
Used for the `jdiag_gabrieldernbach` and `FFD` algorithm.
"""
function frobenius_offdiag_norm(Xm::AbstractArray{T,3})::Real where {T<:Number}
    sum = zero(real(T))
    for i in axes(Xm, 1), j in axes(Xm, 2), k in axes(Xm, 3)
        i == j && continue
        sum += abs2(Xm[i, j, k])
    end
    return sum
end

"""
    get_offdiag_elements(A::Array{<:Number,3})
Input
* A: Vector of matrices

Takes an array of matrices, copies it and sets the diagonal to zero to get the matrix E with only the offdiagonal elements of A.

# Examples
```jldoctest
julia> A = ones(3,3,3);

julia> AJD.get_offdiag_elements(A)
3×3×3 Array{Float64, 3}:
[:, :, 1] =
 0.0  1.0  1.0
 1.0  0.0  1.0
 1.0  1.0  0.0

[:, :, 2] =
 0.0  1.0  1.0
 1.0  0.0  1.0
 1.0  1.0  0.0

[:, :, 3] =
 0.0  1.0  1.0
 1.0  0.0  1.0
 1.0  1.0  0.0```
"""
function get_offdiag_elements(A::Array{<:Number,3})
    rows,_ = size(A)
    #copy matrix A to not overwrite it with zeros
    E = copy(A)
    for row in 1:rows # maybe eachindex better - NG
        E[row,row,:] .= 0 
    end
    return E
end
"""
    get_diag_elements(A::Array)
Input
* A: Vector of matrices

Takes an array of matrices, takes the diagonal elements and returns the diagonal elements as a diagonal matrix D.

# Examples

```jldoctest
julia> A = ones(3,3,3);

julia> AJD.get_diag_elements(A)
3×3×3 Array{Float64, 3}:
[:, :, 1] =
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0

[:, :, 2] =
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0

[:, :, 3] =
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0```
"""
function get_diag_elements(A::Array)

    rows, columns, k = size(A)

    if typeof(A) <: AbstractArray{<:Complex}
        D = complex.(zeros(rows, columns, k))
    else
        D = zeros(rows, columns, k)
    end

    for row in 1:rows
        D[row,row,:] = A[row,row,:]
    end

    return D
end

function sort_offdiag_elements(A::AbstractArray{<:Number,3})
    rows,_,k = size(A)
    no_el = rows*rows-rows #number of offdiagonal elements
    if typeof(A) <: AbstractArray{<:Complex}
        sorted_array = complex.(zeros(no_el*k))
    else
        sorted_array = zeros(no_el*k)
    end

    for matrices in 1:k
        sorted_array[1+(matrices-1)*(no_el):matrices*(no_el)] = sort_offdiag_elements(A[:,:,matrices])
    end
    return sorted_array
end

function sort_offdiag_elements(A::AbstractArray{<:Number,2})
    rows,_ = size(A)

    if typeof(A) <: AbstractArray{<:Complex}
        sorted_array = complex.(zeros(rows*rows-rows))
    else
        sorted_array = zeros(rows*rows-rows)
    end
   
    i = 1
    for row in 1:rows-1,column in row+1:rows
        #TODO: Benchmark if append! or push! is faster
        sorted_array[i] = A[row, column]
        i += 1
        sorted_array[i] = A[column,row]
        i +=1
        
    end
    return sorted_array
end

"""
    is_commuting(A::AbstractMatrix, B::AbstractMatrix) 
Input:
* A: AbstractMatrix of dimension nxn
* B: AbstractMatrix of dimension nxn

Check if two matrices A, B are commuting.
 A * B = B * A must hold. 
 """
function is_commuting(A::AbstractMatrix, B::AbstractMatrix)
    return isapprox(A*B, B*A)
end
"""
    is_same_size(A::AbstractMatrix, B::AbstractMatrix)
Input:
* A: AbstractMatrix of variable size
* B: AbstractMatrix of variable size

Will return true if dimension of A and B is matching otherwise false.
"""
function is_same_size(A::AbstractMatrix, B::AbstractMatrix)
    return size(A) == size(B)
end
"""
    isstrictly_diagonally_dominant(A::AbstractMatrix)
Input:
* A: AbstractMatrix

Used for the FFDiag Algorithm to define whether the Matrix A is strictly diagonally dominant and therefore has an Inverse or not.
A matrix is strictly dominant if:```|a_{ii}| > \\sum |a_{ij}|, i ≠ j```
"""
function isstrictly_diagonally_dominant(A::AbstractMatrix)
    for i in eachindex(A[1:end, 1])
         
        if abs(sum(A[i,:])) - abs(A[i,i]) > abs(A[i,i]) ? true : false
            
            return false

        end
    end

    return true
end
"""
    get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)
Input
* D: Diagonal Matrix with offdiagonal elements set to zero
* i,j: Denotes the indexes the matrix D

Calculates the factor ```z_ij``` which is defined by:
``` ∑_{k} D_{i,i}^{k}D_{j,j}^{k} ```
"""
function get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)
    return sum(D[i,i,:].*D[j,j,:])
end
"""
    get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int, j::Int)
Input
* D: Diagonal Matrix with offdiagonal elements set to zero
* E: Diagonal Matrix with diagonal elements set to zero
* i,j: Denotes the indexes of the matrices D and E

Calculates the factor ```y_ij``` which is defined by:
```math ∑_{k} D_{j,j}^{k}E_{j,i}^{k} ```
"""
function get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int,j::Int)
    return sum(D[j,j,:].*E[i,j,:])
end

"Check for valid input of diagonalize()."
function check_input(A::Vector{<:AbstractMatrix{<:Number}})
    # Input may be empty.
    if length(A) <= 0
        return false
    end
    # All matrices must be commuting and of same size.
    for index in 1:length(A)-1
        if !is_same_size(A[index], A[index+1]) || !is_commuting(A[index], A[index+1])
            return false
        end
    end
    return true
end

function addrandomnoise(A::Vector{M};σ = 0.5,same_noise = true) where {T<:Number, M<:AbstractMatrix{T}}
    k = length(A)
    rows,columns = size(A[1])
    if same_noise == true
        R = randn(rows,columns)
        for index_k = 1:k
            A[index_k] = A[index_k] + σ*R*R'
        end
    else
        for index_k = 1:k
            R = randn(rows,columns)
            A[k] = A[k] + σ*R*R'
        end
    end
    return A
end
function addrandomnoise!(A::Vector{M};σ = 0.5,same_noise = true) where {T<:Number, M<:AbstractMatrix{T}}
    k = length(A)
    rows,columns = size(A[1])
    if same_noise == true
        R = randn(rows,columns)
        for index_k = 1:k
            A[index_k] = A[index_k] + σ*R*R'
        end
    else
        for index_k = 1:k
            R = randn(rows,columns)
            A[k] = A[k] + σ*R*R'
        end
    end
    return A
end
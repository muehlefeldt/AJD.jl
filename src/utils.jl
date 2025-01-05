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

"""
    frobenius_off_diag_normation(A::Array)
Input
* A: Vector of communting matrices with index k
Takes an array namely the Array of matrices A_k and gets the offdiagonal elements and applies the frobenius norm (∑ |a_{i,j}|^{2}). 
Used for the `jdiag_gabrieldernbach` and `FFD` algorithm.
"""
#TODO: Might be better to put back into jdiag_gabrieldernbach since it is only used there? - NG
function frobenius_offdiag_normation(A::Array{<:Number,3})
    #slower method instead of using get_off_diag_elements
    # row,column,k = size(A)
    # non_diag_elements_vector = [A[index_row, index_column,index_k] for index_row = 1:row, index_column = 1:column, index_k = 1:k if index_row != index_column]
    #new method -NG
    return sum(abs.(get_offdiag_elements(A)).^2) #frobenius norm is applied
end
#TODO: Find out if needed
# function frobenius_off_diag_normation(A::AbstractArray{<:Number,2})
#     row,column = size(A)
#     non_diag_elements_vector = [A[index_row, index_column] for index_row = 1:row, index_column = 1:column if index_row != index_column]
#     return sum(abs.(non_diag_elements_vector).^2) #frobenius norm is applied
# end

function get_offdiag_elements(A::Array{<:Number,3})
    rows,_ = size(A)
    E = copy(A)
    for row in 1:rows # maybe eachindex better - NG
        E[row,row,:] .= 0 
    end
    return E
end

function get_diag_elements(A::Array)

    rows, columns, k = size(A)
    D = zeros(rows, columns, k)

    for row in 1:rows
        D[row,row,:] = A[row,row,:]
    end

    return D
end

function sort_offdiag_elements(A::AbstractArray{<:Number,3})
    rows,_,k = size(A)
    no_el = rows*rows-rows #number of offdiagonal elements
    sorted_array = zeros((no_el)*k)

    for matrices in 1:k
        sorted_array[1+(matrices-1)*(no_el):matrices*(no_el)] = sort_offdiag_elements(A[:,:,matrices])
    end
    return sorted_array
end

function sort_offdiag_elements(A::AbstractArray{<:Number,2})
    rows,_ = size(A)
    sorted_array = zeros(rows*rows-rows)
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

#TODO: Unclear if needed for FFDiag
# function inf_off_diag_normation(W::Array{<:Number})
#     row,column = size(W)
#     objective = 0
#     i = 1
#     for row in eachrow(W)
#         objective_new = sum(abs.(row)) - abs(W[i,i])
#         i += 1
#         if objective_new > objective
#             objective = objective_new
#         end
        
#     end
#     return objective
# end

# Check if two matrices A, B are commuting.
# A * B = B * A must hold.
function is_commuting(A::AbstractMatrix, B::AbstractMatrix)
    return isapprox(A*B, B*A)
end

function is_same_size(A::AbstractMatrix, B::AbstractMatrix)
    return size(A) == size(B)
end

function isstrictly_diagonally_dominant(A::AbstractMatrix)
    for i in eachindex(A[1:end, 1])
         
        if abs(sum(A[i,:])) - abs(A[i,i]) > abs(A[i,i]) ? true : false
            
            return false

        end
    end

    return true
end

function get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)
    return sum(D[i,i,:].*D[j,j,:])
end
function get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int,j::Int)
    return sum(D[j,j,:].*E[i,j,:])
end
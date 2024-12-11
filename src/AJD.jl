module AJD

"""
    multiply(x, y)

Multiply inputs `x` and `y`.

Function only introduced to test module function against Python implementation.
"""
function multiply(x, y)
    return x * y
end

export multiply

using LinearAlgebra
# Write your package code here.
function Is_Commuting(A::AbstractMatrix, B::AbstractMatrix)
    return A*B == B*A
end

function Is_Same_size(A::AbstractMatrix, B::AbstractMatrix)
    return size(A) == size(B)
end

function Is_Symmetric(A::AbstractMatrix)
    return size(A,1) == size(A,2)
end

function get_non_Diag_elements(A::AbstractMatrix)
    #Code found under: https://discourse.julialang.org/t/off-diagonal-elements-of-matrix/41169/4
    #best method regarding compilation time found so far
        row, column = size(A)
        non_diag_elements_vector = [A[index_row, index_column] for index_row = 1:row, index_column = 1:column if index_row != index_column]
        return non_diag_elements_vector
end
function Jacobi_Rotation(G::Matrix)

    Eigenvalues,Eigenvector = eigen(G) #sorted by highest value last

    max_eigenvector = Eigenvector[:,end] #get the eigenvector of the corresponding highest eigenvalue
    max_eigenvector = sign(max_eigenvector[1])*max_eigenvector #why is that? i don't know why i need to do that but the code says so?
    
    x = max_eigenvector[1]
    y = max_eigenvector[2]
    z = max_eigenvector[3]

    r = sqrt(x^2+y^2+z^2)

    c = sqrt((x+r)/2*r)

    s = (y - z*im)/(sqrt(2*r(x+r)))
    R = [c conj(s); -s conj(c)]
    return R

end
function JADE(A::AbstractArray;threshhold = 10e-18, max_iter = 1000)
    #A concatenate in third dimension by  A =[[1 2; 1 2];;;[2 3; 4 5]]
    
    rows, columns, k = size(A)

    #initialize the apporximate joint eigenvecotrs as described in Cardoso
    V = (1.0+0.0*im)*I(rows)
    iteration_step = 0
    #Calculate h first for all the entries in the matrix
    while iteration_step <= max_iter
        for row = 1:rows
            for column = 2:columns
            end
        end 

    end

    

    
    return  V
end

export Is_Commuting
export Is_Same_size
export Is_Symmetric
export get_non_Diag_elements
export Jacobi_Rotation
export JADE
end

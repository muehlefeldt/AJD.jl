module AJD

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
function get_Diag_elements(A::AbstractMatrix)
    #Code found under: https://discourse.julialang.org/t/off-diagonal-elements-of-matrix/41169/4
    #best method regarding compilation time found so far
        row, column = size(ma)
        [ma[index_row, index_column] for index_row = 1:row, index_column = 1:column if index_row != index_column]
      end 
end
export Is_Commuting
export Is_Same_size
export Is_Symmetric
export get_Diag_elements
end

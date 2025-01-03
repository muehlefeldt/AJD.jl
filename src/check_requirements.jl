function iscommuting(A::AbstractMatrix, B::AbstractMatrix)
    return A*B == B*A
end

function issamesize(A::AbstractMatrix, B::AbstractMatrix)
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
# function issymmetric(A::AbstractMatrix)
#     return size(A,1) == size(A,2)
# end
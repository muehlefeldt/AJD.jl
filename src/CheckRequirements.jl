module CheckRequirements
function iscommuting(A::AbstractMatrix, B::AbstractMatrix)
    return A*B == B*A
end

function issamesize(A::AbstractMatrix, B::AbstractMatrix)
    return size(A) == size(B)
end

# function issymmetric(A::AbstractMatrix)
#     return size(A,1) == size(A,2)
# end

export iscommuting, issamesize
# , issymmetric

end
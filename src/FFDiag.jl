include("check_requirements.jl")
using LinearAlgebra
using PosDefManifold

function FFD!(A::Vector{M}) where {T <: Number, M<:AbstractMatrix{T}}
    
    if typeof(A) <: AbstractArray{<:Int}
        A = float.(A)
    end
    
    A = cat(A..., dims = 3)
    rows,columns,k = size(A)
    for row in 1:rows
        A[row,row] 
    end

end
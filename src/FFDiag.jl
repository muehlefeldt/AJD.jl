using LinearAlgebra
using PosDefManifold
using SparseArrays

function FFD!(A::Vector{M}; threshold = eps(), max_iter = 1000, norm = "frobenius",θ = 1) where {T <: Number, M<:AbstractMatrix{T}}
    
    if typeof(A) <: AbstractArray{<:AbstractArray{<:Int}}
        A = float.(A)
    end
    if norm == "frobenius"
        normation = X -> norm(X,2) #frobenius norm
    elseif norm == "inf"
        normation = X -> opnorm(X,Inf) #infinity norm
    end

    A = cat(A..., dims = 3)
    rows,columns,k = size(A)
    #initialization
    iteration_step = 0
    active = true
    V = 1.0*Matrix(I, rows,columns)
    W = zeros(rows,columns)
    
    E = get_offdiag_elements(A)
    D = get_diag_elements(A)
    
    w = sort_offdiag_elements(W)
    e = sort_offdiag_elements(E)
    objective = frobenius_offdiag_normation(W.*D+D.*W'+E)
    while active == true && iteration_step <= max_iter
        W = 
        if normation(W) > θ
            W = θ/normation(W)*W
        end
        A = D + W.*D + D.*W'+E
        E = get_offdiag_elements(A)
        D = get_diag_elements(A)
        
        
        objective_new = frobenius_offdiag_normation(W.*D+D.*W'+E)
        diff = objective - objective_new
        objective = objective_new #what about pointer?
        
        if abs(diff) > threshold
            active = true
        else
            active = false
        end
        iteration_step += 1
    end
    

end

function off_diag_sorted()

end
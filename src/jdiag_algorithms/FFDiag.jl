using LinearAlgebra
using PosDefManifold

function FFD!(A::Vector{M}; threshold = eps(), max_iter = 1000, norm_ = "frobenius",θ = 1) where {T <: Number, M<:AbstractMatrix{T}}
    
    if typeof(A) <: AbstractArray{<:AbstractArray{<:Int}}
        A = float.(A)
    end
    if norm_ == "frobenius"
        normation = X -> norm(X,2) #frobenius norm, it says it is frobenius norm but that doesn't make sense!
    elseif norm_ == "inf"
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
    objective = frobenius_offdiag_normation(A)
    
    #while active == true && iteration_step <= max_iter
    while iteration_step <= max_iter
        #TODO: This needs a more elegant solution
        for i = 1:rows-1, j = i+1:rows
            z_ij = get_z_fdiag(D,i,j)
            z_i = get_z_fdiag(D,i,i)
            z_j = get_z_fdiag(D,j,j)
            y_ij = get_y_fdiag(D,E,i,j)
            y_ji = get_y_fdiag(D,E,j,i)
            W[i,j] = (z_ij*y_ji - z_i*y_ij)/(z_j*z_i-z_ij^2)
            W[j,i] = (z_ij*y_ij - z_j*y_ji)/(z_j*z_i-z_ij^2)
        end
        
        if normation(W) > θ
            W = θ/normation(W)*W
        end

        A = D + W.*D + D.*W'+E
        V = (I+W)*V
        #A = (I+W).*A.*(I+W)'      
        
        objective_new = frobenius_offdiag_normation(A)
        diff = objective - objective_new
        objective = copy(objective_new) #what about pointer?
        E = get_offdiag_elements(A)
        D = get_diag_elements(A)
        if abs(diff) > threshold
            active = true
        else
            active = false
        end
        iteration_step += 1
    end
    @info iteration_step
    return A,V
end
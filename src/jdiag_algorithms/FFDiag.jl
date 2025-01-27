using LinearAlgebra
using PosDefManifold
"""
    function ffd(A::Vector{M}; <keyword_arguments>) where {T <: Number, M<:AbstractMatrix{T}}
*`A`: Vector of matrices of dimension ``n × n × k``

*`threshold`: threshold until calculation should stop. default is eps().

*`max_iter`: max number of iterations. default is 100.

*`norm_`: norm by which the update matrix is divided by. can either be "frobenius" or "inf". default is "frobenius"

*`θ`: normation criterion. should be smaller than 1. default is 0.99. See paper on [ffdiag](https://www.jmlr.org/papers/volume5/ziehe04a/ziehe04a.pdf) for more information.

*`plot_convergence`: whether the convergence plot should be shown or not. default is false.

*`initial_guess`: is the matrix by which the matrices can be diagonalized. if close to the solution the calculation gets faster. default is identity matrix with size ``n × n``.

Calculates the diagonalization of a set of matrices proposed in [2](https://www.jmlr.org/papers/volume5/ziehe04a/ziehe04a.pdf). 
Might be faster than the jade algorithms for certain matrices if initial guess is close to the diagonalized solution.

Won't work if input `A` is of size ``n × n × 1 `` since this will lead to NaN values due to how the update matrix is calculated or if all of the entries of `A` are the same.
Will lead to a warning and stops the calculation.
"""
function ffd(
    A::Vector{M};
    threshold = eps(),
    max_iter = 100,
    norm_ = "frobenius",
    θ = 0.99,
    plot_convergence::Bool = false,
    initial_guess = 1.0*Matrix(I,size(A[1])[1],size(A[1])[2])) where {T <: Number, M<:AbstractMatrix{T}}
    
    if typeof(A) <: AbstractArray{<:AbstractArray{<:Int}}
        A = float.(A)
    end
    if norm_ == "frobenius"
        norm_function = X -> norm(X,2) #frobenius norm, it says it is frobenius norm but that doesn't make sense!
    elseif norm_ == "inf"
        norm_function = X -> opnorm(X,Inf) #infinity norm
    end

    A = cat(A..., dims = 3)
    rows,columns,k = size(A)
    #initialization
    iteration_step = 0
    
    V = initial_guess
    W = zeros(rows,columns)
        
    E = get_offdiag_elements(A)
    D = get_diag_elements(A)
    
    objective = frobenius_offdiag_norm(A)
    
    error_array = Float64[] 
    if plot_convergence
        push!(error_array, frobenius_offdiag_norm(A))
    end

   
    while iteration_step <= max_iter
       
        for i = 1:rows-1, j = i+1:rows
            z_ij = get_z_fdiag(D,i,j)
            z_i = get_z_fdiag(D,i,i)
            z_j = get_z_fdiag(D,j,j)
            y_ij = get_y_fdiag(D,E,i,j)
            y_ji = get_y_fdiag(D,E,j,i)
           
            W[i,j] = (z_ij*y_ji - z_i*y_ij)/(z_j*z_i-z_ij^2) 
            W[j,i] = (z_ij*y_ij - z_j*y_ji)/(z_j*z_i-z_ij^2)
        end
        
        if norm_function(W) > θ
            W = θ/norm_function(W)*W
        end
        
        V = (I+W)*V
        for m = 1:k
            A[:,:,m] = (I+W)*A[:,:,m]*(I+W)'
        end
        
        objective_new = frobenius_offdiag_norm(A)
        diff = objective - objective_new
        if !isfinite(diff)
            @warn ("There are NaN values in the matrix. This can occur for a set of a single matrix or if all matrices are the same.")
            break
        end
        abs(diff) < threshold && break

        objective = objective_new 
        
        E = get_offdiag_elements(A)
        D = get_diag_elements(A)
        
        iteration_step += 1

        if plot_convergence
            push!(error_array, frobenius_offdiag_norm(A))
        end
    end
    return Matrix(V'), A, error_array
end
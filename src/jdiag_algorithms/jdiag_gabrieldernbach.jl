"""
    (1) jdiag_gabrieldernbach!(
        A::Vector{M};
        threshold::AbstractFloat = eps(),
        rel_threshold = 1e-3, #not used in function but easier for multiple dispatch
        max_iter = 1000,
        plot_convergence::Bool = false) where {T<:Real, M<:AbstractMatrix{T}}

    (2) jdiag_gabrieldernbach!(
    A::Vector{M};
    threshold = eps(),
    rel_threshold = 1e-3,
    max_iter = 1000,
    plot_convergence::Bool = false) where {T<:Complex, M<:AbstractMatrix{T}}

* `A`: set of matrices with dimension `` n × n × k``
* `threshhold`: absolute error used for stopping the iteration
* `rel_threshold`: relative error used for stopping the iteration

Computes the JDiag Algorithm for an Input of a Vector of Matrices and returns the diagonalization matrix, 
the diagonalized set of matrices and corresponding algortihm parameters like iteration step and corresponding error from frobenius norm calculation.

*Dispatch (1)*
JDiag algorithm based on the implementation by Gabrieldernbach in Python.
Source: [Algorithm](https://github.com/gabrieldernbach/approximate_joint_diagonalization/blob/master/jade/jade_cpu.py)

*Dispatch (2)*
JDiag algorithm for complex matrices based on the implementation by Gabrieldernbach in Python, the Cardoso Paper and the code 
of [Algorithm](https://github.com/edouardpineau/Time-Series-ICA-with-SOBI-Jacobi)

"""
function jdiag_gabrieldernbach!(
        A::Vector{M};
        threshold::AbstractFloat = eps(),
        rel_threshold = 1e-3, #not used in function but easier for multiple dispatch
        max_iter = 1000,
        plot_convergence::Bool = false) where {T<:Real, M<:AbstractMatrix{T}}
    
    
    A = cat(A...,dims = 3)::AbstractArray{<:Real}
    #convert to 3 dimensional matrix and concatenate in the third dimension
    #will also reset dimensions of matrix if OffsetArray
    rows, columns, k = size(A)

    error_array = Float64[] 
    if plot_convergence
        push!(error_array, frobenius_offdiag_norm(A))
    end

    # Initial setup of the progressbar.
    diff = threshold
    progress_bar = ProgressThresh(diff; desc="Minimizing:")

    #initialize the approximate joint eigenvecotrs as described in Cardoso
    V = Matrix((1.0)*I(rows)) 
    
    iteration_step = 0

    active = true #flag if threshold is reached

    while iteration_step < max_iter && active == true
        iteration_step += 1
        active = false

        for row = 1:rows
            for column = row+1:columns #row_index != column_index
                h_diag = A[row,row,:] - A[column,column,:] #first entry of h in Cardoso paper
                h_non_diag = A[row,column,:] + A[column,row,:] #second entry of h in Cardoso paper
                
                #rotational computations
                ton = dot(h_diag,h_diag) - dot(h_non_diag,h_non_diag)
                toff = 2*dot(h_diag,h_non_diag)
                θ = 0.5*atan(toff, ton + sqrt(ton*ton + toff * toff))

                #rotational arguments for real matrices
                c = cos(θ)
                s = sin(θ)
                diff = abs(s)
                R = [ c s; -s c]

                # if threshold for minimize function is reached abort calculations 
                #(or max iter is reached)
                # otherwise rotation is applied to matrix
                active = diff > threshold
                
                if diff > threshold

                    pair = [row, column]
                
                    for n = 1:k
                        A[:,pair,n] = transpose(R*transpose(@view(A[:,pair,n])))
                        A[pair,:,n] = R*@view(A[pair,:,n])
                    end

                    V[:,pair] = transpose(R*transpose(@view(V[:,pair])))

                end
                
            end
        end 

        # Update progress info.
        update!(progress_bar, diff)

        if plot_convergence
            push!(error_array, frobenius_offdiag_norm(A))
        end

        
    end
    
    # Return of the filter, diagonalized matrices and the convergence error (optional).
    return V, A, error_array, iteration_step

end

function jdiag_gabrieldernbach!(
    A::Vector{M};
    threshold = eps(),
    rel_threshold = 1e-3,
    max_iter = 1000,
    plot_convergence::Bool = false) where {T<:Complex, M<:AbstractMatrix{T}}
    
    #fixes type instability of cat
    #with the Output type
    A = cat(A...,dims = 3)::AbstractArray{<:Complex}
    rows, columns, k = size(A)
    #initialize the apporximate joint eigenvecotrs as described in Cardoso
    V = complex.(Matrix((1.0)*I(rows))) 
    #needs to be added otherwise we cannot manipulate the non diag. elements of V
    
    #objective_function to be minimized by algorithm
    objective = frobenius_offdiag_norm(A)
   
    #conditions for abortion initialized
    iteration_step = 0
    active = true

    # Make sure empty error array exists even if not tracked.
    # Track error if plot_convergence is selected.
    error_array = Float64[] 
    if plot_convergence
        push!(error_array, frobenius_offdiag_norm(A))
    end

    # Initial setup of the progressbar.
    progress_bar = ProgressThresh(threshold; desc="Minimizing:")

    while iteration_step < max_iter && active == true
        iteration_step += 1
        active = false
        #calculations of cardoso paper
        for row = 1:rows
        
            for column = row+1:columns
               
                h_diag = A[row,row,:] - A[column, column,:]
                h_non_diag = A[row,column,:] + A[column,row,:]
                h_imag = im*(A[column,row,:] - A[row,column,:])
                #changed the h_imag to be imaginary
            
                h = [h_diag h_non_diag h_imag]
                
                #initialize the matrix G consisting 
                #of h as mentioned in Cardoso [1]
                G = Matrix{Number}[]
               
                for k_index = 1:k
                    if isempty(G) == false
                        G = G + real(adjoint(transpose(@view(h[k,:])))*
                        transpose(@view(h[k,:])))
                    else
                        G = real(adjoint(transpose(@view(h[k,:])))*
                        transpose(@view(h[k,:])))
                    end
                end
                
                R = Jacobi_Rotation(G)
                pair = [row,column]
                #A[:,pair,n] = transpose(R*transpose(A[:,pair,n]))
                #A[pair,:,n] = R*A[pair,:,n]
                for k_index = 1:k                                       
                    A[:,pair,k_index] = @view(A[:,pair,k_index])*R'
                    A[pair,:,k_index] = R*@view(A[pair,:,k_index])
                end
                
                V[:,[row,column]] = @view(V[:,[row,column]])*R'
            end
        
        end
    
        objective_new = frobenius_offdiag_norm(A)
        diff = objective_new - objective

        # Update progress info.
        update!(progress_bar, diff)

        if abs(diff) > threshold || rel_threshold*objective > diff
            active = true
        end

        objective= objective_new
        

        # Add error at the end of the iteration to track error convergence.
        if plot_convergence
            push!(error_array, objective)
        end

    end

    # Return of the filter, diagonalized matrices and the convergence error (optional).
    return V, A, error_array, iteration_step

end

function Jacobi_Rotation(G::Matrix)
    Eigenvalues, Eigenvector = eigen(G) #sorted by highest value last

    max_eigenvector = @view(Eigenvector[:,end]) #get the eigenvector of the corresponding highest eigenvalue
    
    x,y,z = max_eigenvector
    if x < 0.0
        x, y, z = -x, -y, -z
    end
    r = sqrt(abs(x^2+y^2+z^2))

    c = sqrt(abs((x+r)/2*r))

    s = (y - z*im)/(sqrt(abs(2*r*(x+r))))
    R = [c conj(s); -s conj(c)]
    return R

end
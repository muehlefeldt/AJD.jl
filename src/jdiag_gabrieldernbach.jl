"""
    jdiag_gabrieldernbach(A::Vector{Matrix{Float64}}; threshold = 10e-18, max_iter = 1000)

JDiag algorithm based on the implementation by Gabrieldernbach in Python.

Source: https://github.com/gabrieldernbach/approximate_joint_diagonalization/blob/master/jade/jade_cpu.py
"""
function jdiag_gabrieldernbach(A::Vector{Matrix{Float64}}; threshold = eps(), max_iter = 1000)
    #A concatenate in third dimension by  A =[[1 2; 1 2];;;[2 3; 4 5]]
    #only works for Real Matrices of A but not complex
    #A = Float64.(A) #if the Array isn't already of Float64
    A = cat(A...,dims = 3)
    rows, columns, k = size(A)

    #initialize the apporximate joint eigenvecotrs as described in Cardoso
    V = Matrix((1.0)*I(rows)) #needs to be added otherwise we cannot manipulate the non diag. elements of V
    
    iteration_step = 0

    active = true #flag if threshold is reached
    while iteration_step <= max_iter && active == true
   
        active = false

        for row = 1:rows
            for column = row+1:columns #row_index != column_index
                h_diag = A[row,row,:] - A[column,column,:] #first entry of h
                h_non_diag = A[row,column,:] + A[column,row,:] #second entry of h
                
                ton = dot(h_diag,h_diag) - dot(h_non_diag,h_non_diag)
                toff = 2*dot(h_diag,h_non_diag)
                θ = 0.5*atan(toff, ton + sqrt(ton*ton + toff * toff))

                c = cos(θ)
                s = sin(θ)
                R = [ c s; -s c]
                active = active || abs(s) > threshold
                if abs(s) > threshold
                    pair = [row, column]
                
                    for n = 1:k
                        A[:,pair,n] = transpose(R*transpose(A[:,pair,n]))
                        A[pair,:,n] = R*A[pair,:,n]
                        
                    end
                    V[:,pair] = transpose(R*transpose(V[:,pair]))
                end
                
            end
        end 
        iteration_step += 1
    end
    
    return  A,V

end
function off_diag_normation(A::Array)
    row, column,k = size(A)
    non_diag_elements_vector = [A[index_row, index_column,index_k] for index_row = 1:row, index_column = 1:column, index_k = 1:k if index_row != index_column]
    normation = sum(abs.(non_diag_elements_vector).^2)

    return normation
end

function jdiag_gabrieldernbach(A::Vector{Matrix{ComplexF64}}; threshold = eps(), max_iter = 1000)

    A = cat(A...,dims = 3)
    rows, columns, k = size(A)
    #initialize the apporximate joint eigenvecotrs as described in Cardoso
    V = Matrix((1.0)*I(rows)+im*zeros(rows,rows)) #needs to be added otherwise we cannot manipulate the non diag. elements of V

    objective_function = off_diag_normation(A)
   
    iteration_step = 0
    active = true
    


    while iteration_step <= max_iter && active == true
        active = false
        for row = 1:rows
            for column = row+1:columns
                #TODO: throws a NaN for the last values in the Matrix A, unclear why it does that, works well with normal indexing!
                h_diag = A[row,row,:] - A[column, column,:]
                h_non_diag = A[row,column,:] + A[column,row,:]
                h_imag = A[column,row,:] - A[row,column,:]
            
                h = [h_diag h_non_diag h_imag]
                
                #TODO: Make h_diag elements the transposed vecotrs!
                G = Matrix{Number}[]
                for k_index = 1:k
                    if isempty(G) == false
                        G = G + real(adjoint(transpose(h[k,:]))*transpose(h[k,:]))
                    else
                        G = real(adjoint(transpose(h[k,:]))*transpose(h[k,:]))
                    end
                end
                #@info typeof(G)
                R = Jacobi_Rotation(G)
                
                for k_index = 1:k
                    A[[row,column],[row,column],k_index] = R*A[[row,column],[row,column],k_index]*adjoint(R) #might not be correct, maybe use the matrix and multiply like in the python code
                end
            end
        
        end
    
        objective_function_new = off_diag_normation(A)
        diff = objective_function_new - objective_function

        if abs(diff) <= threshold
            active = true
        end
        objective_function = objective_function_new
        iteration_step += 1
    end
    return A
    #TODO: Return Eigenvectors V as well
end

function Is_Commuting(A::AbstractMatrix, B::AbstractMatrix)
    return A*B == B*A
end

function Is_Same_size(A::AbstractMatrix, B::AbstractMatrix)
    return size(A) == size(B)
end

function Is_Symmetric(A::AbstractMatrix)
    return size(A,1) == size(A,2)
end
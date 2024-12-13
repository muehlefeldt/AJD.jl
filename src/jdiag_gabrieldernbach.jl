"""
    jdiag_gabrieldernbach(A::Vector{Matrix{Float64}}; threshold = 10e-18, max_iter = 1000)

JDiag algorithm based on the implementation by Gabrieldernbach in Python.

Source: https://github.com/gabrieldernbach/approximate_joint_diagonalization/blob/master/jade/jade_cpu.py
"""

function jdiag_gabrieldernbach(A::Vector{Matrix{Float64}}; threshold = 10e-18, max_iter = 1000)
    #A concatenate in third dimension by  A =[[1 2; 1 2];;;[2 3; 4 5]]
    #only works for Real Matrices of A but not complex
    #A = Float64.(A) #if the Array isn't already of Float64
    A = cat(A...,dims = 3)
    rows, columns, k = size(A)

    #initialize the apporximate joint eigenvecotrs as described in Cardoso
    V = Matrix((1.0)*I(rows)) #needs to be added otherwise we cannot manipulate the non diag. elements of V
    
    iteration_step = 0

    active = true #flag if threshold is reached
    while iteration_step >= max_iter || active == true
   
        active = false

        for row = 1:rows
            for column = 2:columns
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

function Is_Commuting(A::AbstractMatrix, B::AbstractMatrix)
    return A*B == B*A
end

function Is_Same_size(A::AbstractMatrix, B::AbstractMatrix)
    return size(A) == size(B)
end

function Is_Symmetric(A::AbstractMatrix)
    return size(A,1) == size(A,2)
end
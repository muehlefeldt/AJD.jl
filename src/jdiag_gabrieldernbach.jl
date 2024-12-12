function JADE(A::AbstractArray, threshold = 10e-18, max_iter = 1000)
    #A concatenate in third dimension by  A =[[1 2; 1 2];;;[2 3; 4 5]]
    #only works for Real Matrices of A but not complex
    A = Float64.(A) #if the Array isn't already of Float64
    
    rows, columns, k = size(A)

    #initialize the apporximate joint eigenvecotrs as described in Cardoso
    V = (1.0)*I(rows)+zeros(rows, columns) #needs to be added otherwise we cannot manipulate the non diag. elements of V
    
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
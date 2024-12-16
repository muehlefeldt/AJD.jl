module AJD


    """
        jdiag_cardoso(A::Vector{Matrix{Float64}}; threshold = 10e-8)

    JDiag algorithm based on the implementation by Jean-Fran\c{c}ois Cardoso in Matlab.

    Source: https://www2.iap.fr/users/cardoso/code/Joint_Diag/joint_diag.m
    """

    function jdiag_cardoso(M::Vector{Matrix},jthresh=10e-8)
        # This version only works for matrix with real valued entries
        # Input:
        # M is a vector with n mxm matrices [M1,...,Mn]
        # thresh is a threshold for approximation stop, default = 10e-8
        # Output:
        # V : is a  mxm matrix, which accumulates givens rotations G in each iteration
        # A : is a mxnm matrix, which contains [VM1V',...,VMnV']
        # iter: accumulates the iteration numbers


		
        A = concatenate_to_one_matrix(M)
        m,nm = size(A)
        iter = 0

        V= Matrix{Float64}(I,m,m)

        # Initializing flag
        flag = true
        B = [1 0 0;0 1 1;0 -im im]
        

        while flag
            flag = false

            iter+=1
            for p in 1: m-1
                iter+=1
                Ip = p:m:nm
                
                for q in p+1:m
                    iter+=1
                    Iq = q:m:nm



                    # computing the givens angles base on Cardoso's paper
                    g = [(A[p,Ip]-A[q,Iq])';A[p,Iq]';A[q,Ip]']
                    
                    # Calculating eigenvalue and eigenvector of 
                    F = eigen(real(B*(g*g')B'))
                    sp = sortperm(F.values)
                    # choose the eigenvector with largest eigenvalue
                    angles = F.vectors[:,sp[3]]
                    if angles[1] < 0
                        angles =-angles
                    end
                    @assert angles[3] ==0
                    # calculate the parameters for givens rotation
                    c = sqrt(0.5+angles[1]/2)
                    s = 0.5*(angles[2])/c
                    
                    #update matrices A and V by a givens rotation
                    if abs(s)>jthresh

                        flag = true
                        pair = [p,q]
                        G = [c -conj(s); s c]

                        # update V, which accumulates givens rotations
                        V[:,pair] = V[:,pair]*G



                        
                        # update related two rows, p,q, of Real matrix A by one givens rotation 
                        A[pair,:] = G'*A[pair,:]

                        




                        # update related two columns ,p,q of Real matrix A by one givens rotation
                        A[:,Ip] = c*A[:,Ip]+s*A[:,Iq] 

                        A[:,Iq] = -conj(s)*A[:,Ip]+c*A[:,Iq]

                        
                    end

                end
            end
        end
        return A,V,iter
    end



    function concatenate_to_one_matrix(M::Vector{Matrix})
		return hcat(M)
	end




end



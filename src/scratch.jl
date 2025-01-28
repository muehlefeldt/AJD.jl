function off_diagonal_norm_cardoso(A::AbstractMatrix)
    # Input: m x nm matrix
    #        input has to be n mxm matrices
    # Output: the sum of the euklidian norm of these n matrices with diagonal elements are 0

    m,nm = size(A)
    n = nm ÷ m
    
    norm_sum = 0.0
    # @show n,m
    for i in 1:n
        M = A[:,(i-1)*m + 1:i*m]

        off_diag = M .-Diagonal(diag(M))

        norm_sum += norm(off_diag)^2
    end

    return sqrt(norm_sum)
end


function jdiag_cardoso(N::Vector{M},jthresh=10e-8,max_iter=800) where {M<:AbstractMatrix}
    # This version  works for matrix with both real and complex valued entries
    # Input:
    # M is a vector with n mxm matrices [M1,...,Mn]
    # thresh is a threshold for approximation stop, default = 10e-8
    # Output:
    #       iter: accumulates the iteration numbers
    #       V : is a  mxm matrix, which accumulates givens rotations G in each iteration
    #       A : is a mxnm matrix, which contains [VM1V',...,VMnV']
    #       off_norm_array : a vector which contains every off_norm after each iter

    
    A = ComplexF64.(hcat(N))
    m,nm = size(A)
    iter = 0
    off_norm = off_diagonal_norm_cardoso(A)
    off_norm_array = [off_norm]
    V= Matrix{ComplexF64}(I,m,m)

    # Initializing flag
    flag = true
    B = [1 0 0;0 1 1;0 -im im]

    while flag && iter < max_iter
        flag = false

        iter+=1
        for p in 1: m-1
            Ip = p:m:nm
            
            for q in p+1:m
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
                # @assert angles[3] ==0
                # calculate the parameters for givens rotation
                c = sqrt((1+angles[1])/2)
                s = (angles[2]-angles[3]*im)/sqrt(2.0*(angles[1]+1.0))
                
                #update matrices A and V by a givens rotation
                if off_norm > jthresh

                    flag = true
                    pair = [p,q]
                    G = [c -conj(s); s c]

                    # update V, which accumulates givens rotations
                    V[:,pair] = V[:,pair]*G
                    
                    # update related two rows, p,q, of Real matrix A by one givens rotation 
                    A[pair,:] = G'*A[pair,:]

                    # update related two columns ,p,q of Real matrix A by one givens rotation
                    A[:,[Ip Iq]] = [c*A[:,Ip]+s*A[:,Iq];;-conj(s).*A[:,Ip]+c*A[:,Iq]]

                    # after update matrix A calculate off_diag_norm of A
                    off_norm = off_diagonal_norm_cardoso(A)
                    
                end #if 

            end #for q
        end # for p
        push!(off_norm_array,off_norm)

    end # while
    return iter, A, V, off_norm_array
end
"""
    jdiag_cardoso(M,jthresh)

Only works for matrix with real valued entries. Based on [Matlab Code by Cardoso](https://www2.iap.fr/users/cardoso/jointdiag.html).

Input:
* A is a ``m × m × n`` matrix,(``A_1,...,A_n``), each with dimension ``m × m``
* thresh is a threshold for approximation stop, normally = 10e-8.

Output:
* V : is a  ``m × m`` matrix, which accumulates givens rotations G in each iteration.
* A : is a `` m × m × n`` matrix, which contains [``VA_1V'``,...,``VA_nV'``]
* iter: accumulates the iteration numbers
"""
function jdiag_cardoso(
    M::Vector{<:AbstractMatrix{<:Real}},
    jthresh :: Real;
    max_iter = 800,
    plot_convergence::Bool = false)
    # This version only works for matrix with real valued entries
    # Input:
    # A is a mxnm matrix,(A1,...,An),each with dimension mxm
    # thresh is a threshold for approximation stop, normally = 10e-8
    # Output:
    # iter: iteration number
    # A : is a mxnm matrix, which contains [VA1V',...,VAnV']
    # V : is a  mxm matrix, which accumulates givens rotations G in each iteration
    #off_norm_array: an array of summation of off-diagonal-norm of matrices set M in each iteration
    
    A = copy(hcat(M...))
    A = float.(A)
    m,nm = size(A)
    iter = 0
    V= Matrix{Float64}(I,m,m)
    off_norm = off_diagonal_norm_cardoso(A)
    off_norm_array = []

    
    # Initializing flag
    flag = true
    B = [1 0 0;0 1 1;0 -im im]

    if plot_convergence
        push!(off_norm_array, off_norm)
    end
    
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
                # calculate the parameters for givens rotation
                c = sqrt(0.5+angles[1]/2)
                s = 0.5*(angles[2])/c
                
                #update matrices A and V by a givens rotation
                if off_norm >jthresh
                    flag = true
                    pair = [p,q]
                    G = [c -conj(s); s c]

                    # update V, which accumulates givens rotations
                    V[:,pair] = V[:,pair]*G
                    
                    # update related two rows, p,q, of Real matrix A by one givens rotation 
                    A[pair,:] = G'*A[pair,:]

                    # update related two columns ,p,q of Real matrix A by one givens rotation
                    A[:,[Ip Iq]] = [c*A[:,Ip]+s*A[:,Iq];;-conj(s)*A[:,Ip]+c*A[:,Iq]]
                    
                    # after update matrix A calculate off_diag_norm of A
                    off_norm = off_diagonal_norm_cardoso(A)
                    
                end # if
            end # for q
        end # for p

        if plot_convergence
            push!(off_norm_array,off_norm)
        end
        
    end # while
    return V, A, off_norm_array
end




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
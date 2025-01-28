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
    jthresh::Real;
    plot_convergence::Bool = false)

    A = copy(hcat(M...))
    m,nm = size(A)
    iter = 0

    V= Matrix{Float64}(I,m,m)

    # Initializing flag
    flag = true
    B = [1 0 0;0 1 1;0 -im im]
    
    error_array = [] 
    if plot_convergence
        push!(error_array, frobenius_offdiag_norm(reverse_hcat(A)))
    end

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

        if plot_convergence
            push!(error_array, frobenius_offdiag_norm(reverse_hcat(A)))
        end
    end
    
    return V, A, error_array, iter
end

function reverse_hcat(A)
    n_dim = size(A, 1)
    n_matrix = size(A, 2)
    return cat([reshape(A, n_dim, n_dim, n_matrix ÷ n_dim)[:, :, index] for index in 1:n_matrix ÷ n_dim]..., dims=3)
end
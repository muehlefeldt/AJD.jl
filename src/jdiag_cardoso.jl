function jdiag_cardoso(M,jthresh)
    # Version02, 11.Dec.2024
    # This version only works for matrix with real valued entries
    # Input:
    # A is a mxnm matrix,(A1,...,An),each with dimension mxm
    # thresh is a threshold for approximation stop, normally = 10e-8
    # Output:
    # V : is a  mxm matrix, which accumulates givens rotations G in each iteration
    # A : is a mxnm matrix, which contains [VA1V',...,VAnV']
    # iter: accumulates the iteration numbers
    A = copy(M)
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

# generating test matrices which are positive definite and symmetric
function generate_psd_matrix(n::Int)
    """
    生成一个 n x n 的随机正定对称矩阵
    generate one nxn positive symmetric matrix
    """
    A = randn(n, n)  # 随机矩阵 a random matrix
    return A * A'    # 保证对称性和正定性 return A*A^T to garentee symmetric and positive definite
end

function generate_stacked_psd_matrices(n::Int, count::Int)
    """
    生成 count 个 n x n 的随机正定对称矩阵，并将它们拼接成一个 n x (n * count) 的大矩阵
    generate count-mal n x n positive definite matrices by casting them together
    """
    psd_matrices = [generate_psd_matrix(n) for _ in 1:count]  # 生成矩阵列表
    return hcat(psd_matrices...)  # 按列拼接
end
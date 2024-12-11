module AJD
function testAJD(A,jthresh)
    # Version02, 11.Dec.2024
    # This version only works for matrix with real valued entries
	# Input:
	# A is a mxnm matrix,(A1,...,An),each with dimension nxn
	# thresh is a threshold for approximation stop, normally = 10e-8
    # Output:
    # V : is a  mxm matrix, which accumulates givens rotations G in each iteration
    # A : is a mxnm matrix, which contains [VA1V',...,VAnV']
    
	m,nm = size(A)

	V= Matrix{Float64}(I,m,m)

	# Initializing flag
	flag = true
	B = [1 0 0;0 1 1;0 -im im]
	

	while flag
		flag = false


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
	return A,V
end
end

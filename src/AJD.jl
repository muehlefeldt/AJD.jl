module AJD
function testAJD(A,jthresh)
	# Input:
	# A is a mxnm matrix,(A1,...,An),each with dimension nxn
	# thresh is a threshold for approximation stop, normally = 10e-8
	m,nm = size(A)

	V= Matrix{ComplexF64}(I,m,m)

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
				angles = F.vectors[:,sp[3]]
				if angles[1] < 0
					angles =-angles
				end
				# 计算 givens rotation 的参数
				c = sqrt(0.5+angles[1]/2)
				s = 0.5*(angles[2]-im*angles[3])/c
				
				#update matrices A and V by a givens rotation
				if abs(s)>jthresh

					flag = true
					pair = [p,q]
					G = [c -conj(s); s c]

					
					V[:,pair] = V[:,pair]*G
					

					
					# Real matrix A
					A[pair,:] = G'*A[pair,:]

					# complex matrix
					# M[pair,:] = G'*M[pair,:]
					



					#####REAL
					A[:,Ip] = c*A[:,Ip]+s*A[:,Iq] 
					#
					# M[:,Ip] = c*M[:,Ip]+s*M[:,Iq] 
					#####REAL A
					A[:,Iq] = -conj(s)*A[:,Ip]+c*A[:,Iq]

					#complex matrix M
					# M[:,Iq] = -conj(s)*M[:,Ip]+c*M[:,Iq]
					
				end

			end
		end
	end
	return A,V
end
end

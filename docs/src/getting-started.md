```@meta
CurrentModule = AJD
```

# Getting Started Guide
This module aims to implement the Approximate Joint Diagonalization found in the paper: 

[1] *Jacobi Angles for Simultaneous Diagonalization* by Jean-François Cardoso and Antoine Souloumiac: SIAM journal on matrix analysis and applications 17.1 (1996): 161-164.

At some point the module will also contain the algorithm proposed in the paper: 

[2] *A fast algorithm for joint diagonalization with non-orthogonal transformations and its application to blind source separation* by Andreas Ziehe, Pavel Laskov, Guido Nolte, and Klaus-Robert Müller published in *The Journal of Machine Learning Research*, 5:777– 800, 2004.

For the implementation of the Approximate Joint Diagonalization in [1] three algorithms are used for comparison. 

The `AJD.jdiag_cardoso` function is found under the resource: [Matlab Code Cardoso](https://www2.iap.fr/users/cardoso/jointdiag.html) and works for an input of Array of concatenated real matrices e.g.

```julia 
rows = 3
M = 1.0*[I(rows) I(rows)]
function jdiag_cardoso(M,threshold)
return A,V
end
#with A as the diagonlized Matrices
```

The second implementation `AJD.jdiag_gabrieldernbach` is a python implementation parsed into julia for matrices consisting of real values found under the resource: [Python Code gabrieldernbach](https://github.com/gabrieldernbach/approximate_joint_diagonalization/).

Another method for the function exists which takes matrices with complex values, which uses the paper found in [1] the aformentioned algorithm and loosely the algorithm [Python Code edouardpineau](https://github.com/edouardpineau/Time-Series-ICA-with-SOBI-Jacobi) (also used for the third implementation):

```julia 
function jdiag_gabrieldernbach(A::Vector{Matrix{Float64}}; threshold = eps(), max_iter = 1000)
function jdiag_gabrieldernbach(A::Vector{Matrix{ComplexF64}}; threshold = eps(), max_iter = 1000)
return A,V
#e.g. for real: A = (1.0)* [Matrix(I, 6, 6) , Matrix(I, 6, 6)]
#for complex: A = [[ 1.0 0.0 1.0*im; 0.0 2.0 0.0; 1.0*im 0.0 1.0],[ 1.0 0.0 1.0*im; 0.0 2.0 0.0; 1.0*im 0.0 1.0]] (we know this is diagonalizable)
```

The third code is based on the code found here: [Python Code edouardpineau](https://github.com/edouardpineau/Time-Series-ICA-with-SOBI-Jacobi) and is a bit more versatile since it takes a Vector of AbstractMatrices as an input and works for real and complex matrices, as well as hermitian matrices included in the module `PosDefManifold.jl`, which are used in the `Diagonalizations.jl` package too. (Hermitian Matrices are included due to possible integration of other functions defined in the `Diagonalizations.jl` package)

```julia 
function jdiag_edourdpineau(X::Vector{M}; iter=100, eps=1e-3) where {T<:Union{Real,Complex},M<:AbstractMatrix{T}}
return A,V
```

A minimal working example can be found below, which uses the function `diagonalize(A::Vector{Matrix}; algorithm::String)`, which is the **only** function exported from the module (for now).

```julia 
function diagonalize(A::Vector{Matrix}; algorithm::String)
    if algorithm == "jdiag"
        return jdiag_gabrieldernbach(A)
    end
    if algorithm =="jdiag_cardoso"
        return jdiag_cardoso(A, 10e-8)
    end
    if algorithm == "jdiag_edourdpineau"
        return jdiag_edourdpineau(A)
    end
    return error
end
```
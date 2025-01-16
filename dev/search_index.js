var documenterSearchIndex = {"docs":
[{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"CurrentModule = AJD","category":"page"},{"location":"getting-started/#Getting-Started-Guide","page":"Getting Started","title":"Getting Started Guide","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"This guide provides information on implemented algorithms and ","category":"page"},{"location":"getting-started/#Background","page":"Getting Started","title":"Background","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"This module aims to implement the Approximate Joint Diagonalization (AJD) proposed in the paper: ","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"[1] Jacobi Angles for Simultaneous Diagonalization by Jean-François Cardoso and Antoine Souloumiac: SIAM journal on matrix analysis and applications 17.1 (1996): 161-164.","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"At some point the module will also contain the algorithm proposed in the paper: ","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"[2] A fast algorithm for joint diagonalization with non-orthogonal transformations and its application to blind source separation by Andreas Ziehe, Pavel Laskov, Guido Nolte, and Klaus-Robert Müller published in The Journal of Machine Learning Research, 5:777– 800, 2004.","category":"page"},{"location":"getting-started/#JDiag-Implementations","page":"Getting Started","title":"JDiag Implementations","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"For the implementation of the AJD in [1] three algorithms are used for comparison.","category":"page"},{"location":"getting-started/#AJD.jdiag_cardoso","page":"Getting Started","title":"AJD.jdiag_cardoso","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The AJD.jdiag_cardoso function is based on Matlab Code by Cardoso. Basic usage with two identity matrices as input:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"using AJD\nusing LinearAlgebra\n\n# Generate a vector of matrices with real values as input.\ntest_input = (1.0) * [Matrix(I, 6, 6) , Matrix(I, 6, 6)]\n\n# Calculate diagonlization of the input matrices.\ndiagonalize(test_input, algorithm=\"jdiag_cardoso\")","category":"page"},{"location":"getting-started/#AJD.jdiag_gabrieldernbach","page":"Getting Started","title":"AJD.jdiag_gabrieldernbach","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The second implementation AJD.jdiag_gabrieldernbach is based on a Python implementation by Gabrieldernbach for matrices consisting of real and complex values.","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Another version of the function exists which takes matrices with complex values, which uses the paper found in [1] the aformentioned algorithm and loosely the algorithm Python Code edouardpineau (also used for the third implementation).","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"using AJD\nusing LinearAlgebra\n\n# Generate a vector of matrices with real values as input.\ntest_input = (1.0) * [Matrix(I, 6, 6) , Matrix(I, 6, 6)]\n\n# Calculate diagonlization of the input matrices.\ndiagonalize(test_input, algorithm=\"jdiag_gabrieldernbach\")","category":"page"},{"location":"getting-started/#AJD.jdiag_edourdpineau","page":"Getting Started","title":"AJD.jdiag_edourdpineau","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The third implementation is based on the Python code by Edouardpineau and is a bit more versatile since it takes a Vector of AbstractMatrices as an input and works for real and complex matrices, as well as hermitian matrices included in the module PosDefManifold.jl, which are used in the Diagonalizations.jl package too. (Hermitian Matrices are included due to possible integration of other functions defined in the Diagonalizations.jl package)","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"using AJD\nusing LinearAlgebra\n\n# Generate a vector of matrices with real values as input.\ntest_input = (1.0) * [Matrix(I, 6, 6) , Matrix(I, 6, 6)]\n\n# Calculate diagonlization of the input matrices.\ndiagonalize(test_input, algorithm=\"jdiag_edourdpineau\")","category":"page"},{"location":"getting-started/#Minimal-Working-Example","page":"Getting Started","title":"Minimal Working Example","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"A minimal working example for all algorithms is demonstrated below. The function diagonalize() is the only function exported from the module and supports access to all algorithms described above. If supported, a complex matrix is used as input as well.","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"using AJD\nusing LinearAlgebra\n\n# Generate a vector of matrices with real values as input.\ntestinput = (1.0) * [Matrix(I, 6, 6) , Matrix(I, 6, 6)]\n\n@info \"Jdiag\", diagonalize(testinput; algorithm = \"jdiag\")\n@info \"jdiag_edourdpineau\", diagonalize(testinput; algorithm = \"jdiag_edourdpineau\")\n@info \"jdiag_cardoso\", diagonalize(testinput; algorithm = \"jdiag_cardoso\")\n\ntestinput_imag = [[ 1.0 0.0 1.0*im; 0.0 2.0 0.0; 1.0*im 0.0 1.0],[ 1.0 0.0 1.0*im; 0.0 2.0 0.0; 1.0*im 0.0 1.0]]\n\n@info \"jdiag\",diagonalize(testinput_imag; algorithm = \"jdiag\")\n@info \"jdiag_edourdpineau\",diagonalize(testinput_imag; algorithm = \"jdiag_edourdpineau\")\n","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"For generating further testdata of real matrices the following function can be used (though the function might only work for diagonalize(input; algorithm = \"jdiag_edourdpineau\")):","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"using LinearAlgebra\nfunction random_normal_commuting_matrices(n::Int, m::Int)\n    Q, _ = qr(rand(n,n))\n    Q = Matrix(Q)\n    return [Q*Diagonal(rand(n))*Q' for _ in 1:m]\nend\n","category":"page"},{"location":"getting-started/#Limitations-and-Known-issue","page":"Getting Started","title":"Limitations and Known issue","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Due to the different implementations the algorithm = \"jdiag\" and algorithm = \"jdiag_edourdpineau\" give different results in order of approx. 10^-1. However due to machine precision it is unclear how reliable those values really are.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = AJD","category":"page"},{"location":"#AJD.jl-Documentation","page":"Home","title":"AJD.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for the AJD.jl Package to perform Approximate Joint Diagonalization (AJD) on multiple matrices.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [AJD]","category":"page"},{"location":"#AJD.ALL_ALGORITHMS","page":"Home","title":"AJD.ALL_ALGORITHMS","text":"List of all algorithms keywords. Used in tests.\n\n\n\n\n\n","category":"constant"},{"location":"#AJD.COMLPLEX_ALGORITHMS","page":"Home","title":"AJD.COMLPLEX_ALGORITHMS","text":"List of all algorithms supporting complex matrices.\n\n\n\n\n\n","category":"constant"},{"location":"#AJD.ajd_benchmark-Tuple{Int64, Int64}","page":"Home","title":"AJD.ajd_benchmark","text":"ajd_benchmark(n_dims::Int, n_matrices::Int)\n\nRun benchmark of implemented algorithms with random inputs. Prints basic overview of median execution times. Returns BenchmarkGroup containing detailed results.\n\n\n\n\n\n","category":"method"},{"location":"#AJD.check_input-Tuple{Vector{<:AbstractMatrix{<:Number}}}","page":"Home","title":"AJD.check_input","text":"Check for valid input of diagonalize().\n\n\n\n\n\n","category":"method"},{"location":"#AJD.create_linear_filter-Tuple{Matrix{T} where T<:Number}","page":"Home","title":"AJD.create_linear_filter","text":"Create LinearFilter object as introduced by Diagonalizations.jl. Output of AJD.jl follows convention of Diagonalizations.jl and produces a LinearFilter.\n\n\n\n\n\n","category":"method"},{"location":"#AJD.diagonalize-Tuple{Vector{<:AbstractMatrix{<:Number}}}","page":"Home","title":"AJD.diagonalize","text":"diagonalize(\n    A::Vector{<:AbstractMatrix{<:Number}};\n    algorithm::String = \"jdiag_gabrieldernbach\",\n    max_iter::Int = 1000,\n    threshold::AbstractFloat = eps()\n    )\n\nCalculate joint diagonalization of multiple input matrices using the requested algorithms.\n\nMain function of the AJD package. Implemented algorithms at this point in time are limited to the JDiag algorithm  in different versions. Input of matrices to be diagonalized need to be a vector of matrices. The matrices can be of types Float64 or Complex. Limitations of the different implementations apply.\n\nSupported algorithms are jdiag_gabrieldernbach, jdiag_cardoso and jdiag_edourdpineau. See the Getting Started Guide for information on the algorithms.\n\n\n\n\n\n","category":"method"},{"location":"#AJD.frobenius_offdiag_norm-Union{Tuple{AbstractArray{T, 3}}, Tuple{T}} where T<:Number","page":"Home","title":"AJD.frobenius_offdiag_norm","text":"frobenius_off_diag_norm(A::Array{<:Number,3})\n\nInput\n\nA: Vector of matrices\n\nTakes an array namely the Array of matrices Ak and gets the offdiagonal elements and applies the frobenius norm (∑ |a{i,j}|^{2}).  Used for the jdiag_gabrieldernbach and FFD algorithm.\n\n\n\n\n\n","category":"method"},{"location":"#AJD.get_diag_elements-Tuple{Array}","page":"Home","title":"AJD.get_diag_elements","text":"get_diag_elements(A::Array)\n\nInput\n\nA: Vector of matrices\n\nTakes an array of matrices, takes the diagonal elements and returns the diagonal elements as a diagonal matrix D.\n\nExamples\n\n```jldoctest julia> A = ones(3,3,3);\n\njulia> AJD.getdiagelements(A) 3×3×3 Array{Float64, 3}: [:, :, 1] =  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0\n\n[:, :, 2] =  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0\n\n[:, :, 3] =  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0```\n\n\n\n\n\n","category":"method"},{"location":"#AJD.get_offdiag_elements-Tuple{Array{<:Number, 3}}","page":"Home","title":"AJD.get_offdiag_elements","text":"get_offdiag_elements(A::Array{<:Number,3})\n\nInput\n\nA: Vector of matrices\n\nTakes an array of matrices, copies it and sets the diagonal to zero to get the matrix E with only the offdiagonal elements of A.\n\nExamples\n\n```jldoctest julia> A = ones(3,3,3);\n\njulia> AJD.getoffdiagelements(A) 3×3×3 Array{Float64, 3}: [:, :, 1] =  0.0  1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0\n\n[:, :, 2] =  0.0  1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0\n\n[:, :, 3] =  0.0  1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0```\n\n\n\n\n\n","category":"method"},{"location":"#AJD.get_y_fdiag-Tuple{AbstractArray{<:Number}, AbstractArray{<:Number}, Int64, Int64}","page":"Home","title":"AJD.get_y_fdiag","text":"get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int, j::Int)\n\nInput\n\nD: Diagonal Matrix with offdiagonal elements set to zero\nE: Diagonal Matrix with diagonal elements set to zero\ni,j: Denotes the indexes of the matrices D and E\n\nCalculates the factor y_ij which is defined by: math ∑_{k} D_{j,j}^{k}E_{j,i}^{k}\n\n\n\n\n\n","category":"method"},{"location":"#AJD.get_z_fdiag-Tuple{AbstractArray{<:Number}, Int64, Int64}","page":"Home","title":"AJD.get_z_fdiag","text":"get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)\n\nInput\n\nD: Diagonal Matrix with offdiagonal elements set to zero\ni,j: Denotes the indexes the matrix D\n\nCalculates the factor z_ij which is defined by: ∑_{k} D_{i,i}^{k}D_{j,j}^{k}\n\n\n\n\n\n","category":"method"},{"location":"#AJD.is_commuting-Tuple{AbstractMatrix, AbstractMatrix}","page":"Home","title":"AJD.is_commuting","text":"is_commuting(A::AbstractMatrix, B::AbstractMatrix)\n\nInput:\n\nA: AbstractMatrix of dimension nxn\nB: AbstractMatrix of dimension nxn\n\nCheck if two matrices A, B are commuting.  A * B = B * A must hold. \n\n\n\n\n\n","category":"method"},{"location":"#AJD.is_same_size-Tuple{AbstractMatrix, AbstractMatrix}","page":"Home","title":"AJD.is_same_size","text":"is_same_size(A::AbstractMatrix, B::AbstractMatrix)\n\nInput:\n\nA: AbstractMatrix of variable size\nB: AbstractMatrix of variable size\n\nWill return true if dimension of A and B is matching otherwise false.\n\n\n\n\n\n","category":"method"},{"location":"#AJD.isstrictly_diagonally_dominant-Tuple{AbstractMatrix}","page":"Home","title":"AJD.isstrictly_diagonally_dominant","text":"isstrictly_diagonally_dominant(A::AbstractMatrix)\n\nInput:\n\nA: AbstractMatrix\n\nUsed for the FFDiag Algorithm to define whether the Matrix A is strictly diagonally dominant and therefore has an Inverse or not. A matrix is strictly dominant if:|a_{ii}| > \\sum |a_{ij}|, i ≠ j\n\n\n\n\n\n","category":"method"},{"location":"#AJD.jdiag_cardoso-Tuple{Vector{<:AbstractMatrix{<:Real}}, Real}","page":"Home","title":"AJD.jdiag_cardoso","text":"jdiag_cardoso(M,jthresh)\n\nOnly works for matrix with real valued entries. Based on Matlab Code by Cardoso.\n\nInput:\n\nA is a mxnm matrix,(A1,...,An), each with dimension mxm\nthresh is a threshold for approximation stop, normally = 10e-8.\n\nOutput:\n\nV : is a  mxm matrix, which accumulates givens rotations G in each iteration.\nA : is a mxnm matrix, which contains [VA1V',...,VAnV']\niter: accumulates the iteration numbers\n\n\n\n\n\n","category":"method"},{"location":"#AJD.jdiag_edourdpineau-Union{Tuple{Vector{M}}, Tuple{M}, Tuple{T}} where {T<:Number, M<:AbstractMatrix{T}}","page":"Home","title":"AJD.jdiag_edourdpineau","text":"jdiag_edourdpineau(X::Vector{M}; iter=100, eps=1e-3)\n    where {T<:Union{Real,Complex},M<:AbstractMatrix{T}}\n\nDiagonalize a set of matrices using the Jacobi method (\"Jacobi Angles for Simultaneous Diagonalization\"). Code adapted from Edouardpineaus Python implementation\n\n\n\n\n\n","category":"method"},{"location":"#AJD.jdiag_gabrieldernbach!-Union{Tuple{Vector{M}}, Tuple{M}, Tuple{T}} where {T<:Real, M<:AbstractMatrix{T}}","page":"Home","title":"AJD.jdiag_gabrieldernbach!","text":"(1) jdiag_gabrieldernbach(A::Vector{Matrix{Float64}}; threshold = eps(), max_iter = 1000)\n\nJDiag algorithm based on the implementation by Gabrieldernbach in Python.\n\nSource: https://github.com/gabrieldernbach/approximatejointdiagonalization/blob/master/jade/jade_cpu.py\n\n(2) jdiag_gabrieldernbach(A::Vector{Matrix{ComplexF64}}; threshold = eps(), max_iter = 1000)\n\nJDiag algorithm for complex matrices based on the implementation by Gabrieldernbach in Python, the Cardoso Paper and the code      of https://github.com/edouardpineau/Time-Series-ICA-with-SOBI-Jacobi.\n\n\n\n\n\n","category":"method"},{"location":"#AJD.plot_convergence_lineplot-Tuple{AbstractArray, String}","page":"Home","title":"AJD.plot_convergence_lineplot","text":"plot_convergence_lineplot(error_array::AbstractArray, name::String)\n\nPlot the convergence error as recorded during the algorithm execution.\n\n\n\n\n\n","category":"method"},{"location":"#AJD.plot_matrix_heatmap-Tuple{AbstractArray, AbstractArray}","page":"Home","title":"AJD.plot_matrix_heatmap","text":"plot_matrix_heatmap(filter::AbstractMatrix, diag_matrices::AbstractMatrix)\n\nPlot a heatmap of the calculated filter matrix and the mean of the diagonlized matrices. Complex matrices are reduced to real matrices for the plot.\n\n\n\n\n\n","category":"method"},{"location":"#AJD.random_commuting_matrices-Tuple{Int64, Int64}","page":"Home","title":"AJD.random_commuting_matrices","text":"Generate m random commuting matrices of size n x n These will produce all real rotation matrices using the Jacobi method\n\nMi * Mj = Mj * Mi for all i,j\n\n\n\n\n\n","category":"method"},{"location":"#AJD.random_matrices-Tuple{Int64, Int64}","page":"Home","title":"AJD.random_matrices","text":"random_matrices(n::Int, m::Int)\n\nGenerate m random matrices of size n x n\n\n\n\n\n\n","category":"method"},{"location":"#AJD.random_normal_commuting_matrices-Tuple{Int64, Int64}","page":"Home","title":"AJD.random_normal_commuting_matrices","text":"random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)\n\nGenerate m random normal commuting matrices of size n x n These can be exactly diagonalized\n\nMi * Mj = Mj * Mi for all i,j Mi*Mi' = Mi'*Mi for all i\n\n\n\n\n\n","category":"method"},{"location":"#AJD.random_symmetric_matrices-Tuple{Int64, Int64}","page":"Home","title":"AJD.random_symmetric_matrices","text":"Generate m random symmetric matrices of size n x n\n\n\n\n\n\n","category":"method"}]
}

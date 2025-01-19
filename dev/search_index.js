var documenterSearchIndex = {"docs":
[{"location":"detailed-docs/#Detailed-Docs","page":"Detailed Docs","title":"Detailed Docs","text":"","category":"section"},{"location":"detailed-docs/","page":"Detailed Docs","title":"Detailed Docs","text":"","category":"page"},{"location":"detailed-docs/","page":"Detailed Docs","title":"Detailed Docs","text":"Modules = [AJD]","category":"page"},{"location":"detailed-docs/#AJD.ALL_ALGORITHMS","page":"Detailed Docs","title":"AJD.ALL_ALGORITHMS","text":"List of all algorithms keywords. Used in tests.\n\n\n\n\n\n","category":"constant"},{"location":"detailed-docs/#AJD.COMLPLEX_ALGORITHMS","page":"Detailed Docs","title":"AJD.COMLPLEX_ALGORITHMS","text":"List of all algorithms supporting complex matrices.\n\n\n\n\n\n","category":"constant"},{"location":"detailed-docs/#AJD.ajd_benchmark-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.ajd_benchmark","text":"ajd_benchmark(n_dims::Int, n_matrices::Int)\n\nRun benchmark of implemented algorithms with random inputs. Prints basic overview of median execution times. Returns BenchmarkGroup containing detailed results.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.check_input-Tuple{Vector{<:AbstractMatrix{<:Number}}}","page":"Detailed Docs","title":"AJD.check_input","text":"Check for valid input of diagonalize().\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.create_linear_filter-Tuple{Matrix{T} where T<:Number}","page":"Detailed Docs","title":"AJD.create_linear_filter","text":"Create LinearFilter object as introduced by Diagonalizations.jl. Output of AJD.jl follows convention of Diagonalizations.jl and produces a LinearFilter.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.diagonalize-Tuple{Vector{<:AbstractMatrix{<:Number}}}","page":"Detailed Docs","title":"AJD.diagonalize","text":"diagonalize(\n    A::Vector{<:AbstractMatrix{<:Number}};\n    algorithm::String = \"jdiag_gabrieldernbach\",\n    max_iter::Int = 1000,\n    threshold::AbstractFloat = eps(),\n    plot_matrix::Bool = false,\n    plot_convergence::Bool = false\n    )\n\nCalculate joint diagonalization of multiple input matrices M_k.\n\nMain function of the AJD package. Implemented algorithms are JDiag and FFDiag. Input of matrices M_k need to be a vector of matrices. The matrices can be of types Float64 or Complex.\n\nSupported algorithms are jdiag_gabrieldernbach, jdiag_cardoso and jdiag_edourdpineau. See the Getting Started Guide for information on the algorithms. Test\n\nM_k.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.frobenius_offdiag_norm-Union{Tuple{AbstractArray{T, 3}}, Tuple{T}} where T<:Number","page":"Detailed Docs","title":"AJD.frobenius_offdiag_norm","text":"frobenius_offdiag_norm(A::AbstractArray{T,3})::Real where {T<:Number}\n\nInput\n\nA: Vector of matrices\n\nTakes an array namely the Array of matrices Ak and gets the offdiagonal elements and applies the frobenius norm (``\\sum |a{i,j}|^{2}``). \n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.generate_correlation_matrix-Tuple{Any, Any}","page":"Detailed Docs","title":"AJD.generate_correlation_matrix","text":"generate_correlation_matrix(signal_one_data,signal_two_data)\n\nsignalonedata: Array of dimension (NxM)\nsignaltwodata: Array of dimension (NxM)\n\nCalculates correlation matrix between observations x_i(t) and x_i(t+τ).\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.generate_testdata-Tuple{AbstractArray{<:Function}, AbstractMatrix{<:Number}}","page":"Detailed Docs","title":"AJD.generate_testdata","text":"generate_testdata(signal_sources::AbstractArray{<:Function}, mixing_matrix::AbstractMatrix{<:Number}; <keyword_arguments>)\n\nsignal_sources: Array of anonymous functions for generating time series data of the uncorrelated signals s of BSS e.g. [ s1 = x-> 1.4*sin(2x), s2 = 2.2sin(x)]\nmixing matrix: mixing matrix by which the signals s_j are multiplied to get the measurements/observations x_i\n\nFirst Calculates from a given array of functions resolved in the time domain, which simulate the uncorrelated signals s_j and a mixing matrix A, the measurements x_i with:\n\nx_i(t) = \\sum_{t = 1}^{T} a_{i,j} * s_j(t).\n\nThen time delayed correlation matrix between observations mentioned in [source] for a specified number in no_of_corr. \n\nArguments\n\ndelay::Number = 1: time delay between signals\nsample_time::Number = 10: length of single time series (same for all observations)\nnoofsamples::Int = 100: number of observations in the time series\nnoofcor::Int = 10: number of observations made over the entire measurement\n\nFor numbers to generate test data set of timeresolved function look here: https://doi.org/10.21595/jve.2021.21961 p.1709\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.generate_testdata-Tuple{AbstractArray}","page":"Detailed Docs","title":"AJD.generate_testdata","text":"generate_testdata(signal_sources::AbstractArray; delay::Number = 10, no_of_segments::Int = 10)\n\nsignalsources: Matrix of rowwise signals [``x1x2xn``]\ndelay: Time/index shift between observations to be correlated\nnoof:segments: Number of correlation matrices - 1 to be calculated. Puts signalsources into even segments to be correlated. If the number leads to uneven correlation will throw an error. \n\nGenerate Correlation Matrices for discrete observations x_i.\n\nKnown Issue\n\nIf your data has a lot of zeros inside the observations setting noofsegements too high will lead to NaN values since the variance of a vector of zeros is zero! You might want to manipulate your data or change the number of segments to be less!\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_diag_elements-Tuple{Array}","page":"Detailed Docs","title":"AJD.get_diag_elements","text":"get_diag_elements(A::Array)\n\nInput\n\nA: Vector of matrices\n\nTakes an array of matrices, takes the diagonal elements and returns the diagonal elements as a diagonal matrix D.\n\nExamples\n\n```jldoctest julia> A = ones(3,3,3);\n\njulia> AJD.getdiagelements(A) 3×3×3 Array{Float64, 3}: [:, :, 1] =  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0\n\n[:, :, 2] =  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0\n\n[:, :, 3] =  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0```\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_offdiag_elements-Tuple{Array{<:Number, 3}}","page":"Detailed Docs","title":"AJD.get_offdiag_elements","text":"get_offdiag_elements(A::Array{<:Number,3})\n\nInput\n\nA: Vector of matrices\n\nTakes an array of matrices, copies it and sets the diagonal to zero to get the matrix E with only the offdiagonal elements of A.\n\nExamples\n\n```jldoctest julia> A = ones(3,3,3);\n\njulia> AJD.getoffdiagelements(A) 3×3×3 Array{Float64, 3}: [:, :, 1] =  0.0  1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0\n\n[:, :, 2] =  0.0  1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0\n\n[:, :, 3] =  0.0  1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0```\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_y_fdiag-Tuple{AbstractArray{<:Number}, AbstractArray{<:Number}, Int64, Int64}","page":"Detailed Docs","title":"AJD.get_y_fdiag","text":"get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int, j::Int)\n\nInput\n\nD: Diagonal Matrix with offdiagonal elements set to zero\nE: Diagonal Matrix with diagonal elements set to zero\ni,j: Denotes the indexes of the matrices D and E\n\nCalculates the factor y_ij which is defined by: _k D_jj^kE_ji^k\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_z_fdiag-Tuple{AbstractArray{<:Number}, Int64, Int64}","page":"Detailed Docs","title":"AJD.get_z_fdiag","text":"get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)\n\nInput\n\nD: Diagonal Matrix with offdiagonal elements set to zero\ni,j: Denotes the indexes the matrix D\n\nCalculates the factor z_ij which is defined by: _k D_ii^kD_jj^k\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.is_commuting-Tuple{AbstractMatrix, AbstractMatrix}","page":"Detailed Docs","title":"AJD.is_commuting","text":"is_commuting(A::AbstractMatrix, B::AbstractMatrix)\n\nInput:\n\nA: AbstractMatrix of dimension nxn\nB: AbstractMatrix of dimension nxn\n\nCheck if two matrices A, B are commuting.  A * B = B * A must hold. \n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.is_same_size-Tuple{AbstractMatrix, AbstractMatrix}","page":"Detailed Docs","title":"AJD.is_same_size","text":"is_same_size(A::AbstractMatrix, B::AbstractMatrix)\n\nInput:\n\nA: AbstractMatrix of variable size\nB: AbstractMatrix of variable size\n\nWill return true if dimension of A and B is matching otherwise false.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.isstrictly_diagonally_dominant-Tuple{AbstractMatrix}","page":"Detailed Docs","title":"AJD.isstrictly_diagonally_dominant","text":"isstrictly_diagonally_dominant(A::AbstractMatrix)\n\nInput:\n\nA: AbstractMatrix\n\nUsed for the FFDiag Algorithm to define whether the Matrix A is strictly diagonally dominant and therefore has an Inverse or not. A matrix is strictly dominant if:a_ii  sum a_ij i  j\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.jdiag_cardoso-Tuple{Vector{<:AbstractMatrix{<:Real}}, Real}","page":"Detailed Docs","title":"AJD.jdiag_cardoso","text":"jdiag_cardoso(M,jthresh)\n\nOnly works for matrix with real valued entries. Based on Matlab Code by Cardoso.\n\nInput:\n\nA is a mxnm matrix,(A1,...,An), each with dimension mxm\nthresh is a threshold for approximation stop, normally = 10e-8.\n\nOutput:\n\nV : is a  mxm matrix, which accumulates givens rotations G in each iteration.\nA : is a mxnm matrix, which contains [VA1V',...,VAnV']\niter: accumulates the iteration numbers\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.jdiag_edourdpineau-Union{Tuple{Vector{M}}, Tuple{M}, Tuple{T}} where {T<:Number, M<:AbstractMatrix{T}}","page":"Detailed Docs","title":"AJD.jdiag_edourdpineau","text":"jdiag_edourdpineau(X::Vector{M}; iter=100, eps=1e-3)\n    where {T<:Union{Real,Complex},M<:AbstractMatrix{T}}\n\nDiagonalize a set of matrices using the Jacobi method (\"Jacobi Angles for Simultaneous Diagonalization\"). Code adapted from Edouardpineaus Python implementation\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.jdiag_gabrieldernbach!-Union{Tuple{Vector{M}}, Tuple{M}, Tuple{T}} where {T<:Real, M<:AbstractMatrix{T}}","page":"Detailed Docs","title":"AJD.jdiag_gabrieldernbach!","text":"(1) jdiag_gabrieldernbach(A::Vector{Matrix{Float64}}; threshold = eps(), max_iter = 1000)\n\nJDiag algorithm based on the implementation by Gabrieldernbach in Python.\n\nSource: https://github.com/gabrieldernbach/approximatejointdiagonalization/blob/master/jade/jade_cpu.py\n\n(2) jdiag_gabrieldernbach(A::Vector{Matrix{ComplexF64}}; threshold = eps(), max_iter = 1000)\n\nJDiag algorithm for complex matrices based on the implementation by Gabrieldernbach in Python, the Cardoso Paper and the code      of https://github.com/edouardpineau/Time-Series-ICA-with-SOBI-Jacobi.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.plot_convergence_lineplot-Tuple{AbstractArray, String}","page":"Detailed Docs","title":"AJD.plot_convergence_lineplot","text":"plot_convergence_lineplot(error_array::AbstractArray, name::String)\n\nPlot the convergence error as recorded during the algorithm execution.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.plot_matrix_heatmap-Tuple{AbstractArray, AbstractArray}","page":"Detailed Docs","title":"AJD.plot_matrix_heatmap","text":"plot_matrix_heatmap(filter::AbstractMatrix, diag_matrices::AbstractMatrix)\n\nPlot a heatmap of the calculated filter matrix and the mean of the diagonlized matrices. Complex matrices are reduced to real matrices for the plot.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.random_commuting_matrices-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.random_commuting_matrices","text":"Generate m random commuting matrices of size n x n These will produce all real rotation matrices using the Jacobi method\n\nMi * Mj = Mj * Mi for all i,j\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.random_matrices-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.random_matrices","text":"random_matrices(n::Int, m::Int)\n\nGenerate m random matrices of size n x n\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.random_normal_commuting_matrices-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.random_normal_commuting_matrices","text":"random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)\n\nGenerate m random normal commuting matrices of size n x n These can be exactly diagonalized\n\nMi * Mj = Mj * Mi for all i,j Mi*Mi' = Mi'*Mi for all i\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.random_symmetric_matrices-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.random_symmetric_matrices","text":"Generate m random symmetric matrices of size n x n\n\n\n\n\n\n","category":"method"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"CurrentModule = AJD","category":"page"},{"location":"getting-started/#Getting-Started-Guide","page":"Getting Started","title":"Getting Started Guide","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"This guide provides information on the basic usage of the AJD.jl package.","category":"page"},{"location":"getting-started/#Installation","page":"Getting Started","title":"Installation","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"To install the Package follow these instructions to add the package to a basic Julia environment or use the package in a Pluto notebook.","category":"page"},{"location":"getting-started/#Julia","page":"Getting Started","title":"Julia","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Create temporary environment in the Julia REPL:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"] activate --temp","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Add the package:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"] add https://github.com/muehlefeldt/AJD.jl","category":"page"},{"location":"getting-started/#Pluto","page":"Getting Started","title":"Pluto","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The use in Pluto notebooks is supported. Add the package in one cell of the notebook:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"begin\n    using Pkg\n    Pkg.activate(\"MyTestEnv\")\n    Pkg.add(url=\"https://github.com/muehlefeldt/AJD.jl\")\n    using AJD\nend","category":"page"},{"location":"getting-started/#Basic-usage","page":"Getting Started","title":"Basic usage","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The most basic usage of the package given a vector M of k square matrices (n times n):","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"using AJD\nfilter = diagonalize(M)","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The returned object filter of diagonalize() is a LinearFilter containing the matrix F as well as F. Assuming A in M these are used to diagonlize A by calculating D = F * A * F. D being the diagonalized matrix. In Julia you calculate:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"D = filter.iF * A * filter.F","category":"page"},{"location":"getting-started/#Advanced-options","page":"Getting Started","title":"Advanced options","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The main function diagonlaize() provides several options to be choosen by the user.","category":"page"},{"location":"getting-started/#Algorithms","page":"Getting Started","title":"Algorithms","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The package allows to choose diffrent algorithms to calculate the AJD.","category":"page"},{"location":"getting-started/#JDiag","page":"Getting Started","title":"JDiag","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Three implementations of the JDiag algorithm are avalaible [1]:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The Cardoso implementation is based on Matlab Code by Cardoso. Use the keyword algorithm=\"jdiag_cardoso\".\nBased on a Python implementation by Gabrieldernbach the second implementation is suitable for matrices consisting of real and complex values. Use the keyword algorithm=\"jdiag_gabrieldernbach\".\nThe third implementation of JDiag is based on the Python code by Edouardpineau also supports real and complex matrices, as well as hermitian matrices included in the module PosDefManifold.jl, which are used in the Diagonalizations.jl package too. Use the keyword algorithm=\"jdiag_edourdpineau\".","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"For example execute:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Generate 1000 matrices of 10 x 10 size.\nM = AJD.random_normal_commuting_matrices(10, 1000)\n# Diagonalize M using selected algorithm.\ndigonalize(M, algorithm=\"jdiag_edourdpineau\")","category":"page"},{"location":"getting-started/#FFDiag","page":"Getting Started","title":"FFDiag","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"One implementation of the FFDiag algorithm is available through the keyword algorithm=\"ffdiag\" [2].","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"For example execute:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Generate 1000 matrices of 10 x 10 size.\nM = AJD.random_normal_commuting_matrices(10, 1000)\n# Diagonalize M using selected algorithm.\ndigonalize(M, algorithm=\"ffdiag\")","category":"page"},{"location":"getting-started/#Plotting","page":"Getting Started","title":"Plotting","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Visual feedback is available by optional plots of the convergence behaviour and heatmaps of the filter matrix as well as the diagonlized matrices.","category":"page"},{"location":"getting-started/#Plot-Convergence","page":"Getting Started","title":"Plot Convergence","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"To plot the convergence behaviour use the keyword plot_convergence=true. An example output of random matrices can be generated by","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Generate 1000 matrices of 10 x 10 size.\nM = AJD.random_normal_commuting_matrices(10, 1000)\n# Diagonalize M.\ndiagonalize(M, plot_convergence=true)","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The generate plot:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"(Image: )","category":"page"},{"location":"getting-started/#Plot-Filter-and-Diagonalized-matrix","page":"Getting Started","title":"Plot Filter and Diagonalized matrix","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"To plot the filter matrix and the diagonalized matrices use the keyword plot_matrix=true. The diagonalized matrices D_k are summarised and only mean(D_k) is plotted to proivde an overview. Execute:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Generate 1000 matrices of 10 x 10 size.\nM = AJD.random_normal_commuting_matrices(10, 1000)\n# Diagonalize M.\ndiagonalize(M, plot_matrix=true)","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The generated plot:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"(Image: )","category":"page"},{"location":"getting-started/#Benchmarking","page":"Getting Started","title":"Benchmarking","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"All implementations of the algorithms can be benchmarked by executing:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Using matrices of size 10 x 10. Use 1000 random matrices for each benchmarked run.\najd_benchmark(10, 1000)","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Large input sizes require a lot compute time du to repeated execution.","category":"page"},{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"[1] J.-F. Cardoso and A. Souloumiac, \"Jacobi Angles for Simultaneous Diagonalization\", SIAM Journal on Matrix Analysis and Applications, vol. 17, no. 1, pp. 161-164. 1996. doi:10.1137/S0895479893259546.","category":"page"},{"location":"references/","page":"References","title":"References","text":"[2] A. Ziehe, P. Laskov, G. Nolte, K.-R. Müller, \"A Fast Algorithm for Joint Diagonalization with Non-orthogonal Transformations and its Application to Blind Source Separation\", Journal of Machine Learning Research, vol. 5, pp. 777-800. 07, 2004. Published: jmlr.org","category":"page"},{"location":"references/","page":"References","title":"References","text":"For further code examples and references please refer to the task document detailing the project.","category":"page"},{"location":"theoretical-background/#Theoretical-Background","page":"Theory","title":"Theoretical Background","text":"","category":"section"},{"location":"theoretical-background/","page":"Theory","title":"Theory","text":"This module implements Approximate Joint Diagonalization (AJD) of multiple matrices.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = AJD","category":"page"},{"location":"#AJD.jl-Documentation","page":"Home","title":"AJD.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for the AJD.jl Julia Package to perform Approximate Joint Diagonalization (AJD) on multiple matrices.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Getting Started Guide\nDetailed documentation\nTheory\nReferences","category":"page"},{"location":"#Main-AJD.jl-functions","page":"Home","title":"Main AJD.jl functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Below the main function of the AJD.jl package diagonlaize() is introduced below. For the entire documentation, including non exported functions, of the package refer to the detailed docs page.","category":"page"},{"location":"","page":"Home","title":"Home","text":"diagonalize(\n        A::Vector{<:AbstractMatrix{<:Number}};\n        algorithm::String = \"jdiag_gabrieldernbach\",\n        max_iter::Int = 1000,\n        threshold::AbstractFloat = eps(),\n        plot_matrix::Bool = false,\n        plot_convergence::Bool = false\n        )","category":"page"},{"location":"#AJD.diagonalize-Tuple{Vector{<:AbstractMatrix{<:Number}}}","page":"Home","title":"AJD.diagonalize","text":"diagonalize(\n    A::Vector{<:AbstractMatrix{<:Number}};\n    algorithm::String = \"jdiag_gabrieldernbach\",\n    max_iter::Int = 1000,\n    threshold::AbstractFloat = eps(),\n    plot_matrix::Bool = false,\n    plot_convergence::Bool = false\n    )\n\nCalculate joint diagonalization of multiple input matrices M_k.\n\nMain function of the AJD package. Implemented algorithms are JDiag and FFDiag. Input of matrices M_k need to be a vector of matrices. The matrices can be of types Float64 or Complex.\n\nSupported algorithms are jdiag_gabrieldernbach, jdiag_cardoso and jdiag_edourdpineau. See the Getting Started Guide for information on the algorithms. Test\n\nM_k.\n\n\n\n\n\n","category":"method"}]
}

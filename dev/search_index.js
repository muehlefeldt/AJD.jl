var documenterSearchIndex = {"docs":
[{"location":"detailed-docs/#Detailed-Docs","page":"Detailed Docs","title":"Detailed Docs","text":"","category":"section"},{"location":"detailed-docs/","page":"Detailed Docs","title":"Detailed Docs","text":"","category":"page"},{"location":"detailed-docs/","page":"Detailed Docs","title":"Detailed Docs","text":"Modules = [AJD]","category":"page"},{"location":"detailed-docs/#AJD.ALL_ALGORITHMS","page":"Detailed Docs","title":"AJD.ALL_ALGORITHMS","text":"List of all algorithms keywords. Used in tests.\n\n\n\n\n\n","category":"constant"},{"location":"detailed-docs/#AJD.COMLPLEX_ALGORITHMS","page":"Detailed Docs","title":"AJD.COMLPLEX_ALGORITHMS","text":"List of all algorithms supporting complex matrices.\n\n\n\n\n\n","category":"constant"},{"location":"detailed-docs/#AJD.ajd_benchmark-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.ajd_benchmark","text":"ajd_benchmark(n_dims::Int, n_matrices::Int)\n\nRun benchmark of implemented algorithms with random inputs. Prints basic overview of median execution times. Returns BenchmarkGroup containing detailed results.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.check_input-Tuple{Vector{<:AbstractMatrix{<:Number}}}","page":"Detailed Docs","title":"AJD.check_input","text":"Check for valid input of diagonalize().\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.create_linear_filter-Tuple{Matrix{T} where T<:Number}","page":"Detailed Docs","title":"AJD.create_linear_filter","text":"Create LinearFilter object as introduced by Diagonalizations.jl. Output of AJD.jl follows convention of Diagonalizations.jl and produces a LinearFilter.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.diagonalize-Tuple{Vector{<:AbstractMatrix{<:Number}}}","page":"Detailed Docs","title":"AJD.diagonalize","text":"diagonalize(\n    A::Vector{<:AbstractMatrix{<:Number}};\n    algorithm::String = \"jdiag_gabrieldernbach\",\n    max_iter::Int = 1000,\n    threshold::AbstractFloat = eps(),\n    plot_matrix::Bool = false,\n    plot_convergence::Bool = false\n    )\n\nCalculate joint diagonalization of multiple input matrices M_k.\n\nMain function of the AJD package. Implemented algorithms are JDiag and FFDiag. Input of matrices M_k need to be a vector of matrices. The matrices can be of types Float64 or Complex.\n\nSupported algorithms are jdiag_gabrieldernbach, jdiag_cardoso and jdiag_edourdpineau. See the Getting Started Guide for information on the algorithms. Test\n\nM_k.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.frobenius_offdiag_norm-Union{Tuple{AbstractArray{T, 3}}, Tuple{T}} where T<:Number","page":"Detailed Docs","title":"AJD.frobenius_offdiag_norm","text":"frobenius_offdiag_norm(A::AbstractArray{T,3})::Real where {T<:Number}\n\nInput\n\nA: Vector of matrices with size n x n x k \n\nTakes the offdiagonal elements of an Array of matrices A^k and applies the frobenius norm (sum a_ij^2). \n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.generate_correlation_matrix-Tuple{Any, Any}","page":"Detailed Docs","title":"AJD.generate_correlation_matrix","text":"generatecorrelationmatrix(signalonedata,signaltwodata)\n\nsignal_one_data: Array of dimension n x m\nsignal_two_data: Array of dimension n x m\n\nCalculates correlation matrix between observations x_i(t) and x_i(t+τ).\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.generate_testdata-Tuple{AbstractArray{<:Function}, AbstractMatrix{<:Number}}","page":"Detailed Docs","title":"AJD.generate_testdata","text":"generate_testdata(signal_sources::AbstractArray{<:Function}, mixing_matrix::AbstractMatrix{<:Number}; <keyword_arguments>)\n\nsignal_sources: Array of anonymous functions for generating time series data of the uncorrelated signals s_j of BSS e.g. [ s1 = x-> 1.4*sin(2x), s2 = 2.2sin(x)]\nmixing_matrix: mixing matrix by which the signals s_j are multiplied to get the measurements/observations x_i\n\nFirst Calculates from a given array of functions resolved in the time domain, which simulate the uncorrelated signals s_j and a mixing matrix A, the measurements x_i with:\n\nx_i(t) = sum_t = 1^T a_ij * s_j(t).\n\nThen time delayed correlation matrix between observations mentioned in [source] for a specified number in no_of_corr. \n\nArguments\n\ndelay::Number = 1: time delay between signals\nsample_time::Number = 10: length of single time series (same for all observations)\nno_of_samples::Int = 100: number of observations in the time series\nno_of_cor::Int = 10: number of observations made over the entire measurement\n\nFor numbers to generate test data set of timeresolved function look here: source p.1709\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.generate_testdata-Tuple{AbstractArray}","page":"Detailed Docs","title":"AJD.generate_testdata","text":"generate_testdata(signal_sources::AbstractArray; delay::Number = 10, no_of_segments::Int = 10)\n\nsignal_sources: Matrix of rowwise signals [x_1;x_2;...;x_n`]\ndelay`: Time/index shift between observations to be correlated\nno_of_segments: Puts signal_sources` into even segments to be correlated. If the number leads to uneven correlation will throw an error. \n\nGenerate Correlation Matrices for discrete observations x_i.\n\nKnown Issue\n\nIf your data has a lot of zeros inside the observations setting no_of_segements too high will lead to NaN values since the variance of a vector of zeros is zero! You might want to manipulate your data or change the number of segments to be less!\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_diag_elements-Tuple{Array}","page":"Detailed Docs","title":"AJD.get_diag_elements","text":"get_diag_elements(A::Array)\n\nA: Vector of matrices\n\nTakes an array of matrices and returns the diagonal elements as a diagonal matrix D.\n\nExamples\n\n```jldoctest julia> A = ones(3,3,3);\n\njulia> AJD.getdiagelements(A) 3×3×3 Array{Float64, 3}: [:, :, 1] =  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0\n\n[:, :, 2] =  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0\n\n[:, :, 3] =  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0```\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_diagonalization-Tuple{Vector{<:AbstractMatrix{<:Number}}}","page":"Detailed Docs","title":"AJD.get_diagonalization","text":"function get_diagonalization(\n    A::Vector{<:AbstractMatrix{<:Number}};\n    algorithm::String = \"jdiag_gabrieldernbach\",\n    max_iter::Int = 1000,\n    threshold::AbstractFloat = eps(),\n    only_plot::Symbol = :no_plot\n    )\n\nGet the actual diagonalization. Function is seperated from diagonalize() to facilitate plotting functionality in the REPL and Pluto. All implemented algorithms are called from this function. To generate the error histories of the algorithm runs, as used for the plots, select only_plot=:plot. Input is checked here as well.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_offdiag_elements-Tuple{Array{<:Number, 3}}","page":"Detailed Docs","title":"AJD.get_offdiag_elements","text":"get_offdiag_elements(A::Array{<:Number,3})\n\nInput\n\nA: Vector of matrices\n\nTakes an array of matrices and returns the offdiagonal elements of A.\n\nExamples\n\n```jldoctest julia> A = ones(3,3,3);\n\njulia> AJD.getoffdiagelements(A) 3×3×3 Array{Float64, 3}: [:, :, 1] =  0.0  1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0\n\n[:, :, 2] =  0.0  1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0\n\n[:, :, 3] =  0.0  1.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0```\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_plot-Tuple{AbstractArray, AbstractArray, AbstractArray, String}","page":"Detailed Docs","title":"AJD.get_plot","text":"get_plot(\n    filter::AbstractArray,\n    diag_matrices::AbstractArray, \n    error_array::AbstractArray, \n    name::String)\n\nIn case of plot is user selected this generates heatmap plot of the filter matrix and the mean of diagonlaised matrices. Also a lineplot of the error history of the algorithm calculation is created. A combined plot is returned.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_y_fdiag-Tuple{AbstractArray{<:Number}, AbstractArray{<:Number}, Int64, Int64}","page":"Detailed Docs","title":"AJD.get_y_fdiag","text":"get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int, j::Int)\n\nD: Diagonal Matrix with offdiagonal elements set to zero\nE: Diagonal Matrix with diagonal elements set to zero\ni,j: Denotes the indexes of the matrices D and E\n\nCalculates the factor y_ij which is defined by: _k D_jj^kE_ji^k\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.get_z_fdiag-Tuple{AbstractArray{<:Number}, Int64, Int64}","page":"Detailed Docs","title":"AJD.get_z_fdiag","text":"get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)\n\nD: Diagonal Matrix with offdiagonal elements set to zero\ni,j: Denotes the indexes the matrix D\n\nCalculates the factor z_ij which is defined by: _k D_ii^kD_jj^k\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.is_commuting-Tuple{AbstractMatrix, AbstractMatrix}","page":"Detailed Docs","title":"AJD.is_commuting","text":"is_commuting(A::AbstractMatrix, B::AbstractMatrix)\n\nA: AbstractMatrix of dimension n x n\nB: AbstractMatrix of dimension n x n\n\nCheck if two matrices A, B are commuting.  A * B = B * A must hold. \n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.is_same_size-Tuple{AbstractMatrix, AbstractMatrix}","page":"Detailed Docs","title":"AJD.is_same_size","text":"is_same_size(A::AbstractMatrix, B::AbstractMatrix)\n\nA: AbstractMatrix of variable size\nB: AbstractMatrix of variable size\n\nWill return true if dimension of A and B is matching otherwise false.\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.isstrictly_diagonally_dominant-Tuple{AbstractMatrix}","page":"Detailed Docs","title":"AJD.isstrictly_diagonally_dominant","text":"isstrictly_diagonally_dominant(A::AbstractMatrix)\n\nA: AbstractMatrix\n\nUsed for the FFDiag Algorithm to define whether the Matrix A is strictly diagonally dominant and therefore has an Inverse or not. A matrix is strictly dominant if: a_ii  sum a_ij i  j\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.jdiag_cardoso-Tuple{Vector{<:AbstractMatrix{<:Real}}, Real}","page":"Detailed Docs","title":"AJD.jdiag_cardoso","text":"jdiag_cardoso(M,jthresh)\n\nOnly works for matrix with real valued entries. Based on Matlab Code by Cardoso.\n\nInput:\n\nA is a mxnm matrix,(A1,...,An), each with dimension mxm\nthresh is a threshold for approximation stop, normally = 10e-8.\n\nOutput:\n\nV : is a  mxm matrix, which accumulates givens rotations G in each iteration.\nA : is a mxnm matrix, which contains [VA1V',...,VAnV']\niter: accumulates the iteration numbers\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.jdiag_edourdpineau-Union{Tuple{Vector{M}}, Tuple{M}, Tuple{T}} where {T<:Number, M<:AbstractMatrix{T}}","page":"Detailed Docs","title":"AJD.jdiag_edourdpineau","text":"jdiag_edourdpineau(X::Vector{M}; iter=100, eps=1e-3)\n    where {T<:Union{Real,Complex},M<:AbstractMatrix{T}}\n\nDiagonalize a set of matrices using the Jacobi method (\"Jacobi Angles for Simultaneous Diagonalization\"). Code adapted from Edouardpineaus Python implementation\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.jdiag_gabrieldernbach!-Union{Tuple{Vector{M}}, Tuple{M}, Tuple{T}} where {T<:Real, M<:AbstractMatrix{T}}","page":"Detailed Docs","title":"AJD.jdiag_gabrieldernbach!","text":"(1) jdiag_gabrieldernbach(A::Vector{Matrix{Float64}}; threshold = eps(), max_iter = 1000)\n\nJDiag algorithm based on the implementation by Gabrieldernbach in Python.\n\nSource: Algorithm\n\n(2) jdiag_gabrieldernbach(A::Vector{Matrix{ComplexF64}}; threshold = eps(), max_iter = 1000)\n\nJDiag algorithm for complex matrices based on the implementation by Gabrieldernbach in Python, the Cardoso Paper and the code  of Algorithm\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.random_commuting_matrices-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.random_commuting_matrices","text":"Generate m random commuting matrices of size n x n These will produce all real rotation matrices using the Jacobi method\n\nM_i * M_j = M_j * M_i for all ij\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.random_matrices-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.random_matrices","text":"random_matrices(n::Int, m::Int)\n\nGenerate m random matrices of size n x n\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.random_normal_commuting_matrices-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.random_normal_commuting_matrices","text":"random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)\n\nGenerate m random normal commuting matrices of size n x n These can be exactly diagonalized\n\nM_i * M_j = M_j * M_i for all ij M_i*M_i^T = M_i^T*M_i for all i\n\n\n\n\n\n","category":"method"},{"location":"detailed-docs/#AJD.random_symmetric_matrices-Tuple{Int64, Int64}","page":"Detailed Docs","title":"AJD.random_symmetric_matrices","text":"Generate m random symmetric matrices of size n x n\n\n\n\n\n\n","category":"method"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"CurrentModule = AJD","category":"page"},{"location":"getting-started/#Getting-Started-Guide","page":"Getting Started","title":"Getting Started Guide","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"This guide provides information on the basic usage of the AJD.jl package.","category":"page"},{"location":"getting-started/#Installation","page":"Getting Started","title":"Installation","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"To install the Package follow these instructions to add the package to a basic Julia environment or use the package in a Pluto notebook.","category":"page"},{"location":"getting-started/#Julia","page":"Getting Started","title":"Julia","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Create temporary environment in the Julia REPL:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"] activate --temp","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Add the package:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"] add https://github.com/muehlefeldt/AJD.jl","category":"page"},{"location":"getting-started/#Pluto","page":"Getting Started","title":"Pluto","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The use in Pluto notebooks is supported. Add the package in one cell of the notebook:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"begin\n    using Pkg\n    Pkg.activate(\"MyTestEnv\")\n    Pkg.add(url=\"https://github.com/muehlefeldt/AJD.jl\")\n    using AJD\nend","category":"page"},{"location":"getting-started/#Basic-usage","page":"Getting Started","title":"Basic usage","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The most basic usage of the package given a vector M of k square matrices (n times n):","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"using AJD\nfilter = diagonalize(M)","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The returned object filter of diagonalize() is a LinearFilter containing the matrix F as well as F. Assuming A in M these are used to diagonlize A by calculating D = F * A * F. D being the diagonalized matrix. In Julia you calculate:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"D = filter.iF * A * filter.F","category":"page"},{"location":"getting-started/#Advanced-options","page":"Getting Started","title":"Advanced options","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The main function diagonlaize() provides several options to be choosen by the user.","category":"page"},{"location":"getting-started/#Algorithms","page":"Getting Started","title":"Algorithms","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The package allows to choose diffrent algorithms to calculate the AJD.","category":"page"},{"location":"getting-started/#JDiag","page":"Getting Started","title":"JDiag","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Three implementations of the JDiag algorithm are avalaible [1]:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The Cardoso implementation is based on Matlab Code by Cardoso. Use the keyword algorithm=\"jdiag_cardoso\".\nBased on a Python implementation by Gabrieldernbach the second implementation is suitable for matrices consisting of real and complex values. Use the keyword algorithm=\"jdiag_gabrieldernbach\".\nThe third implementation of JDiag is based on the Python code by Edouardpineau also supports real and complex matrices, as well as hermitian matrices included in the module PosDefManifold.jl, which are used in the Diagonalizations.jl package too. Use the keyword algorithm=\"jdiag_edourdpineau\".","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"For example execute:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Generate 1000 matrices of 10 x 10 size.\nM = AJD.random_normal_commuting_matrices(10, 1000)\n# Diagonalize M using selected algorithm.\ndigonalize(M, algorithm=\"jdiag_edourdpineau\")","category":"page"},{"location":"getting-started/#FFDiag","page":"Getting Started","title":"FFDiag","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"One implementation of the FFDiag algorithm is available through the keyword algorithm=\"ffdiag\" [2].","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"For example execute:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Generate 1000 matrices of 10 x 10 size.\nM = AJD.random_normal_commuting_matrices(10, 1000)\n# Diagonalize M using selected algorithm.\ndigonalize(M, algorithm=\"ffdiag\")","category":"page"},{"location":"getting-started/#Plotting","page":"Getting Started","title":"Plotting","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Visual feedback is available by optional plots of the convergence behaviour and heatmaps of the filter matrix as well as the diagonlized matrices.","category":"page"},{"location":"getting-started/#Plot-Convergence","page":"Getting Started","title":"Plot Convergence","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"To plot the convergence behaviour use the keyword plot_convergence=true. An example output of random matrices can be generated by","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Generate 1000 matrices of 10 x 10 size.\nM = AJD.random_normal_commuting_matrices(10, 1000)\n# Diagonalize M.\ndiagonalize(M, plot_convergence=true)","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The generate plot:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"(Image: )","category":"page"},{"location":"getting-started/#Plot-Filter-and-Diagonalized-matrix","page":"Getting Started","title":"Plot Filter and Diagonalized matrix","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"To plot the filter matrix and the diagonalized matrices use the keyword plot_matrix=true. The diagonalized matrices D_k are summarised and only mean(D_k) is plotted to proivde an overview. Execute:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Generate 1000 matrices of 10 x 10 size.\nM = AJD.random_normal_commuting_matrices(10, 1000)\n# Diagonalize M.\ndiagonalize(M, plot_matrix=true)","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"The generated plot:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"(Image: )","category":"page"},{"location":"getting-started/#Benchmarking","page":"Getting Started","title":"Benchmarking","text":"","category":"section"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"All implementations of the algorithms can be benchmarked by executing:","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"# Using matrices of size 10 x 10. Use 1000 random matrices for each benchmarked run.\najd_benchmark(10, 1000)","category":"page"},{"location":"getting-started/","page":"Getting Started","title":"Getting Started","text":"Large input sizes require a lot compute time du to repeated execution.","category":"page"},{"location":"references/#References","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"[1] J.-F. Cardoso and A. Souloumiac, \"Jacobi Angles for Simultaneous Diagonalization\", SIAM Journal on Matrix Analysis and Applications, vol. 17, no. 1, pp. 161-164. 1996. doi:10.1137/S0895479893259546.","category":"page"},{"location":"references/","page":"References","title":"References","text":"[2] A. Ziehe, P. Laskov, G. Nolte, K.-R. Müller, \"A Fast Algorithm for Joint Diagonalization with Non-orthogonal Transformations and its Application to Blind Source Separation\", Journal of Machine Learning Research, vol. 5, pp. 777-800. 07, 2004. Published: jmlr.org","category":"page"},{"location":"references/","page":"References","title":"References","text":"For further code examples and references please refer to the task document detailing the project.","category":"page"},{"location":"theoretical-background/#Theoretical-Background","page":"Theory","title":"Theoretical Background","text":"","category":"section"},{"location":"theoretical-background/","page":"Theory","title":"Theory","text":"This module implements Approximate Joint Diagonalization (AJD) of multiple matrices.","category":"page"},{"location":"testdata_generation/#Generating-Test-Data","page":"Test Data Generation","title":"Generating Test Data","text":"","category":"section"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"The function AJD.generate_testdata is currently not exported but will give you the opportunity to create your own testsets to be diagonalized, consisting of an implementation of time delayed correlation matrices. The implementation of the time delayed correlation matrices used are explained here (p.16). There are currently two ways to generate your own testdata, either using discrete signals or continous signals.","category":"page"},{"location":"testdata_generation/#Continous-signals","page":"Test Data Generation","title":"Continous signals","text":"","category":"section"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"The function AJD.generate_testdata(signal_sources::AbstractArray{<:Function}, mixing_matrix::AbstractMatrix{<:Number}; delay::Number = 1, sample_time::Number = 10, no_of_samples::Int = 100, no_of_cor::Int = 10) for continous signals is a bit complicated to be honest. The input signal_sources is an array of anonymous function like:","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"# used signals in signal_sources are similar to signals found in https://doi.org/10.21595/jve.2021.21961 p.1709\nsignal_sources = [t->1.6sin(2pi*5t+5)+2sin(2pi*20t+27)+0.5sin(2pi*100t)+1 , t->1.2(2pi*11t)+sin(2pi*2t)+0.7sin(2pi*111t+10)]","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"and relates to the unmixed signal used in BSS problems. For more information on BSS problems see: Ziehe (p.5 ff.). Source Signals are time dependent functions but do not necessarily have to be harmonic functions (at least for testing) but can also include t^t. Be careful though, NaN values can occur if the signal converges to zero over time!","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"The input mixing_matrix can be a random matrix with the only constraint being that it has to have the same column size as signals being recorded e.g. if s_j with j = 1,2 A has to have two columns.","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"# mixing matrix is the same as in https://doi.org/10.21595/jve.2021.21961 p.1709\nmixing_matrix =  [0.32 -0.43; -1.31 0.34]","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"Multiplying the Mixing Matrix and the signal sources leads to the observations x_i(t) used for the calculated time correlation matrix. (See documentation of function for further explanation of mathematics behind.)","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"The arguments sample_time is used as the ending time of the signal in seconds and  no_of_samples is the frequency/sampletime of the signal. no_of_cor is used as an indicator how many time delayed correlation matrices are to be calculated. The argument delay is used to indicate the time shift between two observations x(t) and x(t+tau).","category":"page"},{"location":"testdata_generation/#Example-continous-signal","page":"Test Data Generation","title":"Example continous signal","text":"","category":"section"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"using AJD\nsignal_sources = [x->1.6sin(2pi*5x+5)+2sin(2pi*20x+27)+0.5sin(2pi*100x)+1,x->1.2(2pi*11x)+sin(2pi*2x)+0.7sin(2pi*111x+10)]\nmixing_matrix = [0.32 -0.43; -1.31 0.34]\ntestset_data = AJD.generate_testdata(signal_sources,mixing_matrix)","category":"page"},{"location":"testdata_generation/#Discrete-signals","page":"Test Data Generation","title":"Discrete signals","text":"","category":"section"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"For realising discrete test data the files channels2_room69_mix.wav, channels3_room69_mix.wav and  channels3_room69_mix_shortened.wav are included in the repository. The data was generated using the github code found here - data consists of different soundfiles, which including a differing number of source signals. The function generate_testdata(signal_sources::AbstractArray; delay::Number = 10, no_of_segments::Int = 10) behaves differently for discrete signals than the continous counterpart. signal_sources is a matrix with the observations x_i(t) ordered rowwise [x_1(t); x_2(t); x_3(t)]. delay is the shift index by which the time delay is emulated. delay= 1000 means x[currcolumn + 1000]. `noofsegmentsdivides thesignalsourcesinto even segments to be correlated. Future implementation will also include the keywordpointspersegmentwhich will behave similar tosample_time` of analog counterpart, to be more intuitive.","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"If the segments are not equally divided the method will throw an error. Make sure the length of the observations is dividable by the number of segments or trim your data to be able to!","category":"page"},{"location":"testdata_generation/#Example-discrete-signal","page":"Test Data Generation","title":"Example discrete signal","text":"","category":"section"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"using WAV\ndata,fs = wavread(\"channels3_room69_mix.wav\")\ndata = data' \nusing AJD\ntestset_data = AJD.generate_testdata(data, delay = 1000,no_of_segments=6)\ndiagonalize(testset_data)\n\nusing WAV\ndata,fs = wavread(\"channels2_room69_mix.wav\")\ndata = data' \nusing AJD\ntestset_data = AJD.generate_testdata(data[:,1:320300], delay = 1000,no_of_segments=100)\ndiagonalize(testset_data)","category":"page"},{"location":"testdata_generation/#Generating-testset-from-[https://github.com/fakufaku/bss*speech*dataset](https://github.com/fakufaku/bss_speech_dataset)","page":"Test Data Generation","title":"Generating testset from https://github.com/fakufaku/bssspeechdataset","text":"","category":"section"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"If you want to generate your own test data or additional datasets refer to the github repository linked above and clone it. You'll need python!","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"For the code of the repo to work the following adaptations have to be made:","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"the generate_samples.py has to be changed on line 120. Change the deprecated np.floatto float","category":"page"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"The generation will take a while and (probably) uses around 10GB of storage space. You can always cancel the generation of testdata during the process and get some testsets generated until then. After generating your testset you might want to edit the soundfiles since some of the testsets have ","category":"page"},{"location":"testdata_generation/#Known-Issues","page":"Test Data Generation","title":"Known Issues","text":"","category":"section"},{"location":"testdata_generation/","page":"Test Data Generation","title":"Test Data Generation","text":"The discrete case won't work sometimes if the soundfile has a vector only consisting of zeros since the correlation of this vector will become NaN since variance is zero! Make sure in that case to decrease number of segments.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = AJD","category":"page"},{"location":"#AJD.jl-Documentation","page":"Home","title":"AJD.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for the AJD.jl Julia Package to perform Approximate Joint Diagonalization (AJD) on multiple matrices.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Getting Started Guide\nDetailed documentation\nTheory\nReferences","category":"page"},{"location":"#Main-AJD.jl-functions","page":"Home","title":"Main AJD.jl functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Below the main function of the AJD.jl package diagonlaize() is introduced below. For the entire documentation, including non exported functions, of the package refer to the detailed docs page.","category":"page"},{"location":"","page":"Home","title":"Home","text":"diagonalize(\n        A::Vector{<:AbstractMatrix{<:Number}};\n        algorithm::String = \"jdiag_gabrieldernbach\",\n        max_iter::Int = 1000,\n        threshold::AbstractFloat = eps(),\n        plot_matrix::Bool = false,\n        plot_convergence::Bool = false\n        )","category":"page"},{"location":"#AJD.diagonalize-Tuple{Vector{<:AbstractMatrix{<:Number}}}","page":"Home","title":"AJD.diagonalize","text":"diagonalize(\n    A::Vector{<:AbstractMatrix{<:Number}};\n    algorithm::String = \"jdiag_gabrieldernbach\",\n    max_iter::Int = 1000,\n    threshold::AbstractFloat = eps(),\n    plot_matrix::Bool = false,\n    plot_convergence::Bool = false\n    )\n\nCalculate joint diagonalization of multiple input matrices M_k.\n\nMain function of the AJD package. Implemented algorithms are JDiag and FFDiag. Input of matrices M_k need to be a vector of matrices. The matrices can be of types Float64 or Complex.\n\nSupported algorithms are jdiag_gabrieldernbach, jdiag_cardoso and jdiag_edourdpineau. See the Getting Started Guide for information on the algorithms. Test\n\nM_k.\n\n\n\n\n\n","category":"method"}]
}

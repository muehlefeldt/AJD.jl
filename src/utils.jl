using LinearAlgebra
using Diagonalizations
using Statistics: cor

"""
    random_matrices(n::Int, m::Int)

Generate m random matrices of size n x n
"""
function random_matrices(n::Int, m::Int)
    return [rand(n,n) for i in 1:m]
end

"""
Generate m random symmetric matrices of size n x n
"""
function random_symmetric_matrices(n::Int, m::Int)
    return [Symmetric(rand(n,n)) for i in 1:m]
end

"""
Generate m random commuting matrices of size n x n
These will produce all real rotation matrices using the Jacobi method

``M_i * M_j = M_j * M_i for all i,j``
"""
function random_commuting_matrices(n::Int, m::Int)
    P = rand(n,n)
    return [P*Diagonal(rand(n))*inv(P) for _ in 1:m] #change inv to cholesky if possible?
end

"""
    random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)
    
Generate m random normal commuting matrices of size n x n
These can be exactly diagonalized

``M_i * M_j = M_j * M_i for all i,j``
``M_i*M_i^{T} = M_i^{T}*M_i for all i``
"""
function random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)
    # Just like the function below, this produces hermitian an symmetric matrices, they are
    # just not annotated as such
    # Q, _ = qr(rand(complex ? ComplexF64 : Float64, n,n))
    Q, _ = qr(rand(complex ? ComplexF64 : Float64, n,n))
    Q = Matrix(Q)
    # if complex
    #     return [Q*Diagonal(rand(ComplexF64, n))*Q' for _ in 1:m]
    # end
    return [Q*Diagonal(rand(n))*Q' for _ in 1:m]
end

function random_normal_commuting_symmetric_matrices(n::Int, m::Int; complex::Bool=false)
    Q, _ = qr(rand(complex ? ComplexF64 : Float64, n,n))
    Q = Matrix(Q)
    if complex
        return [Hermitian(Q*Diagonal(rand(n))*Q') for _ in 1:m]
    end
    return [Symmetric(Q*Diagonal(rand(n))*Q') for _ in 1:m]
end

function get_test_data_complex_real(n::Int, m::Int)
    Q, _ = qr(rand(n,n))
    Q = Matrix(Q)
    
    C = [Q * Diagonal(rand(ComplexF64, n)) * Q' for _ in 1:m]
    R = [Q * Diagonal(rand(n)) * Q' for _ in 1:m]
    return [C..., R...]
end

"""
Create LinearFilter object as introduced by Diagonalizations.jl.
Output of AJD.jl follows convention of Diagonalizations.jl and produces a LinearFilter.
"""
function create_linear_filter(A::Matrix{T} where {T<:Number})
    args=("Approximate Joint Diagonalization", false)
    return LinearFilter(A, Matrix(A'), nothing, nothing, nothing, nothing, args...)
end

"""
    frobenius_offdiag_norm(A::AbstractArray{T,3})::Real where {T<:Number}
Input
* A: Vector of matrices with size n x n x k 

Takes the offdiagonal elements of an Array of matrices ``A^k`` and applies the frobenius norm (``\\sum |a_{i,j}|^{2}``). 
"""
function frobenius_offdiag_norm(A::AbstractArray{T,3})::Real where {T<:Number}
    norm = zero(real(T))
    for row in axes(A, 1), column in axes(A, 2), k in axes(A, 3)
        row == column && continue
        norm += abs2(A[row, column, k])
    end
    return norm
end
"""
    get_offdiag_elements(A::Array{<:Number,3})
Input
* A: Vector of matrices

Takes an array of matrices and returns the offdiagonal elements of A.

# Examples
```jldoctest
julia> A = ones(3,3,3);

julia> AJD.get_offdiag_elements(A)
3×3×3 Array{Float64, 3}:
[:, :, 1] =
 0.0  1.0  1.0
 1.0  0.0  1.0
 1.0  1.0  0.0

[:, :, 2] =
 0.0  1.0  1.0
 1.0  0.0  1.0
 1.0  1.0  0.0

[:, :, 3] =
 0.0  1.0  1.0
 1.0  0.0  1.0
 1.0  1.0  0.0```
"""
function get_offdiag_elements(A::Array{<:Number,3})
    #copy matrix A to not overwrite it with zeros
    E = copy(A)
    for row in axes(E,1) # maybe eachindex better - NG
        E[row,row,:] .= 0 
    end
    return E
end
"""
    get_diag_elements(A::Array)
* A: Vector of matrices

Takes an array of matrices and returns the diagonal elements as a diagonal matrix D.

# Examples

```jldoctest
julia> A = ones(3,3,3);

julia> AJD.get_diag_elements(A)
3×3×3 Array{Float64, 3}:
[:, :, 1] =
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0

[:, :, 2] =
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0

[:, :, 3] =
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0```
"""
function get_diag_elements(A::Array)
    rows,columns,k = size(A)
    if typeof(A) <: AbstractArray{<:Complex}
        D = complex.(zeros(rows, columns, k))
    else
        D = zeros(rows, columns, k)
    end
    for row in 1:rows
        D[row,row,:] = A[row,row,:]
    end
    return D
end

function sort_offdiag_elements(A::AbstractArray{<:Number,3})
    rows,_,k = size(A)
    no_el = rows*rows-rows #number of offdiagonal elements
    if typeof(A) <: AbstractArray{<:Complex}
        sorted_array = complex.(zeros(no_el*k))
    else
        sorted_array = zeros(no_el*k)
    end

    for matrices in 1:k
        sorted_array[1+(matrices-1)*(no_el):matrices*(no_el)] = sort_offdiag_elements(A[:,:,matrices])
    end
    return sorted_array
end

function sort_offdiag_elements(A::AbstractArray{<:Number,2})
    rows,_ = size(A)
    if typeof(A) <: AbstractArray{<:Complex}
        sorted_array = complex.(zeros(rows*rows-rows))
    else
        sorted_array = zeros(rows*rows-rows)
    end
    i = 1
    for row in 1:rows-1,column in row+1:rows
        #TODO: Benchmark if append! or push! is faster
        sorted_array[i] = A[row, column]
        i += 1
        sorted_array[i] = A[column,row]
        i +=1
        
    end
    return sorted_array
end

"""
    is_commuting(A::AbstractMatrix, B::AbstractMatrix) 
* A: AbstractMatrix of dimension n x n
* B: AbstractMatrix of dimension n x n

Check if two matrices A, B are commuting.
 A * B = B * A must hold. 
 """
function is_commuting(A::AbstractMatrix, B::AbstractMatrix)
    return isapprox(A*B, B*A)
end
"""
    is_same_size(A::AbstractMatrix, B::AbstractMatrix)
* A: AbstractMatrix of variable size
* B: AbstractMatrix of variable size

Will return true if dimension of A and B is matching otherwise false.
"""
function is_same_size(A::AbstractMatrix, B::AbstractMatrix)
    return size(A) == size(B)
end
"""
    isstrictly_diagonally_dominant(A::AbstractMatrix)
* A: AbstractMatrix

Used for the FFDiag Algorithm to define whether the Matrix A is strictly diagonally dominant and therefore has an Inverse or not.
A matrix is strictly dominant if: ``|a_{ii}| > \\sum |a_{ij}|, i ≠ j``
"""
function isstrictly_diagonally_dominant(A::AbstractMatrix)
    for i in eachindex(A[1:end, 1])
        if abs(sum(A[i,:])) - abs(A[i,i]) > abs(A[i,i]) ? true : false
            return false
        end
    end
    return true
end
"""
    get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)
* D: Diagonal Matrix with offdiagonal elements set to zero
* i,j: Denotes the indexes the matrix D

Calculates the factor ``z_{ij}`` which is defined by: `` ∑_{k} D_{i,i}^{k}D_{j,j}^{k} ``
"""
function get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)
    return sum(D[i,i,:].*D[j,j,:])
end
"""
    get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int, j::Int)
* D: Diagonal Matrix with offdiagonal elements set to zero
* E: Diagonal Matrix with diagonal elements set to zero
* i,j: Denotes the indexes of the matrices D and E

Calculates the factor ``y_{ij}`` which is defined by:
`` ∑_{k} D_{j,j}^{k}E_{j,i}^{k} ``
"""
function get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int,j::Int)
    return sum(D[j,j,:].*E[i,j,:])
end

"Check for valid input of diagonalize()."
function check_input(A::Vector{<:AbstractMatrix{<:Number}})
    # Input may be empty.
    if length(A) <= 0
        return false
    end
    # All matrices must be commuting and of same size.
    for index in 1:length(A)-1
        if !is_same_size(A[index], A[index+1])
            return false
        end
    end
    return true
end

function addrandomnoise(A::Vector{M};σ = 0.5,same_noise = true) where {T<:Number, M<:AbstractMatrix{T}}
    k = length(A)
    rows,columns = size(A[1])
    if same_noise == true
        R = randn(rows,columns)
        for index_k = 1:k
            A[index_k] = A[index_k] + σ*R*R'
        end
    else
        for index_k = 1:k
            R = randn(rows,columns)
            A[k] = A[k] + σ*R*R'
        end
    end
    return A
end
function addrandomnoise!(A::Vector{M};σ = 0.5,same_noise = true) where {T<:Number, M<:AbstractMatrix{T}}
    k = length(A)
    rows,columns = size(A[1])
    if same_noise == true
        R = randn(rows,columns) #noise matrix
        for index_k = 1:k
            A[index_k] = A[index_k] + σ*R*R'
        end
    else
        for index_k = 1:k
            R = randn(rows,columns)
            A[k] = A[k] + σ*R*R'
        end
    end
    return A
end

"""
  generate_correlation_matrix(signal_one_data,signal_two_data)
                                                                                          
* `signal_one_data`: Array of dimension n x m
* `signal_two_data`: Array of dimension n x m

Calculates correlation matrix between observations ``x_i(t)`` and ``x_i(t+τ)``.
"""
function generate_correlation_matrix(signal_one_data,signal_two_data)

    if size(signal_one_data) != size(signal_two_data)
        throw(ArgumentError("Signals have different sizes!"))
    end

    C = cor(signal_one_data,signal_two_data,dims = 2)
    return C
end

"""
    generate_testdata(signal_sources::AbstractArray{<:Function}, mixing_matrix::AbstractMatrix{<:Number}; <keyword_arguments>)

* `signal_sources`: Array of anonymous functions for generating time series data of the uncorrelated signals ``s_j`` of `BSS` e.g. [ s1 = x-> 1.4*sin(2x), s2 = 2.2sin(x)]
* `mixing_matrix`: mixing matrix by which the signals ``s_j`` are multiplied to get the measurements/observations ``x_i``

First Calculates from a given array of functions resolved in the time domain, which simulate the uncorrelated signals ``s_j`` and a mixing matrix A, the measurements ``x_i`` with:

``x_i(t) = \\sum_{t = 1}^{T} a_{i,j} * s_j(t) ``.

Then time delayed correlation matrix between observations mentioned in [source] for a specified number in `no_of_corr`. 

# Arguments

* `delay`::Number = 1: time delay between signals
* `sample_time`::Number = 10: length of single time series (same for all observations)
* `no_of_samples`::Int = 100: number of observations in the time series
* `no_of_cor`::Int = 10: number of observations made over the entire measurement

For numbers to generate test data set of timeresolved function look here: [source p.1709](https://doi.org/10.21595/jve.2021.21961)
"""
function generate_testdata(signal_sources::AbstractArray{<:Function}, mixing_matrix::AbstractMatrix{<:Number}; 
    delay::Number = 1, sample_time::Number = 10, 
    no_of_samples::Int = 100, no_of_cor::Int = 10)

    rows,columns = size(mixing_matrix)
    if columns != length(signal_sources)
        throw(ArgumentError("Signal source array and mixing matrix have different dimensions."))
    end

    #initialize the matrix to be diagonalized
    C = zeros(rows,rows,no_of_cor)

    for k in 0:no_of_cor-1
       
        x = zeros(rows,no_of_samples)
        x_delay = zeros(rows,no_of_samples)
        for row in 1:rows # axes won't work on Array of Functions
            for source in 1:length(signal_sources)
                #needs to be done since broadcasting on vector of anonymous functions doesn't seem to work
                #tried invoke and .|> but couldn't get it to work which is why iteration is necessary

                #observations at starting point t = 1+(k*T)
                x[row,:] = x[row,:] + mixing_matrix[row,source]*signal_sources[source].(range(k*sample_time+1,(k+1)*sample_time,length = no_of_samples))
                
                #time delayed observations
                x_delay[row,:] = x_delay[row,:] + mixing_matrix[row,source]*signal_sources[source].(range((k*sample_time+1+delay),(k+1)*sample_time+delay, length = no_of_samples))
            end
        end
        C[:,:,k+1] = generate_correlation_matrix(x,x_delay)
    end

    return C
end
"""
    generate_testdata(signal_sources::AbstractArray; delay::Number = 10, no_of_segments::Int = 10)
* `signal_sources``: Matrix of rowwise signals [``x_1``; ``x_2``;...; ``x_n``]
* `delay``: Time/index shift between observations to be correlated
* `no_of_segments`: Puts `signal_sources`` into even segments to be correlated. If the number leads to uneven correlation will throw an error. 

Generate Correlation Matrices for discrete observations ``x_i``.

# Known Issue

If your data has a lot of zeros inside the observations setting `no_of_segements` too high will lead to NaN values since the variance of a vector of zeros is zero!
You might want to manipulate your data or change the number of segments to be less!

"""
function generate_testdata(signal_sources::AbstractArray; 
    delay::Number = 10, no_of_segments::Int = 10)
    
    x = signal_sources
    rows,columns = size(signal_sources) 

    # if signal has length 100 and delay would be 99 there wouldn't be any data after observation 100 to correlate
    if columns/delay < 2
        throw(ArgumentError("Delay too big. Length of signals divided by delay less than 2. Delay shift would lead to array entry for non existent data."))
    end
    
    segmentation = columns/no_of_segments

    #in case segmentation leads to uneven segments i.e. 200 points and 
    #3 segments -> last observation wouldn't be considered since index will be 66

    if isinteger(segmentation)
        segmentation = Int(segmentation)
    else
        throw(ArgumentError("Number of Segments leads to segments of different sizes!"))
    end
    
    # C needs to be declared before the for loop otherwise won't be part of local scope 
    # if statement therefore declaring C inside for loop won't work
    # and declaring C as Matrix[] or Array[] will push first entry to be
    # Array[...] or Matrix[...] with all other entries being of type Matrix or Array

    # initialize the matrix set with the type of signal_sources
    C = Matrix{typeof(x[1])}[]

    for k in 1:no_of_segments-1
        x_new = x[:,(k-1)*segmentation+1:k*segmentation]
        x_delay = x[:,(k-1)*segmentation+1+delay:k*segmentation+delay]
        push!(C,generate_correlation_matrix(x_new,x_delay))
    end
    if isnan.(sum(C)) != zeros(rows,rows)
        throw(ArgumentError("Number of segments leads to NaN inside of correlation matrix. See Documentation for further Info."))
    end

    return C
end

"""
    function get_diagonalization(
        A::Vector{<:AbstractMatrix{<:Number}};
        algorithm::String = "jdiag_gabrieldernbach",
        max_iter::Int = 1000,
        threshold::AbstractFloat = eps(),
        only_plot::Symbol = :no_plot
        )

Get the actual diagonalization. Function is seperated from `diagonalize()` to facilitate plotting functionality in the REPL and Pluto.
All implemented algorithms are called from this function.
To generate the error histories of the algorithm runs, as used for the plots, select `only_plot=:plot`.
Input is checked here as well.
"""
function get_diagonalization(
    A::Vector{<:AbstractMatrix{<:Number}};
    algorithm::String = "jdiag_gabrieldernbach",
    max_iter::Int = 1000,
    threshold::AbstractFloat = eps(),
    only_plot::Symbol = :no_plot
    )
    if !check_input(A)
        throw(ArgumentError("Invalid input."))
    end

    # Do we need to get the error history from the algorithms?
    # If so selected, a performance panelty is to be expected.
    plot_convergence = only_plot == :plot

    if algorithm in ["jdiag", "jdiag_gabrieldernbach"]
        F, B, error_array = jdiag_gabrieldernbach!(A, max_iter = max_iter, threshold = threshold, plot_convergence=plot_convergence)

    elseif algorithm == "jdiag_edourdpineau"
        F, B, error_array = jdiag_edourdpineau(A, iter = max_iter)

    elseif algorithm == "jdiag_cardoso"
        if typeof(A) <: AbstractArray{<:AbstractArray{<:Real}} 
            F, B, error_array = jdiag_cardoso(A, threshold, plot_convergence=plot_convergence)
        else
            throw(ArgumentError("Not supported for set of Matrices containing imaginary values!"))
        end

    elseif algorithm in ["FFD", "ffd", "ffdiag"]
        F, B, error_array = FFD!(
            A,
            threshold=threshold,
            max_iter=max_iter,
            plot_convergence=plot_convergence
        )

    else
        # If no vaild algorithm selected, throw an error.
        throw(ArgumentError(
            "No valid algorithm selected. Available options:" * join(AJD.ALL_ALGORITHMS, ", ")
        ))
    end

    return F, B, error_array
end

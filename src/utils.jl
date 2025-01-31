using LinearAlgebra: qr, Diagonal, Hermitian, Symmetric
using Diagonalizations: LinearFilter
using Statistics: cor
using Random: rand, randn

"""
    random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)

Generate m random normal commuting matrices of size ``n × n``
These can be exactly diagonalized

``M_i M_j = M_j  M_i`` for all i,j
``M_i M_i' = M_i' M_i`` for all i
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
Create [LinearFilter](https://marco-congedo.github.io/Diagonalizations.jl/dev/Diagonalizations/#LinearFilter) object as introduced by [Diagonalizations.jl](https://marco-congedo.github.io/Diagonalizations.jl/dev/).
Output of AJD.jl follows convention of Diagonalizations.jl and produces a LinearFilter.
"""
function create_linear_filter(A::Matrix{T} where {T<:Number})
    args=("Approximate Joint Diagonalization", false)
    return LinearFilter(A, Matrix(A'), nothing, nothing, nothing, nothing, args...)
end

"""
    frobenius_offdiag_norm(A::AbstractArray{T,3})::Real where {T<:Number}
* A: Vector of matrices with size ``n × n × k``

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
* A: Vector of matrices

Takes an array of matrices and returns the offdiagonal elements of A.
"""
function get_offdiag_elements(A::Array{<:Number,3})
    #copy matrix A to not overwrite it with zeros
    E = copy(A)
    #give rows as eachindex in cartesian indexing instead of linear
    #eachrow would be nice however it is slower if eachrow.axes is used
    #according to @btime with replacing eachindex with eachrow.axes
    iterator = eachindex(IndexCartesian(),A[:,begin,begin])
    for row in iterator
        E[row,row,:] .= 0
    end
    return E
end

"""
    get_diag_elements(A::Array)
* A: Vector of matrices

Takes an array of matrices and returns the diagonal elements as a diagonal matrix D.
"""
function get_diag_elements(A::Array{<:Number,3})
    D = zeros(axes(A))
    iterator = eachindex(IndexCartesian(),A[:,begin,begin])
    #eachrow.axes might be nicer however it is slower than using eachindex
    #another way is axes(A,1) which might be faster
    for row in iterator
        D[row,row,:] = A[row,row,:]
    end
    return D
end

"""
    is_commuting(A::AbstractMatrix, B::AbstractMatrix)
* A: AbstractMatrix of dimension ``n × n``
* B: AbstractMatrix of dimension ``n × n``

Check if two matrices A, B are commuting.
 ``AB = BA`` must hold.
 """
function is_commuting(A::AbstractMatrix, B::AbstractMatrix)
    return isapprox(A*B, B*A)
end


"""
    get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)
* D: Diagonal Matrix with offdiagonal elements set to zero
* i,j: Denotes the indexes of matrix D

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

"""
    check_input(
        A::Vector{<:AbstractMatrix{<:Number}},
        max_iter::Int,
        threshold::AbstractFloat,
        )

Check input of diagonalize(). Validate matrices, threshold and max iteration as selected.
"""
function check_input(
    A::Vector{<:AbstractMatrix{<:Number}},
    max_iter::Int,
    threshold::AbstractFloat,
)

    # Input vector of matrices may be empty.
    if length(A) <= 0
        throw(ArgumentError("Invalid input: Vector of matrices must have length > 0."))
    end
    # All matrices must be of same size.
    if !allequal(size.(A))
        throw(ArgumentError("Invalid input: Vector of matrices must be of same size."))
    end

    # Last matrix needs to be checked for all zeros as earlier loop does not touch the last matrix.
    if any(iszero.(A))
        throw(ArgumentError("Invalid input: All zero matrix not allowed."))
    end

    # Max iteration must be 1 or higher.
    if max_iter <= 0
        throw(ArgumentError("Invalid input: Max iteration must be 1 or larger."))
    end

    # Too small threshold is non-sensical to procede with.
    if threshold < eps()
        throw(
            ArgumentError(
                "Invalid input: Threshold too small. Minimum: " * string(eps()),
            ),
        )
    end

    # Warning in case of very high threshold.
    if threshold > 1.0e-1
        @warn "Threshold very high. Recommend threshold of 1e-5 or smaller. Consider machine precision of your system."
    end
end

"""
    addrandomnoise(
        A::Vector{M};
        σ::AbstractFloat = 0.5,
        same_noise::Bool = true,
    ) where {T<:Number,M<:AbstractMatrix{T}}

Add ranom noise to a vector of matrices `A`.
If same noise selected same random matrix R used for all (!) matrices.
Completley random noise added to each matrix in `A` with `same_noise=false`.
"""
function addrandomnoise(
    A::Vector{M};
    σ::AbstractFloat = 0.5,
    same_noise::Bool = true,
) where {T<:Number,M<:AbstractMatrix{T}}
    k = length(A)
    rows, columns = size(A[1])

    if same_noise
        R = randn(rows, columns)
        for index_k = 1:k
            A[index_k] = A[index_k] + σ * R * R'
        end

    else
        for index_k in eachindex(A)
            R = randn(rows,columns)
            A[index_k]= A[index_k] + σ*R*R'
        end
    end
    return A
end

"""
    generate_correlation_matrix(
        signal_one_data::AbstractArray,
        signal_two_data::AbstractArray,
    )

Calculates correlation matrix between observations ``x_i(t)`` and ``x_i(t+τ)``.

Inputs:
* `signal_one_data`: Array of dimension ``n \\times m``
* `signal_two_data`: Array of dimension ``n \\times m``
"""
function generate_correlation_matrix(
    signal_one_data::AbstractArray,
    signal_two_data::AbstractArray,
)

    if size(signal_one_data) != size(signal_two_data)
        throw(ArgumentError("Signals have different sizes!"))
    end

    C = cor(signal_one_data, signal_two_data, dims = 2)
    return C
end

"""
    generate_testdata(signal_sources::AbstractArray{<:Function},
    mixing_matrix::AbstractMatrix{<:Number}; <keyword_arguments>)

* `signal_sources`: Array of anonymous functions for generating time series data of the uncorrelated signals ``s_j`` of `BSS` e.g. [ s1 = x-> 1.4*sin(2x), s2 = 2.2sin(x)]
* `mixing_matrix`: mixing matrix by which the signals ``s_j`` are multiplied to get the measurements/observations ``x_i``

First Calculates from a given array of functions, which simulate the uncorrelated signals ``s_j(t)`` and a mixing matrix A, the observations ``x_i``:

``x_i(t) = \\sum_{t = 1}^{T} a_{i,j} s_j(t) ``.

Then a number of time delayed correlation matrices specified by `no_of_corr` is calculated.
# Arguments

* `delay`::Number = 1: time delay between signals
* `sample_time`::Number = 10: length of single time series (same for all observations)
* `no_of_samples`::Int = 100: number of observations made during `sample_time`
* `no_of_cor`::Int = 10: number of observations made over the entire measurement
# Example for `signal_sources` and `mixing_matrix`
```julia
signal_sources = [x->1.6sin(2pi*5x+5)+2sin(2pi*20x+27)+0.5sin(2pi*100x)+1,x->1.2(2pi*11x)+sin(2pi*2x)+0.7sin(2pi*111x+10)]
mixing_matrix = [0.32 -0.43; -1.31 0.34]
```
"""
function generate_testdata(signal_sources::AbstractArray{<:Function}, mixing_matrix::AbstractMatrix{<:Number};
    delay::Number = 1, sample_time::Number = 10,
    no_of_samples::Int = 100, no_of_cor::Int = 10)

    rows,columns = size(mixing_matrix)
    if columns != length(signal_sources)
        throw(ArgumentError("Signal source array and mixing matrix have different dimensions (columns of matrix don't match signals in signal_sources)."))
    end

    #initialize the matrix to be diagonalized
    C = Matrix{}[]

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
        push!(C,generate_correlation_matrix(x,x_delay))
    end

    #convert C to appropriate type
    if all(isa.(C, Matrix{Float64})) == true
        C = convert(Vector{Matrix{Float64}},C)
    else
        C = convert(Vector{Matrix{ComplexF64}},C)
    end

    return C
end

"""
    generate_testdata(signal_sources::AbstractArray;
    delay::Number = 10, no_of_segments::Int = 10)
* `signal_sources`: Matrix of rowwise signals [``x_1``; ``x_2``;...; ``x_n``]
* `delay`: Time/index shift between observations to be correlated
* `no_of_segments`: Puts `signal_sources` into even segments to be correlated. If the number leads to uneven correlation will throw a warning if show_warning is true
* `show_warning`: If true will show a warning in case segments are uneven. Will lead to one less correlation matrix

Generate Correlation Matrices for discrete observations ``x_i``.

# Known Issue

If your data has a segment with variance close to 0 (e.g. due to all of the values being the same) the correlation matrix will have NaN values inside. Setting the number of segments to a lower value might help.
"""
function generate_testdata(signal_sources::AbstractArray;
    delay::Number = 10, no_of_segments::Int = 10, show_warning::Bool = true)

    x = copy(signal_sources)
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
        # for testing should be false! could be also done with Logging module and
        # disable_logging(Logging.Warn) but that would supress all warnings for the test case
        # which might not be desirable!
        if show_warning == true
            @warn "Number of Segments leads to segments of different sizes! Will skip last segment.\n"
        end
        no_of_segments -= 1
        segmentation = floor(Int64,segmentation)
    end

    # C needs to be declared before the for loop otherwise won't be part of local scope
    # if statement therefore declaring C inside for loop won't work
    # and declaring C as Matrix[] or Array[] will push first entry to be
    # Array[...] or Matrix[...] with all other entries being of type Matrix or Array

    # initialize the matrix set with the type of signal_sources
    # this might lead to an error if the type of signals differs from eachother however
    # if only waveform files are used for testdata generation this shouldn't fail. could
    #
    C = Matrix{}[]

    for k in 1:no_of_segments-1
        #won't work for offset arrays however since this is only for testdata generation
        #of wav files which are read in shouldn't matter
        x_t = x[:,(k-1)*segmentation+1:k*segmentation]
        x_delay = x[:,(k-1)*segmentation+1+delay:k*segmentation+delay]
        push!(C,generate_correlation_matrix(x_t,x_delay))
    end
    if isnan.(sum(C)) != zeros(rows,rows)
        throw(ArgumentError("Number of segments leads to NaN inside of correlation matrix. See Documentation for further Info."))
    end
    # convert C to appropriate type
    if all(isa.(C, Matrix{Float64})) == true
        C = convert(Vector{Matrix{Float64}},C)
    else
        C = convert(Vector{Matrix{ComplexF64}},C)
    end
    return C
end

"""
    generate_random_signals(
        no_of_signals::Int,
        no_of_samples::Int;
        seed = Xoshiro(),
        signal_type::DataType = Float64,
    )

*`no_of_signals`: used to declare how many signals ``s_{j}`` are inside of the testdata
*`no_of_samples`: used to declare how many samples each signal ``s_j`` has
*`seed`: seed for the `rand` function, might need to import the module Random
*`signal_type`: default is Float64 but ComplexF64 is possible as well

Similar to `random_commuting_matrices` function. Gives the opportunity to generate testdata randomly and pluck it into the discrete version of generate_testdata.
"""
function generate_random_signals(
    no_of_signals::Int,
    no_of_samples::Int;
    seed = Xoshiro(),
    signal_type::DataType = Float64,
)
    signals = Matrix{signal_type}(undef, no_of_signals, no_of_samples)
    for signal = 1:no_of_signals
        signals[signal, :] = rand(seed, signal_type, no_of_samples)
    end
    return signals
end

"""
    get_diagonalization(
        A::Vector{<:AbstractMatrix{<:Number}};
        algorithm::AbstractDiagonalization = JDiagGabrielDernbach(),
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
    algorithm::AbstractDiagonalization = JDiagGabrielDernbach(),
    max_iter::Int = 1000,
    threshold::AbstractFloat = eps(),
    only_plot::Symbol = :no_plot
    )

    plot_convergence = only_plot == :plot

    return _get_diagonalization(A, algorithm, max_iter, threshold, plot_convergence)
end

# Multiple dispatch for each algorithm
function _get_diagonalization(A, ::JDiagGabrielDernbach, max_iter, threshold, plot_convergence)
    jdiag_gabrieldernbach!(A, max_iter=max_iter, threshold=threshold, plot_convergence=plot_convergence)
end

function _get_diagonalization(A, ::JDiagEdourdPineau, max_iter, threshold, plot_convergence)
    jdiag_edourdpineau(A, iter=max_iter)
end

function _get_diagonalization(A, ::JDiagCardoso, max_iter, threshold, plot_convergence)
    jdiag_cardoso(A, threshold, plot_convergence=plot_convergence, max_iter=max_iter)
end

function _get_diagonalization(A, ::FFDiag, max_iter, threshold, plot_convergence)
    ffd(A, threshold=threshold, max_iter=max_iter, plot_convergence=plot_convergence)
end

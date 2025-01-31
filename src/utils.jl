using LinearAlgebra: qr, Diagonal, Hermitian, Symmetric
using Diagonalizations: LinearFilter
# using Statistics: cor
using Random: rand, randn

"""
    random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)

Generate m random normal commuting matrices of size ``n × n``
These can be exactly diagonalized

``M_i M_j = M_j  M_i`` for all i,j
``M_i M_i' = M_i' M_i`` for all i
"""
function random_normal_commuting_matrices(n::Int, m::Int; complex::Bool=false)
    # Just like the function below, this produces 
    #hermitian an symmetric matrices, they are
    # just not annotated as such
    Q, _ = qr(rand(complex ? ComplexF64 : Float64, n,n))
    Q = Matrix(Q)
    
    #this is type unstable, since no Seed is given however
    #is probably desirable
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
* `A``: Vector of matrices with size ``n × n × k``

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
* `A`: Vector of matrices

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
* `A`: Vector of matrices

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
* `A`: AbstractMatrix of dimension ``n × n``
* `B`: AbstractMatrix of dimension ``n × n``

Check if two matrices A, B are commuting.
 ``AB = BA`` must hold.
 """
function is_commuting(A::AbstractMatrix, B::AbstractMatrix)
    return isapprox(A*B, B*A)
end


"""
    get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)
* `D`: Diagonal Matrix with offdiagonal elements set to zero
* `i,j`: Denotes the indexes of matrix D

Calculates the factor ``z_{ij}`` which is defined by: `` ∑_{k} D_{i,i}^{k}D_{j,j}^{k} ``
"""
function get_z_fdiag(D::AbstractArray{<:Number}, i::Int, j::Int)
    return sum(@view(D[i,i,:]).*@view(D[j,j,:]))
end

"""
    get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int, j::Int)
* `D`: Diagonal Matrix with offdiagonal elements set to zero
* `E`: Diagonal Matrix with diagonal elements set to zero
* `i,j`: Denotes the indexes of the matrices D and E

Calculates the factor ``y_{ij}`` which is defined by:
`` ∑_{k} D_{j,j}^{k}E_{j,i}^{k} ``
"""
function get_y_fdiag(D::AbstractArray{<:Number}, E::AbstractArray{<:Number}, i::Int,j::Int)
    return sum(@view(D[j,j,:]).*@view(E[i,j,:]))
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

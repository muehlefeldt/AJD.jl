using Base:OneTo

function rotation(aii::Array{T}, ajj::Array{T}, aij::Array{T}, aji::Array{T})::Matrix{T} where {T<:Complex}
    h = hcat(aii .- ajj, aij .+ aji, (aji .- aij) .* 1im)
    G = real(h' * h)
    _, vecs = eigen(G)
    x, y, z = vecs[:, end]
    if x < 0.0
        x, y, z = -x, -y, -z
    end
    r = sqrt(x^2 + y^2 + z^2)  # julia's eigen returns normalized eigenvectors
    @assert isapprox(r, 1.0)
    c = sqrt((x + 1.0) / 2.0)
    s = (y - z * im) / sqrt(2.0 * (x + 1.0))
    return [c conj(s); -s conj(c)]
end
#shouldn't rotation_symmetric only apply to matrices of type real?
function rotation_symmetric(aii::Array{T}, ajj::Array{T}, aij::Array{T})::Matrix{T} where {T<:Union{Real,Complex}}
    h = hcat(aii .- ajj, 2.0 .* aij)
    G = real(h' * h)
    _, vecs = eigen(G)
    x, y = vecs[:, end]
    if x < 0.0
        x, y = -x, -y
    end
    r = sqrt(x^2 + y^2)  # julia's eigen returns normalized eigenvectors
    @assert isapprox(r, 1.0)
    c = sqrt((x + 1.0) / 2.0)
    s = y / sqrt(2.0 * (x + 1.0))
    return [c s; -s c]
end

"""
    jdiag_edourdpineau(X::Vector{M}; iter=100, eps=1e-3) 
        where {T<:Union{Real,Complex},M<:AbstractMatrix{T}}

Diagonalize a set of matrices using the Jacobi method ("Jacobi Angles for Simultaneous Diagonalization").
Code adapted from [Edouardpineaus Python implementation](https://github.com/edouardpineau/Time-Series-ICA-with-SOBI-Jacobi)
"""
function jdiag_edourdpineau(X::Vector{M}; iter=100, eps=1e-3) where {T<:Number,M<:AbstractMatrix{T}}
    
    Xm = cat(X..., dims=3) .+ 0.0im # change to Xm = complex.(cat(X..., dims = 3))? -NG
    m = length(X)
    n = size(X[1], 1)
    @assert n == size(X[1], 2)

    if !(M <: Symmetric) && (T <: Real)
        Xm .+= 0.0im # is this necessary? we already did that in line 51 no? - NG
        V = Matrix{Complex{T}}(I, n, n)
    else
        V = Matrix{T}(I, n, n)
    end

    diag_err = frobenius_offdiag_norm(Xm)
    err_array = [diag_err]
    diff = Inf
    current_iter = 0
    while (diff > eps) && (current_iter < iter)
        current_iter += 1
        # println("Current iteration: ", current_iter, " with error: ", diag_err)
        for i in OneTo(n - 1), j in (i+1):n

            if M <: Symmetric || M <: Hermitian
                R = rotation_symmetric(Xm[i, i, :], Xm[j, j, :], Xm[i, j, :])
            else
                R = rotation(Xm[i, i, :], Xm[j, j, :], Xm[i, j, :], Xm[j, i, :])
            end

            for k in OneTo(m)
                Xm[[i, j], :, k] = R * Xm[[i, j], :, k]
                Xm[:, [i, j], k] = Xm[:, [i, j], k] * R'
            end
            V[:, [i, j]] = V[:, [i, j]] * R'
        end

        new_diag_err = frobenius_offdiag_norm(Xm)
        push!(err_array, new_diag_err)
        # TODO: this is relative error, should also add an absolute error
        diff = abs(new_diag_err - diag_err) / diag_err
        diag_err = new_diag_err
    end
    return V, Xm, err_array #added return statement for clearer code-NG
end

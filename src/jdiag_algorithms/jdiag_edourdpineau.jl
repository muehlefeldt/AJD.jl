
function rotation(
    aii::Array{T},
    ajj::Array{T},
    aij::Array{T},
    aji::Array{T},
)::Matrix{T} where {T<:Complex}
    h = hcat(aii .- ajj, aij .+ aji, (aji .- aij) .* 1im)
    G = real(h' * h)
    _, vecs = eigen(G)
    x, y, z = vecs[:, end]
    if x < 0.0
        # we use abs on x so julia knows that the values in the sqrt are positive
        x, y, z = -x, -y, -z
    end
    c = sqrt(abs((x + 1.0) / 2.0))
    s = (y - z * im) / sqrt(abs(2.0 * (x + 1.0)))
    return [c conj(s); -s conj(c)]
end


function rotation_symmetric(
    aii::Array{T},
    ajj::Array{T},
    aij::Array{T},
)::Matrix{Real} where {T<:Union{Real,Complex}}
    h = hcat(aii .- ajj, 2.0 .* aij)
    G = real(h' * h)
    # G is now a 2x2 symmetric matrix
    a, b = G[1, 1], G[1, 2]
    # Since G is symmetric, G[2,1] = G[1,2] and G[2,2] = c
    c = G[2, 2]

    # Calculate theta for the largest eigenvector
    theta = 0.5 * atan(2b, a - c)
    x, y = cos(theta), sin(theta)

    # x will always be positive, we add the abs so the compiler know this too
    c = sqrt(abs((x + 1.0) / 2.0))
    s = y / sqrt(abs(2.0 * (x + 1.0)))
    return [c s; -s c]
end

"""
    jdiag_edourdpineau(X::Vector{M}; iter=100, eps=1e-3)
        where {T<:Union{Real,Complex},M<:AbstractMatrix{T}}

Diagonalize a set of matrices using the Jacobi method ("Jacobi Angles for Simultaneous Diagonalization").
Code adapted from [Edouardpineaus Python implementation](https://github.com/edouardpineau/Time-Series-ICA-with-SOBI-Jacobi)
"""
function jdiag_edourdpineau(
    X::Vector{M};
    iter = 100,
    rtol = 1e-3,
    atol = eps(),
) where {T<:Number,M<:AbstractMatrix{T}}

    Xm = cat(X..., dims = 3)
    m = length(X)
    n = size(X[1], 1)
    @assert n == size(X[1], 2)

    if !(M <: Symmetric) && (T <: Real)
        Xm = complex.(Xm)
        V = Matrix{Complex{T}}(I, n, n)
    else
        V = Matrix{T}(I, n, n)
    end

    norm = frobenius_offdiag_norm(Xm)
    norm_history = [norm]

    # Initial setup of the progressbar.
    progress_bar = ProgressThresh(atol; desc="Minimizing:")

    # Initial setup of the progressbar.
    progress_bar = ProgressThresh(atol; desc="Minimizing:")

    # Iteration counter.
    n_iteration = 0
    
    for _ = 1:iter
        n_iteration += 1
        for i = 1:(n-1), j = (i+1):n
            if M <: Symmetric
                R = rotation_symmetric(Xm[i, i, :], Xm[j, j, :], Xm[i, j, :])
            else
                R = rotation(Xm[i, i, :], Xm[j, j, :], Xm[i, j, :], Xm[j, i, :])
            end

            for k = 1:m
                Xm[[i, j], :, k] = R * Xm[[i, j], :, k]
                Xm[:, [i, j], k] = Xm[:, [i, j], k] * R'
            end
            V[:, [i, j]] = V[:, [i, j]] * R'
        end

        new_norm = frobenius_offdiag_norm(Xm)
        push!(norm_history, new_norm)

        diff = abs(new_norm - norm)

        # Update progress info.
        update!(progress_bar, diff)

        if diff < atol || diff < rtol * norm
            break
        end
        norm = new_norm
    end
    finish!(progress_bar)
    return V, Xm, norm_history, n_iteration
end

module JDiag
using LinearAlgebra
export testJDiag

function off_diag_norm(Xm::AbstractArray{T,3})::Real where {T}
  sum = zero(real(T))
  for k in axes(Xm, 1), i in axes(Xm, 2), j in axes(Xm, 3)
    i == j && continue
    sum += abs2(Xm[k, i, j])
  end
  return sum
end

function rotation(aii::Array{T}, ajj::Array{T}, aij::Array{T}, aji::Array{T})::Matrix{Complex{T}} where {T<:Union{Real,Complex}}
  # println("Matrix rotation")

  h = hcat(aii .- ajj, aij .* aji, (aji .- aij) .* im)
  G = real(transpose(h) * h)
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

function rotation_symmetric(aii::Array{T}, ajj::Array{T}, aij::Array{T})::Matrix{T} where {T<:Union{Real,Complex}}
  # println("Symmetric rotation")

  h = hcat(aii .- ajj, 2.0 .* aij)
  G = real(transpose(h) * h)
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

function testJDiag(X::Vector{M}; iter=100, eps=1e-3) where {T<:Union{Real,Complex},M<:AbstractMatrix{T}}

  Xm = cat(X..., dims=3)

  m = length(X)
  n = size(X[1], 1)
  k = size(X[1], 2)
  @assert n == k

  V = Matrix{T}(I, n, n)
  println(V)

  diag_err = off_diag_norm(Xm)
  err_array = [diag_err]
  diff = Inf
  current_iter = 0
  while (diff > eps) && (current_iter < iter)
    current_iter += 1
    println("Current iteration: ", current_iter, " with error: ", diag_err)
    for i in Base.OneTo(n - 1), j in (i+1):n
      R = M <: Symmetric ? rotation_symmetric(Xm[i, i, :], Xm[j, j, :], Xm[i, j, :]) :
          rotation(Xm[i, i, :], Xm[j, j, :], Xm[i, j, :], Xm[j, i, :])
      Xm[[i, j], :, :] = reshape(R * reshape(Xm[[i, j], :, :], 2, :), 2, n, m)
      Xm[:, [i, j], :] = permutedims(reshape(reshape(permutedims(Xm[:, [i, j], :], (1, 3, 2)), :, 2) * R', n, m, 2), (1, 3, 2))
      V[:, [i, j]] = V[:, [i, j]] * R'
    end
    new_diag_err = off_diag_norm(Xm)
    push!(err_array, new_diag_err)
    diff = abs(new_diag_err - diag_err)
    diag_err = new_diag_err
  end

  V, Xm, err_array
end
end

Xsset = [Symmetric(rand(Float64, 3, 3)) for _ in 1:3]
Xset = [rand(Float64, 5, 5) for _ in 1:3]

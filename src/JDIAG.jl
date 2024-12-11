module JDiag
using LinearAlgebra
export testJDiag

function rotation(aii::Array{T}, ajj::Array{T}, aij::Array{T}, aji::Array{T})::Matrix{Complex{T}} where {T<:Union{Real,Complex}}
    println("Matrix rotation")

    h = hcat(aii .- ajj, aij * aji, (aji .- aij) .* im)
    G = real(transpose(h) * h)
    _, vecs = eigen(G)
    x, y, z = vecs[:, end]
    r = sqrt(x^2 + y^2 + z^2)
    c = sqrt((x+r)/(2r))
    s = (y - z*im) / sqrt(2r*(x+r))
    return [c conj(s);
            -s conj(c)]
end

function rotation_symmetric(aii::Array{T}, ajj::Array{T}, aij::Array{T})::Matrix{T} where {T<:Union{Real,Complex}}
    println("Symmetric rotation")

    h = hcat(aii .- ajj, 2 .* aij)
    G = real(transpose(h) * h)
    _, vecs = eigen(G)
    x, y = vecs[:, end]
    r = sqrt(x^2 + y^2)
    c = sqrt((x+r)/(2r))
    s = y / sqrt(2r*(x+r))
    return [c s;
            -s c]
end

function testJDiag(X::Vector{M})::T where {T<:Union{Real,Complex}, M<:AbstractMatrix{T}}
    Xm = cat(X..., dims=3)
    m = length(X)
    n = size(X[1], 1)
    println("m: ", m)
    println("n: ", n)
    i = 1
    j = 2
    if M <: Symmetric
        @assert Xm[i,j,:] == Xm[j,i,:]
        rotation_symmetric(Xm[i,i,:],Xm[j,j,:],Xm[i,j,:])
    else
        rotation(Xm[i,i,:],Xm[j,j,:],Xm[i,j,:],Xm[j,i,:])
    end

    1.0
end

end

Xsset = [Symmetric(rand(Float64, 5,5)) for _ in 1:3]
Xset = [rand(Float64, 5,5) for _ in 1:3]

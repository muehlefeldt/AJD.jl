include("AJD.jl")
using .AJD
using Plots
using LinearAlgebra

all_errors = Vector{Vector{Float64}}()
#
m = 200  # number of matrices
n = 4  # matrix size
n_runs = 10

for i in 1:10
    # First we generate a random orthonormal matrix Q
    P = rand(n,n)
    Q, _ = qr(P)
    Q = Matrix(Q)
    # 
    D_list = [Diagonal(rand(n)) for i in 1:m]
    A = [Q*D*Q' for D in D_list]

    V, _, err = AJD.jdiag_edourdpineau(A)
    norm(V)
    push!(all_errors, err)
end

p = plot(all_errors[5], yaxis=:log, xlabel="Iteration", ylabel="Error",
         title="Convergence of Jacobi method", label="Run 1", alpha=0.7)

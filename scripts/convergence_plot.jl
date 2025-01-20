include("AJD.jl")
using .AJD
using Plots
using LinearAlgebra


symm = false
m = 15  # number of matrices
n = 10  # matrix size
n_runs = 10

all_errors = Vector{Vector{Float64}}()
for i in 1:10
    if symm
        A = [Symmetric(rand(n,n)) for i in 1:m]
    else
        A = [rand(n,n) for i in 1:m]
    end

    _, _, err = AJD.jdiag_edourdpineau(A)
    push!(all_errors, err)

    # len = length(err)
    # s = min((len+1) % 8, len)

    plot(err, yaxis=:log, xlabel="Iteration", ylabel="Error", title="Convergence of Jacobi method")
end


# Add remaining trajectories
p = plot(all_errors[5], yaxis=:log, xlabel="Iteration", ylabel="Error",
         title="Convergence of Jacobi method", label="Run 1", alpha=0.7)

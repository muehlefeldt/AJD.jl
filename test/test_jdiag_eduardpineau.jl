using AJD

include("test_utils.jl")

@testset "jdiag_exact" begin
    n = 10
    m = 20
    X = random_normal_commuting_matrices(n, m)
    V, Xnew, e = AJD.jdiag_edourdpineau(X)
    @assert isreal(V)
    @assert e[end] < 1e-15
    @assert V'*V ≈ I
    for k in 1:m
        @assert V*Xnew[:,:,k]*V' ≈ X[k]
    end
end



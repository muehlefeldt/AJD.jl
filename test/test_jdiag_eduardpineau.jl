
@testset "jdiag_exact" begin
    n = 10
    m = 20
    X = AJD.random_normal_commuting_matrices(n, m)
    V, Xnew, e = AJD.jdiag_edourdpineau(X)
    @assert isreal(V)
    @assert e[end] < 1e-15
    @assert V'*V ≈ I
    for k in 1:m
        @assert V*Xnew[:,:,k]*V' ≈ X[k]
    end
end
@testset "jdiag_edourdpineau" begin
    testinput = [[ 1.0 0.0 1.0*im; 0.0 2.0 0.0; 1.0*im 0.0 1.0],[ 1.0 0.0 1.0*im; 0.0 2.0 0.0; 1.0*im 0.0 1.0]]
    @info diagonalize(testinput, algorithm = "jdiag_edourdpineau")
end


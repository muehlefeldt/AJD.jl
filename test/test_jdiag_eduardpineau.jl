# Test implementation of AJD.jdiag_edourdpineau().

"Accepted error level."
accepted_error = 1e-15

@testset "AJD.jdiag_edourdpineau() Basic" begin
    n = 10
    m = 20
    X = AJD.random_normal_commuting_matrices(n, m)
    V, Xnew, e = AJD.jdiag_edourdpineau(X)
    @test isreal(V)
    @test e[end] < accepted_error
    @test V' * V ≈ I

    for k in 1:m
        @test V * Xnew[:,:,k] * V' ≈ X[k]
    end
end

# Single test case to verify complex matrices handling.
# Checks nonDiagonality of F' * A * F with F being the filter and A the complex matrix.
# Compare to test_nondiagonality.jl.
@testset "AJD.jdiag_edourdpineau() Complex" begin
    test_input = [[ 1.0 0.0 1.0*im; 0.0 2.0 0.0; 1.0*im 0.0 1.0],[ 1.0 0.0 1.0*im; 0.0 2.0 0.0; 1.0*im 0.0 1.0]]
    result = diagonalize(test_input, algorithm = AJD.JDiagEdourdPineau())

    @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
end

# Test the nonDiagonality of the algorithm outputs when working on
# time correlated covriance matrices.

# Base testset.
@testset "Nondiagonality Time Correlated" begin
    # Define the acceptable error level.
    # How far the diagonalised matrices can be away from a perfect diagonal matrix.
    accepted_error = 0.8
    for alg in AJD.ALL_ALGORITHMS
        test_input = AJD.get_test_data(:approx_diag,"")
        result = diagonalize(test_input, algorithm=alg)
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end
end

@testset "Nondiagonality Time Correlated Large" begin
    accepted_error = 0.8
    for alg in AJD.ALL_ALGORITHMS
        test_input = AJD.get_test_data(:approx_diag_large,"")
        result = diagonalize(test_input, algorithm=alg)
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end
end

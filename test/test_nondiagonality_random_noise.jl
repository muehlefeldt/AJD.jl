@testset "Nondiagonality Random Noise" begin
    # Define the acceptable error level.
    # How far the diagonalised matrices can be away from a perfect diagonal matrix.
    accepted_error = 0.5
    for name in AJD.ALL_ALGORITHMS
        test_input = AJD.get_test_data(:random_noise)
        result = diagonalize(test_input, algorithm=name)
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end
end
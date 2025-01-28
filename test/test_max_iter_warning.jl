# Test warning on reaching max_iter during calculations.

# Warning issued on reaching the max iteration count.
# Assumption is that further iterations may yield better result. 
# As such warning needs to be issued.
@testset "Warning max_iter reached" begin
    for name in AJD.ALL_ALGORITHMS
        test_input = AJD.random_normal_commuting_matrices(20, 6)
        @test_logs (
            :warn,
            "Max iteration was reached. Consider increasing max_iter: diagonalize(M, max_iter=...).",
        ) diagonalize(test_input, algorithm = name, max_iter = 1)
    end
end

# Error thrown if invalid iteration count requested by the user.
@testset "Error max_iter too low" begin
    for name in AJD.ALL_ALGORITHMS
        test_input = AJD.random_normal_commuting_matrices(20, 6)
        @test_throws ArgumentError diagonalize(test_input, algorithm = name, max_iter = 0)
        @test_throws ArgumentError diagonalize(test_input, algorithm = name, max_iter = -1)
    end
end
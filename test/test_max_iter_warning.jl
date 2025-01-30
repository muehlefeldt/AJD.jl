# Test warning on reaching max_iter during calculations.

# Warning issued on reaching the max iteration count.
# Assumption is that further iterations may yield better result.
# As such warning needs to be issued.
@testset "Warning max_iter reached" begin
    for alg in AJD.ALL_ALGORITHMS
        test_input = AJD.random_normal_commuting_matrices(20, 6)
        @test_logs (
            :warn,
            "Max iteration was reached. Consider increasing max_iter: diagonalize(M, max_iter=...).",
        ) diagonalize(test_input, algorithm = alg, max_iter = 1)
    end
end

# Error thrown if invalid iteration count requested by the user.
@testset "Error max_iter too low" begin
    for alg in AJD.ALL_ALGORITHMS
        test_input = AJD.random_normal_commuting_matrices(20, 6)
        @test_throws ArgumentError diagonalize(
            test_input,
            algorithm = alg,
            max_iter = 0,
        )
        @test_throws ArgumentError diagonalize(
            test_input,
            algorithm = alg,
            max_iter = -1,
        )
    end
end

# Make sure the counting of the iterations is consitent between the algorithms.
@testset "Number of iterations correct" begin
    for alg in AJD.ALL_ALGORITHMS
        test_input = AJD.random_normal_commuting_matrices(20, 6)
        requested_iter = 1
        _, _, _, n_iter = AJD.get_diagonalization(
            test_input,
            algorithm = alg,
            max_iter = requested_iter,
        )
        @test n_iter == requested_iter

        requested_iter = 2
        _, _, _, n_iter = AJD.get_diagonalization(
            test_input,
            algorithm = alg,
            max_iter = requested_iter,
        )
        @test n_iter == requested_iter
    end
end

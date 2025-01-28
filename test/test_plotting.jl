# Test the plotting of diagonalized matrices and convergence behaviour.
# Verify all plot combinations and the return type.

using Plots: Plot

accepted_error = 0.1

@testset "Plotting" begin
    # Check if the normal operation is still performed even if plot of filter and
    # input matrices is selected.
    # Only using very small inputs.
    # Check all implemented algorithms.
    for name in AJD.ALL_ALGORITHMS
        test_input = AJD.random_normal_commuting_matrices(10, 6)
        result = diagonalize(test_input, :plot, algorithm=name)
        @test typeof(result) <: Plot

        test_input = AJD.random_normal_commuting_matrices(10, 6)
        @test_throws ArgumentError diagonalize(test_input, :test, algorithm=name)

        # Verify get_diagonalization() works as well when plot is requested.
        # Expected are the returns:
        # Filter, diagonalized matrices and the error array.
        result = AJD.get_diagonalization(test_input, algorithm=name, only_plot=:plot)
        @test length(result) == 4
        @test length(result[3]) > 0
        @test mean([nonDiagonality((result[1]') * A * (result[1])) for A in test_input]) < accepted_error
    end
end

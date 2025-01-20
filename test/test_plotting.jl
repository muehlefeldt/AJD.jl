# Test the plotting of diagonalized matrices and convergence behaviour.
# Verify all plot combinations and the return type.

using Plots: Plot
using PosDefManifold: mean
using AJD

accepted_error = 0.1

@testset "Plotting" begin
    # Check if the normal operation is still performed even if plot of filter and
    # input matrices is selected.
    # Only using very small inputs.
    # Check all implemented algorithms.
    for name in AJD.ALL_ALGORITHMS
        test_input = AJD.random_normal_commuting_matrices(10, 6)
        result = diagonalize(test_input, algorithm=name, plot_matrix=true)
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error

        # Check if the normal operation is still performed even if plot of filter and
        # input matrices is selected in addition to convergence plot.
        test_input = AJD.random_normal_commuting_matrices(10, 6)
        result = diagonalize(test_input, algorithm=name, plot_matrix=true, plot_convergence=true)
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error

        # Only convergence plot selected.
        test_input = AJD.random_normal_commuting_matrices(10, 6)
        result = diagonalize(test_input, algorithm=name, plot_convergence=true)
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end

    # Directly test plot functions to generate lineplot and heatmaps with complex matrices.
    # For now only the return type can be checked.
    test_input = AJD.random_normal_commuting_matrices(10, 6)
    # F, B are complex.
    F, B, error_array = AJD.jdiag_edourdpineau(test_input)
    @test typeof(AJD.plot_convergence_lineplot(error_array, "test")) <: Plot
    @test typeof(AJD.plot_matrix_heatmap(F, B)) <: Plot
    @test length(error_array) > 2
    
    # Directly test plot functions to generate lineplot and heatmaps.
    # For now only the return type can be checked.
    test_input = AJD.random_normal_commuting_matrices(10, 6)
    F, B, error_array = AJD.jdiag_cardoso(test_input, 1e-6, plot_convergence=true)
    @test typeof(AJD.plot_convergence_lineplot(error_array, "test")) <: Plot
    @test typeof(AJD.plot_matrix_heatmap(F, B)) <: Plot
    @test length(error_array) > 2
end

# Test the plotting of diagonalize.

@testset "Plotting" begin
    test_input = AJD.random_normal_commuting_matrices(10, 6)
    result = diagonalize(test_input, plot_matrix=true)
    #@test
end
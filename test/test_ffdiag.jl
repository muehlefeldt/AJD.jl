# FFDiag specific tests.

# using Diagonalizations: nonDiagonality
# using PosDefManifold: mean

accepted_error = 1e-6

@testset "FFdiag Functionality" begin
    # Test default FFDiag implementation.
    test_input = AJD.random_normal_commuting_matrices(3, 2)
    result = diagonalize(test_input, algorithm="FFD")
    @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    
    # Test use of infinity norm. Default always uses frobenius norm.
    result = AJD.create_linear_filter(AJD.ffd(test_input, norm_ = "inf")[1])
    @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
end
using LinearAlgebra
using Diagonalizations
using PosDefManifold

accepted_error = 1e-6

@testset "FFdiag Functionality" begin
    test_input = AJD.random_normal_commuting_matrices(3, 2)
    result = diagonalize(test_input, algorithm="FFD")
    @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    #@info result
    #@info test_input
    #@info mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) 
    #@info "FFDiag",diagonalize(test_input, algorithm = "FFD")
    #test_input = [Matrix(1.0I,3,3), Matrix(1.0I,3,3)]
    #result = diagonalize(test_input, algorithm="FFD")
    #@info result
    #@info test_input
    #@info mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) 
end
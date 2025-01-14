using LinearAlgebra
using Diagonalizations
using PosDefManifold
@testset "ffdiag_functionality" begin
    test_input = AJD.random_normal_commuting_matrices(3, 2; complex=false)
    result = diagonalize(test_input, algorithm="FFD")
    @info result
    @info test_input
    @info mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) 
    #@info "FFDiag",diagonalize(test_input, algorithm = "FFD")
    test_input = [Matrix(1.0I,3,3), Matrix(1.0I,3,3)]
    result = diagonalize(test_input, algorithm="FFD")
    @info result
    @info test_input
    @info mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) 
end
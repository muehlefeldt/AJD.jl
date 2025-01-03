using AJD
using Diagonalizations

# Define the acceptable error level.
# How far the diagonalised matrices can be away from a perfect diagonal matrix.
const accepted_error = 1e-6

@testset "Nondiagonality" begin
    for name in ["jdiag_gabrieldernbach", "jdiag_cardoso", "jdiag_edourdpineau"]
        test_input = AJD.random_normal_commuting_matrices(10, 6)
        result = diagonalize(test_input, algorithm=name)
        
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end
end

@testset "Nondiagonality Complex" begin
    for name in ["jdiag_edourdpineau", "jdiag_gabrieldernbach"]
        test_input = AJD.random_normal_commuting_matrices(10, 6; complex=true)
        result = diagonalize(test_input, algorithm=name)
        
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end
end
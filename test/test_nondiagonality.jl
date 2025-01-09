# Testing the Nondiagonality of F' * A * F.
# F being the calculated filter and A being an input matrix.
# The result nonDiagonality() indicates the "distance" from a perfect diagonal matrix.

using AJD
using Diagonalizations
using PosDefManifold
# Define the acceptable error level.
# How far the diagonalised matrices can be away from a perfect diagonal matrix.
# Diagonalizations.jl set the error level to 1e-6.
const accepted_error = 1e-6
#"jdiag_cardoso"
@testset "Nondiagonality" begin
    for name in ["jdiag_gabrieldernbach" ,"jdiag_edourdpineau"]
        test_input = AJD.random_normal_commuting_matrices(10, 6)
        result = diagonalize(test_input, algorithm=name)
        
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end
end

@testset "Nondiagonality Complex" begin
    # Select algorithms supporting complex matrices.
    for name in ["jdiag_edourdpineau", "jdiag_gabrieldernbach"]
        # Generate complex test input and calculate filter.
        test_input = AJD.random_normal_commuting_matrices(6, 6; complex=true)
        result = diagonalize(test_input, algorithm=name)
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end
end
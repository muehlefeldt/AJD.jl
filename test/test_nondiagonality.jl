# Testing the Nondiagonality of F' * A * F.
# F being the calculated filter and A being an input matrix.
# The result of nonDiagonality() indicates the "distance" from a perfect diagonal matrix.
# nonDiagonality(I) = 0 as I is a diagonal matrix.

using AJD
using Diagonalizations
using PosDefManifold
# Define the acceptable error level.
# How far the diagonalised matrices can be away from a perfect diagonal matrix.
# Diagonalizations.jl set the error level to 1e-6.
accepted_error = 1e-6

# Base testset with real matrices.
@testset "Nondiagonality" begin
    for name in ["jdiag_gabrieldernbach" ,"jdiag_edourdpineau", "FFD"]
        test_input = AJD.random_normal_commuting_matrices(10, 6)
        result = diagonalize(test_input, algorithm=name)
        # Cardoso implementation shows very high error level.
        if name == "jdiag_cardoso"
            @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < 0.1
        else 
            @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
        end
    end
end

@testset "Nondiagonality Symmetric" begin
    # jdiag_edourdpineau has a function that is overloaded for symmetric matrices.
    for name in ["jdiag_edourdpineau"]
        test_input = AJD.random_normal_commuting_symmetric_matrices(10, 6)
        result = diagonalize(test_input, algorithm=name)
        # Cardoso implementation shows very high error level.
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end
end

# Test the algorithms with a single matrix as input.
@testset "Nondiagonality Single Matrix" begin
    for name in ["jdiag_gabrieldernbach", "jdiag_cardoso", "jdiag_edourdpineau"]
        test_input = AJD.random_normal_commuting_matrices(10, 1)
        result = diagonalize(test_input, algorithm=name)
        # TODO: Cardoso implementation shows very high error level.
        if name == "jdiag_cardoso"
            @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < 0.1
        else 
            @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
        end
    end
end

# Test nonDiagonality of complex inputs.
@testset "Nondiagonality Complex" begin
    # Select algorithms supporting complex matrices.
    for name in ["jdiag_edourdpineau", "jdiag_gabrieldernbach"]
        # Generate complex test input and calculate filter.
        test_input = AJD.random_normal_commuting_matrices(6, 6; complex=true)
        result = diagonalize(test_input, algorithm=name)
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end

    # Tests to check combination of real and complex matrices.
    for name in ["jdiag_edourdpineau", "jdiag_gabrieldernbach"]
        #A = AJD.random_normal_commuting_matrices(10, 6)
        #B = AJD.random_normal_commuting_matrices(10, 1; complex=true)
        #test_input = [A..., B...]
        test_input = AJD.get_test_data_complex_real(10, 10)
        result = diagonalize(test_input, algorithm=name)
        
        # TODO: Error level very high with both implementations.
        # Surprisingly was jdiag_edourdpineau few days ago.
        @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
    end
end
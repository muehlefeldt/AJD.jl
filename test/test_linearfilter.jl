# Test the output of all algorithms for LinearFilter.

@testset "LinearFilter" begin
    # Iterate over all implemented algorithms.
    for name in ["jdiag_gabrieldernbach", "jdiag_cardoso", "jdiag_edourdpineau"]
        test_input = AJD.random_normal_commuting_matrices(10, 6)
        result = diagonalize(test_input, algorithm=name)
        
        # Check for LinearFilter.
        @test typeof(result) == LinearFilter

        # Check Linearfilter contents.
        @test typeof(result.F) <: AbstractMatrix{T} where {T<:Number}
        @test typeof(result.iF) <: AbstractMatrix{T} where {T<:Number}
    end
end
# Verify an ArgumentError is thrown when algorithm = "jdiag_cardoso" is used with
# complex matrices.
@testset "AJD.jdiag_cardoso() No Complex" begin
    test_input = (1.0) * [Matrix(I, 6, 6), Matrix(I, 6, 6)]
    @test_throws ArgumentError diagonalize(test_input*im, algorithm = "jdiag_cardoso")
end
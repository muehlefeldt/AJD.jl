#check if int conversion in AJD works
#could do a check for all algorithms however due to 
#constraints of FFDiag
#only edouardpineau is tested
accepted_error = 1e-6
@testset "type_conversion" begin
    test_input = 1*[I(3),I(3)]
    result = diagonalize(test_input)
    @test mean([nonDiagonality(result.iF * A * result.F) for A in test_input]) < accepted_error
end
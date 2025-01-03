using PosDefManifold
using LinearAlgebra: opnorm
using AJD: iscommuting, issamesize, isstrictly_diagonally_dominant, get_offdiag_elements, get_diag_elements,frobenius_offdiag_normation, sort_offdiag_elements

@testset "iscommuting" begin
    A = [1 1; 1 1]
    C = [0 2; 1 2]
    @test iscommuting(A,A) == true
    @test iscommuting(A,C) == false
end

@testset "issamesize" begin
    A = ones(2,2)
    B = ones(1,2)
    @test issamesize(A,A) == true
    @test issamesize(A,B) == false
end
@testset "strictly_dominant" begin
    A = 1.0*Matrix(I,3,3)
    @test isstrictly_diagonally_dominant(A) == true
end
@testset "get_off_diag_elements" begin
    input = [[1 2 3; 4 5 6; 7 8 9];;;[0 2 4; 4 2 1; 9 0 1]]
    @test get_offdiag_elements(input) == [[0 2 3; 4 0 6; 7 8 0];;; [0 2 4; 4 0 1; 9 0 0]]
    input = reshape(repeat(Matrix(1.0I,3,3), outer = (1,2)),3,3,2)
    @test get_offdiag_elements(input) == zeros(size(input))
end
@testset "frobenius_off_diag_normation" begin
    input = [[1 2 3; 4 5 6; 7 8 9];;;[0 2 4; 4 2 1; 9 0 1]]
    @test frobenius_offdiag_normation(input) == 296
    input = reshape(repeat(Matrix(1.0I,3,3), outer = (1,2)),3,3,2)
    @test frobenius_offdiag_normation(input) == 0
end
@testset "get_diag_elements" begin
    input = reshape(repeat(Matrix(1.0I,3,3), outer = (1,2)),3,3,2)
    @test get_diag_elements(input) == input
end
@testset "sort_offdiag_elements" begin
    input_1D = [1 2 3; 4 5 6; 7 8 9]
    @test sort_offdiag_elements(input_1D) == float.([2,4,3,7,6,8])
    input = [[1 2 3; 4 5 6; 7 8 9];;;[0 2 4; 4 2 1; 9 0 1]]
    @test sum(sort_offdiag_elements(input)) == sum(get_offdiag_elements(input))
    input = reshape(repeat(Matrix(1.0I,2,2)))
end
#TODO: Old code, if not needed should be deleted and the functions as well
# @testset "inf_off_diag_normation" begin
#     W = [1 2 3; 4 5 6; 7 8 9]
#     @test inf_off_diag_normation(W) == 15
# end

# @testset "issymmetric" begin
#     A = ones(2,2)
#     B = ones(2,1)
#     @test issymmetric(A) == true
#     @test issymmetric(B) == false
# end

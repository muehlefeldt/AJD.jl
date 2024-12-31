using PosDefManifold
using LinearAlgebra
using AJD: iscommuting, issamesize, isstrictly_diagonally_dominant

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

# @testset "issymmetric" begin
#     A = ones(2,2)
#     B = ones(2,1)
#     @test issymmetric(A) == true
#     @test issymmetric(B) == false
# end

# @testset "get_non_Diag_elements" begin
#     for i = 2:10
#         A = ones(i,i)
#         rows, columns = size(A)
#         @test get_non_Diag_elements(A) == reshape(ones(1,rows*columns - columns),rows*columns - columns)
#     end
# end

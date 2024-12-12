using AJD
using Test
using LinearAlgebra
using PosDefManifold
@testset "Is_Commuting" begin
        A = [1 1; 1 1]
        C = [0 2; 1 2]
        @test Is_Commuting(A,A) == true
        @test Is_Commuting(A,C) == false
end
@testset "Is_Same_size" begin
    A = ones(2,2)
    B = ones(1,2)
    @test Is_Same_size(A,A) == true
    @test Is_Same_size(A,B) == false
end
@testset "Is_Symmetric" begin
    A = ones(2,2)
    B = ones(2,1)
    @test Is_Symmetric(A) == true
    @test Is_Symmetric(B) == false
end
@testset "get_non_Diag_elements" begin
    for i = 2:10
        A = ones(i,i)
        rows, columns = size(A)
        @test get_non_Diag_elements(A) == reshape(ones(1,rows*columns - columns),rows*columns - columns)
    end
end
@testset "JADE" begin
    A = [[1 0; 0 1];;;[2 2; 2 2]]
    
    @info JADE(A)
end
#@testset "convert_Hermitian" begin
    
#end
#@testset "Jacobi_Rotation" begin
#    A = 1.0*I(3)
#    @test Jacobi_Rotation
#end
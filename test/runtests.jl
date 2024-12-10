using AJD
using Test

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

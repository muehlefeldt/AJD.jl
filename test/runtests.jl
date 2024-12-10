using AJD
using Test

@testset "Is_Commuting" begin
    @test Is_Commuting([1 0; 1 0], [1 0; 1 0]) == true
    @test Is_Commuting([1 2; 1 2], [2 2; 2 2]) == false
 end    

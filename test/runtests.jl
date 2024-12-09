using AJD
using Test
using Diagonalizations

@testset "AJD.jl" begin
    @testset "pythoncomparison" begin
        @test 1+2 == 1 + 2
    end 
    # Write your tests here.
    include("test_python_function.jl")
end

using AJD
using Test
using Diagonalizations

@testset "AJD.jl" begin
    # Write your tests here.
    include("test_python_function.jl")
    include("test_diagonalizations_against_python.jl")
end

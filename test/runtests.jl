using AJD
using Test
using Diagonalizations
using LinearAlgebra

@testset "AJD.jl" begin
    # Write your tests here.
    include("test_python_function.jl")
    include("test_diagonalizations_against_python.jl")
    include("test_jdiag_against_python.jl")
    include("test_jdiag_jade_algorithm.jl")
    include("test_jdiag_cardoso.jl")
end

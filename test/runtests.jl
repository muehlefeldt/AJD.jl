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
end



#@testset "Jacobi_Rotation" begin
#    A = 1.0*I(3)
#    @test Jacobi_Rotation
#end
using Test
using PyCall
using AJD

py"""
def py_add(x, y):
    return x + y

def py_multiply(x, y):
    return x * y
"""

# Example implementation of tests comparing Julia code to Python implementation.
@testset "pythoncomparison" begin
    @test py"py_add"(1, 2) == 1 + 2
    @test py"py_multiply"(3, 4) == multiply(3, 4)
end

using Test
using PyCall

py"""
def py_add(x, y):
    return x + y
"""

@testset "pythoncomparison" begin
    @test py"py_add"(1, 2) == 1 + 2
end 

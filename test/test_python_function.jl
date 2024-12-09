using Test
using PyCall

py"""
    def add(x, y):
        return a + y
"""
@testset begin
    @test py"add"(1, 2) == 1 + 2
end 

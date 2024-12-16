using AJD
using Test
using Diagonalizations
using LinearAlgebra

@testset "AJD.jl" begin
    include("test_jdiag_gabrieldernbach.jl")
    include("test_jdiag_gabrieldernbach_utils.jl")
    include("test_jdiag_cardoso.jl")
    include("test_jdiag_eduardpineau.jl")
end

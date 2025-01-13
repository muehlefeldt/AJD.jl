using AJD
using Test
using Diagonalizations
using LinearAlgebra
@testset "AJD.jl" begin
    #include("test_jdiag_gabrieldernbach.jl")
    include("test_jdiag_gabrieldernbach_utils.jl")
    #include("test_jdiag_cardoso.jl")
    #include("test_jdiag_eduardpineau.jl")
    #include("test_linearfilter.jl")
    #include("test_nondiagonality.jl")
    #include("test_utils.jl")
    #include("test_benchmark.jl")
    #include("test_ffdiag.jl")
end

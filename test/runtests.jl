using AJD
using Test
# using Diagonalizations
using LinearAlgebra
using Diagonalizations: nonDiagonality,LinearFilter
using PosDefManifold: mean
using BenchmarkTools: BenchmarkGroup
# using PyCall
using Plots: Plot

@testset "AJD.jl" begin
    # Removed test_jdiag_gabrieldernbach_utils.jl due to an unresolved issue with PyCall
    # and the use in Github CI test runs.
    include("test_jdiag_gabrieldernbach.jl")

    include("test_jdiag_cardoso.jl")
    include("test_jdiag_eduardpineau.jl")
    include("test_linearfilter.jl")

    # Test the actual result of the algorithm by cheking how nondiagonal the matrices are.
    include("test_nondiagonality.jl")
    include("test_nondiagonality_time_correlated.jl")

    include("test_utils.jl")
    include("test_benchmark.jl")
    include("test_ffdiag.jl")
    include("test_plotting.jl")

    include("test_max_iter_warning.jl")
end

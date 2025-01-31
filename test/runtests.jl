using AJD
using Test

using LinearAlgebra
using Diagonalizations: nonDiagonality, LinearFilter
using PosDefManifold: mean
using BenchmarkTools: BenchmarkGroup
using Plots: Plot
using WAV

@testset "AJD.jl" begin
    include("test_jdiag_cardoso.jl")
    include("test_jdiag_eduardpineau.jl")
    include("test_linearfilter.jl")

    # Test the actual result of the algorithm by checking how nondiagonal the matrices are.
    # Tests use nonDiagonality() of Diagonalizations.jl.
    include("test_nondiagonality.jl")
    include("test_nondiagonality_time_correlated.jl")
    include("test_nondiagonality_random_noise.jl")

    include("test_utils.jl")
    include("test_benchmark.jl")
    include("test_ffdiag.jl")
    include("test_plotting.jl")
    include("test_ajd_type_conversion.jl")
    # Make sure diffrent selected algorithms show comparable error handling.
    include("test_max_iter_warning.jl")
    include("test_AJDExtTimedTestdata_functions.jl")
end

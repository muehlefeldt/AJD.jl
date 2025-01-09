# Tests to verfiy benchmark functionality.
using AJD
using BenchmarkTools: BenchmarkGroup

@testset "Benchmark" begin
    # Run the benchmark.
    result = ajd_benchmark(10, 100) 
    # Check returned type.
    @test typeof(result) == BenchmarkGroup
end

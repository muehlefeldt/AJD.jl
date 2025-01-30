# Tests to verfiy benchmark functionality.

@testset "Benchmark" begin
    # Run the benchmark.
    result = diagonalize(:benchmark, 10, 10) 
    # Check returned type.
    @test typeof(result) == BenchmarkGroup
end

# using AJDExtTimedTestdata: generate_correlation_matrix, generate_random_signals
using Random: Xoshiro

#seems like the best way to get functions inside of Extension
ext = Base.get_extension(AJD,:AJDExtTimedTestdata)

@testset "generate_correlation_matrix_erroring" begin
    signal1 = [1 2 3 4]
    signal2 = [1 2 3]
    #check if function throws error if both signals have different sizes
    @test_throws ArgumentError ext.generate_correlation_matrix(signal1,signal2)
    @test_logs (
        :warn,
        "This is a scalar value due to inputs being size nx1. Should be matrix! Increase dimension!"
    ) ext.generate_correlation_matrix(signal2,signal2)
    signal1 = [1 0.5; 1 0.5]
    signal2 = signal1
    @test ext.generate_correlation_matrix(signal1,signal2) ≈ ones(2,2)

end

@testset"generate_random_signals" begin
    test_input = ext.generate_random_measurements(3,3,seed = Xoshiro(123),signal_type = Float64)
    @test typeof(test_input) == Matrix{Float64}
    @test size(test_input) == (3,3)
    @test test_input[1,:] == rand(Xoshiro(123),Float64,3)
end

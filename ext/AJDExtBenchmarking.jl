module AJDExtBenchmarking

using AJD
using BenchmarkTools: @benchmarkable, BenchmarkGroup, tune!, run

function AJD.diagonalize(
    benchmark::Symbol,
    n_dims::Int,
    n_matrices::Int,
)

    if benchmark == :benchmark
        return ajd_benchmark(n_dims, n_matrices)
    else
        throw(ArgumentError("Please use symbol :benchmark to generate plots."))
    end
end

"""
    ajd_benchmark(n_dims::Int, n_matrices::Int)

Run benchmark of implemented algorithms with random inputs.
Prints basic overview of median execution times.
Returns BenchmarkGroup containing detailed results.
"""
function ajd_benchmark(n_dims::Int, n_matrices::Int)
    # Define a parent BenchmarkGroup to contain our suite
    suite = BenchmarkGroup()
    algorithms = [JDiagEdourdPineau(), FFDiag()]

    for alg in algorithms
        name = string(typeof(alg))
        suite[name] = BenchmarkGroup([name])
        suite[name]["exact_diag"] = begin
            @benchmarkable diagonalize(data, algorithm = $alg) setup = (
                data = AJD.get_test_data(
                    :exact_diag,
                    n_dims = $n_dims,
                    n_matrices = $n_matrices,
                )
            )
        end
        suite[name]["approx_diag_large"] = begin
            @benchmarkable diagonalize(data, algorithm = $alg) setup = (
                data = AJD.get_test_data(
                    :approx_diag_large,
                    n_dims = $n_dims,
                    n_matrices = $n_matrices,
                )
            )
        end
        suite[name]["random"] = begin
            @benchmarkable diagonalize(data, algorithm = $alg) setup = (
                data = AJD.get_test_data(
                    :random_noise,
                    n_dims = $n_dims,
                    n_matrices = $n_matrices,
                )
            )
        end
    end

    # Run the actual benchmark.
    tune!(suite)
    results = run(suite, verbose = true)

    # Return BenchmarkGroup for further evaluation.
    return results
end

end # module

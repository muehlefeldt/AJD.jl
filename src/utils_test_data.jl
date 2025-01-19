using WAV

"""
    function get_test_data(
    type::Symbol,
    n_dim::Int = 10,
    n_matrices::Int = 10
    )

Use Symbols `:exact_diag`, `:approx_diag`` and ...
"""
function get_test_data(
    type::Symbol,
    n_dim::Int = 10,
    n_matrices::Int = 10
    )
    if type == :exact_diag
        return random_normal_commuting_matrices(n_dim, n_matrices)
    elseif type == :approx_diag
        data, _ = wavread("test_data/channels3_room69_mix.wav")
        return AJD.generate_testdata(data', delay = 1000, no_of_segments=6)
    elseif type == :approx_diag_large
        data, _ = wavread("test_data/channels3_room69_mix.wav")
        return AJD.generate_testdata(data', delay=100, no_of_segments=10)
    else
        throw(ArgumentError("No valid type of test data selected."))
    end
end
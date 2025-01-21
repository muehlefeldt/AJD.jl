using WAV
"""
    function get_test_data(
    type::Symbol,
    n_dim::Int = 10,
    n_matrices::Int = 10
    )

Usage:
* `:exact_diag`
* `:approx_diag`
* `:random_noice`
"""
#not pretty but if the files should be included the directory is needed to make abs. path
directory = dirname(@__DIR__)
function get_test_data(type::Symbol; n_dim::Int = 10, n_matrices::Int = 10)
    if type == :exact_diag
        return random_normal_commuting_matrices(n_dim, n_matrices)
    elseif type == :approx_diag 
        file = "/test_data/channels3_room69_mix.wav"
        data, _ = wavread(directory * file)
        return AJD.generate_testdata(data', delay = 1000, no_of_segments = 6, 
        show_warning = false)
    elseif type == :approx_diag_large
        file = "/test_data/channels3_room69_mix.wav"
        data, _ = wavread(directory * file)
        return AJD.generate_testdata(data', delay = 1000, no_of_segments = n_matrices+1,
        show_warning = false)
    elseif type ==:random_noice
        A = random_normal_commuting_matrices(n_dim, n_matrices)
        return addrandomnoise(A, same_noise=false)
    #TODO: include this but will need additional arguments in get_test_data
    #look for best way to do this!
    # elseif type == :random_signals
    #     A = AJD.generate_random_signals(n_dim,n_dim)
    #     return AJD.generate_testdata(A,delay = 1000, no_of_segements = n_matrices + 1,
    #     show_warning = false)
    else
        throw(ArgumentError("No valid type of test data selected."))
    end
end
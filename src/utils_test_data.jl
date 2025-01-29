using WAV
"""
    get_test_data(type::Symbol; n_dims::Int = 10, n_matrices::Int = 10)

Get test data to be used as example data or test data.

Usage:
* `get_test_data(:exact_diag)` returns exactly diagonalizable matrices. 
* `get_test_data(:approx_diag)` returns time correlated correlation matrices that are approx. diagonalizable. 
* `get_test_data(:random_noise)` returns complete random matrices.

Parameter `n_dims` and `n_matrices` used to select size of matrices and number of matrices respectivley.
"""
#not pretty but if the files should be included the directory is needed to make abs. path
directory = dirname(@__DIR__)
function get_test_data(type::Symbol; n_dims::Int = 10, n_matrices::Int = 10)
    if type == :exact_diag
        return random_normal_commuting_matrices(n_dims, n_matrices)

    elseif type == :approx_diag
        file = "/test_data/channels3_room69_mix.wav"
        data, _ = wavread(directory * file)
        return AJD.generate_testdata(
            data',
            delay = 1000,
            no_of_segments = 6,
            show_warning = false,
        )
    elseif type == :approx_diag_large
        file = "/test_data/channels3_room69_mix.wav"
        data, _ = wavread(directory * file)
        return AJD.generate_testdata(
            data',
            delay = 1000,
            no_of_segments = n_matrices + 1,
            show_warning = false,
        )
    elseif type == :random_noise
        A = random_normal_commuting_matrices(n_dims, n_matrices)
        return addrandomnoise(A, same_noise = false)
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
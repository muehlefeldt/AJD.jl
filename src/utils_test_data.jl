
"""
    get_test_data(type::Symbol; n_dims::Int = 10, n_matrices::Int = 10)

Get test data to be used as example data or test data.

Usage:
* `get_test_data(:exact_diag)` returns exactly diagonalizable matrices. 
* `get_test_data(:approx_diag)` returns time correlated correlation matrices that are approx. diagonalizable. 
* `get_test_data(:random_noise)` returns complete random matrices.

Parameter `n_dims` and `n_matrices` used to select size of matrices and number of matrices respectively.
"""

function get_test_data(type::Symbol; n_dims::Int = 10, n_matrices::Int = 10)
    if type == :exact_diag
        return random_normal_commuting_matrices(n_dims, n_matrices)

    elseif type in [:random_noise, :random]
        A = random_normal_commuting_matrices(n_dims, n_matrices)
        return addrandomnoise(A, same_noise = false)
    else
        throw(ArgumentError("No valid type of test data selected."))
    end
end
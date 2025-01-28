using AJD

function profile_myfunction(n)
    for i = 1:n
        diagonalize(AJD.get_test_data(:exact_diag, n_dims = 10, n_matrices = 10), algorithm="ffd")# call with desired arguments
    end
end

@profview profile_myfunction(1) # run once to trigger compilation & ignore
@profview profile_myfunction(10) # measure runtime
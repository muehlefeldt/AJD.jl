using AJD

function profile_myfunction(n)
    for i = 1:n
        diagonalize(AJD.get_test_data(:exact_diag, 10, 10), algorithm="jade")# call with desired arguments
    end
end

@profview profile_myfunction(1) # run once to trigger compilation & ignore
@profview profile_myfunction(10) # measure runtime
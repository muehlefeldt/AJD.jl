# Test utility functions.

# Test the generation of random, commuting matrices.
# Used to create test data for the actual diagonalization.
# Real and complex matrices are tested.
@testset "Generate commuting matrices" begin
    matrices = AJD.random_normal_commuting_matrices(4, 10)
    for index in 1:length(matrices)-1
        @test AJD.is_commuting(matrices[index], matrices[index+1])
    end

    matrices = AJD.random_normal_commuting_matrices(4, 10; complex=true)
    for index in 1:length(matrices)-1
        @test AJD.is_commuting(matrices[index], matrices[index+1])
    end
end

# Verify basic function of is_commuting().
@testset "AJD.is_commuting()" begin
    A = [1 1; 1 1]
    C = [0 2; 1 2]
    @test AJD.is_commuting(A,A)
    @test !AJD.is_commuting(A,C)
end

@testset "AJD.get_z_fdiag" begin
    A = reshape(repeat(Matrix(1.0I(3)),outer = (1,2)),3,3,2)
    m,n,_ = size(A)
    for i = 1:m, j = 1:n
        @test AJD.get_z_fdiag(A,i,j) == 2
    end
    A = ones(3,3,2)
    m,n,_ = size(A)
    for i = 1:m, j = 1:n
        @test AJD.get_z_fdiag(A,i,j) == 2
    end

end

@testset "AJD.get_y_fdiag" begin
    D = reshape(repeat(Matrix(1.0I(3)),outer = (1,2)),3,3,2)
    E = (-1.0).*(reshape(repeat(Matrix(1.0I(3)),outer = (1,2)),3,3,2) .-ones(3,3))
    @test AJD.get_y_fdiag(D,E,1,2) == 2
end

# Input verification using AJD.check_input().
# Test of AJD.check_input() directly and of diagonalize().
# diagonalize() calls the verification during normal ops.
# We want to make sure error as correctly handled / passed during all ops situations.
@testset "AJD.check_input()" begin
    # Input of same size matrices, ok max_iter and ok threshold.
    A = AJD.get_test_data(:exact_diag, n_dims=10, n_matrices=10)
    @test_nowarn AJD.check_input(A, 1000, eps())

    # Diffrent size of the matrices.
    # Verify check_input() and diagonalize() itself.
    B = AJD.random_normal_commuting_matrices(9, 10)
    @test_throws ArgumentError AJD.check_input([A..., B...], 1000, eps())
    @test_throws ArgumentError diagonalize([A..., B...])

    # Not acceptable inputs should throw error.
    B = Vector{Matrix}()
    @test_throws MethodError diagonalize(B)
    B = Vector{Matrix}(undef,3)
    @test_throws MethodError diagonalize(B)
    B = Vector{Matrix{Number}}
    @test_throws MethodError diagonalize(B)

    A = AJD.get_test_data(:exact_diag, n_dims=10, n_matrices=10)
    @test_nowarn AJD.check_input(A, 1000, eps())

    # Input of same size matrices, threshold too low.
    A = AJD.get_test_data(:exact_diag, n_dims=10, n_matrices=10)
    @test_throws ArgumentError AJD.check_input(A, 1000, eps()/10)
    @test_throws ArgumentError diagonalize(A, threshold=eps()/10)

    # Input of same size matrices, threshold too high.
    @test_logs (
        :warn,
        "Threshold very high. Recommend threshold of 1e-5 or smaller. Consider machine precision of your system.",
    ) diagonalize(A, threshold=1.0)
    @test_logs (
        :warn,
        "Threshold very high. Recommend threshold of 1e-5 or smaller. Consider machine precision of your system.",
    ) AJD.check_input(A, 1000, 0.2)

    # Input of zero matrices must throw ArgumentError.
    A = AJD.get_test_data(:exact_diag, n_dims=10, n_matrices=10)
    A = [A..., zeros(Float32, 10, 10)]
    @test_throws ArgumentError diagonalize(A)
    @test_throws ArgumentError AJD.check_input(A, 1000, eps())

    # Input of zero matrices must throw ArgumentError.
    A = AJD.get_test_data(:exact_diag, n_dims=10, n_matrices=10)
    A = [zeros(Float32, 10, 10), A...]
    @test_throws ArgumentError diagonalize(A)
    @test_throws ArgumentError AJD.check_input(A, 1000, eps())

    # Input of zero matrices must throw ArgumentError.
    A = [zeros(Float32, 10, 10), zeros(Float32, 10, 10)]
    @test_throws ArgumentError diagonalize(A)
    @test_throws ArgumentError AJD.check_input(A, 1000, eps())

end

@testset "get_off_diag_elements" begin
    input = [[1 2 3; 4 5 6; 7 8 9];;;[0 2 4; 4 2 1; 9 0 1]]
    @test AJD.get_offdiag_elements(input) == [[0 2 3; 4 0 6; 7 8 0];;; [0 2 4; 4 0 1; 9 0 0]]
    input = reshape(repeat(Matrix(1.0I,3,3), outer = (1,2)),3,3,2)
    @test AJD.get_offdiag_elements(input) == zeros(size(input))
end

@testset "frobenius_off_diag_norm" begin
    input = [[1 2 3; 4 5 6; 7 8 9];;;[0 2 4; 4 2 1; 9 0 1]]
    @test AJD.frobenius_offdiag_norm(input) == 296
    input = reshape(repeat(Matrix(1.0I,3,3), outer = (1,2)),3,3,2)
    @test AJD.frobenius_offdiag_norm(input) == 0
end

@testset "get_diag_elements" begin
    input = reshape(repeat(Matrix(1.0I,3,3), outer = (1,2)),3,3,2)
    @test AJD.get_diag_elements(input) == input
end


@testset "random_noise" begin
    input = [Matrix(1.0I,2,2),Matrix(1.0I,2,2)]
    input_noise = AJD.addrandomnoise(input,same_noise = true)
    @test input_noise[1] == input_noise[2]
    #check if noise is the same on both matrices
    #highely unlikely they will ever be the same but test could potentially fail
    input_noise = AJD.addrandomnoise(input,same_noise = false)
    @test input_noise[1] != input_noise[2]
end

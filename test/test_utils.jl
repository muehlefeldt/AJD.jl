# Test utility functions.
# using LinearAlgebra
# using PosDefManifold
using AJD: isstrictly_diagonally_dominant, get_offdiag_elements, get_diag_elements,frobenius_offdiag_norm, sort_offdiag_elements, addrandomnoise,addrandomnoise!

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

# Verify basic function of is_samesize().
@testset "AJD.is_same_size()" begin
    A = ones(2,2)
    B = ones(1,2)
    @test AJD.is_same_size(A,A) == true
    @test AJD.is_same_size(A,B) == false
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
# Input verification.
@testset "AJD.check_input()" begin
    # Input of communting and same size matrices.
    A = AJD.random_normal_commuting_matrices(10, 10)
    @test AJD.check_input(A)

    # Diffrent size of the matrices.
    B = AJD.random_normal_commuting_matrices(9, 10)
    @test !AJD.check_input([A..., B...])
    @test_throws ArgumentError diagonalize([A..., B...])

    # Not acceptable inputs should throw error.
    B = Vector{Matrix}()
    @test_throws MethodError diagonalize(B)
    B = Vector{Matrix}(undef,3)
    @test_throws MethodError diagonalize(B)
    B = Vector{Matrix{Number}}
    @test_throws MethodError diagonalize(B)
end

# Invalid algorithm should lead to a error.
@testset "Invalid Algorithm" begin
    A = AJD.random_normal_commuting_matrices(10, 10)
    @test_throws ArgumentError diagonalize(A, algorithm="hello")
end

@testset "strictly_dominant" begin
    A = 1.0*Matrix(I,3,3)
    @test isstrictly_diagonally_dominant(A) == true
end

@testset "get_off_diag_elements" begin
    input = [[1 2 3; 4 5 6; 7 8 9];;;[0 2 4; 4 2 1; 9 0 1]]
    @test get_offdiag_elements(input) == [[0 2 3; 4 0 6; 7 8 0];;; [0 2 4; 4 0 1; 9 0 0]]
    input = reshape(repeat(Matrix(1.0I,3,3), outer = (1,2)),3,3,2)
    @test get_offdiag_elements(input) == zeros(size(input))
end

@testset "frobenius_off_diag_norm" begin
    input = [[1 2 3; 4 5 6; 7 8 9];;;[0 2 4; 4 2 1; 9 0 1]]
    @test frobenius_offdiag_norm(input) == 296
    input = reshape(repeat(Matrix(1.0I,3,3), outer = (1,2)),3,3,2)
    @test frobenius_offdiag_norm(input) == 0
end

@testset "get_diag_elements" begin
    input = reshape(repeat(Matrix(1.0I,3,3), outer = (1,2)),3,3,2)
    @test get_diag_elements(input) == input
end

@testset "sort_offdiag_elements" begin
    input_1D = [1 2 3; 4 5 6; 7 8 9]
    @test sort_offdiag_elements(input_1D) == float.([2,4,3,7,6,8])
    input = [[1 2 3; 4 5 6; 7 8 9];;;[0 2 4; 4 2 1; 9 0 1]]
    @test sum(sort_offdiag_elements(input)) == sum(get_offdiag_elements(input))
    input = reshape(repeat(Matrix(1.0I,2,2), outer = (1,2)),2,2,2) # has dimension 2x2x2 -> number of offdiag_el = 2*2
    @test sort_offdiag_elements(input) == zeros(4)
end

@testset "random_noise" begin
    input = [Matrix(1.0I,2,2),Matrix(1.0I,2,2)]
    input = addrandomnoise!(input)
    
end

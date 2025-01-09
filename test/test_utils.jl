# Test utility functions.
using AJD

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

    # Not commuting matrices.
    B = AJD.random_normal_commuting_matrices(10, 10)
    @test !AJD.check_input([A..., B...])
    @test_throws ArgumentError diagonalize([A..., B...])

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
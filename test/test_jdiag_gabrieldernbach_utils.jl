using PosDefManifold
using LinearAlgebra
using AJD: isstrictly_diagonally_dominant, get_offdiag_elements, get_diag_elements,frobenius_offdiag_norm, sort_offdiag_elements, addrandomnoise,addrandomnoise!


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

@testset "frobenius_off_diag_normation" begin
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
    @info input
    @info "Jdiag_linear",diagonalize(input)
    @info "FFD_Linear",diagonalize(input,algorithm = "FFD")
    @info "Linear_Filter Eduard", diagonalize(input, algorithm = "jdiag_edourdpineau")
end

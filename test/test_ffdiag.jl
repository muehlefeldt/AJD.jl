using LinearAlgebra

@testset "ffdiag_functionality" begin
    input = AJD.random_normal_commuting_matrices(3,2)
    @info "FFDiag",diagonalize(input, algorithm = "FFD")
end
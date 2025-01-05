using LinearAlgebra
function generate_stacked_psd_matrices(n::Int, count::Int)
    """
    生成 count 个 n x n 的随机正定对称矩阵，并将它们拼接成一个 n x (n * count) 的大矩阵
    generate count-mal n x n positive definite matrices by casting them together
    """
    psd_matrices = [generate_psd_matrix(n) for _ in 1:count]  # 生成矩阵列表
    return hcat(psd_matrices...)  # 按列拼接
end

function generate_psd_matrix(n::Int)
    """
    生成一个 n x n 的随机正定对称矩阵
    generate one nxn positive symmetric matrix
    """
    A = randn(n, n)  # 随机矩阵 a random matrix
    return A * A'    # 保证对称性和正定性 return A*A^T to garentee symmetric and positive definite
end

# Verify an ArgumentError is thrown when algorithm = "jdiag_cardoso" is used with
# complex matrices.
@testset "AJD.jdiag_cardoso() No Complex" begin
    test_input = (1.0) * [Matrix(I, 6, 6), Matrix(I, 6, 6)]
    @test_throws ArgumentError diagonalize(test_input*im, algorithm = "jdiag_cardoso")
end
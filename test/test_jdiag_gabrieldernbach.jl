using PyCall
using AJD
using LinearAlgebra
using Diagonalizations

# Python implementation of the JDiag algorithm.
# The Python code is used to compare implementations.
# PyCall used to execute the Python code.
# Source: https://github.com/gabrieldernbach/approximate_joint_diagonalization/blob/master/jade/jade_cpu.py
py"""
import numpy as np

def jade(A, threshold=10e-16):
    #computes the joint diagonal basis V (m x m) of a set of
    #k square matrices provided in A (k x m x m).

    #It returns the basis V as well as the remaining diagonal
    #and possible residual terms of A.

    #Args:
    #    A: np.ndarray
    #        Tensor of shape (k x m x m) to be diagonalized
    #    threshold: float
    #        stopping criterion, stops is update angle is less than threshold

    #Returns:
    #    A: np.ndarray
    #        Tensor of shape (k x m x m) in approximate diagonal form (approximate eigenvalues)
    #    V: np.ndarray
    #        Matrix of shape (m x m) that contains the approximate joint eigenvectors
    #
    A = np.copy(A)
    m = A.shape[1]
    V = np.eye(m)
    active = 1
    while active == 1:
        active = 0
        for p in range(0, m):
            for q in range(p + 1, m):
                # computation of rotations           
                vecp = A[:, p, p] - A[:, q, q]
                vecm = A[:, p, q] + A[:, q, p]
                ton = vecp @ vecp - vecm @ vecm
                toff = 2 * vecp @ vecm
                theta = 0.5 * np.arctan2(toff, ton + np.sqrt(ton * ton + toff * toff))
                c = np.cos(theta)
                s = np.sin(theta)
                J = np.array([[c, s], [-s, c]])

                active = active | (np.abs(s) > threshold)
                # update of A and V matrices
                if abs(s) > threshold:
                    pair = np.array((p, q))
                    A[:, :, pair] = np.einsum('ij,klj->kli', J, A[:, :, pair])
                    A[:, pair, :] = np.einsum('ij,kjl->kil', J, A[:, pair, :])
                    V[:, pair] = np.einsum('ij,kj->ki', J, V[:, pair])
    return A, V
"""
# End of the Python code.

@testset "JDiag Gabrieldernbach vs Python with I" begin
    # Get basic identity matrix input.
    testinput = (1.0)* [Matrix(I, 6, 6) , Matrix(I, 6, 6)]
    
    # # Get the result of the Julia implementation.
    result = diagonalize(testinput, algorithm = "jdiag")

    # Get the result of Python implementation.
    _, V = py"jade"(testinput)

    @test isapprox(result.F, V)
    @test isapprox(result.F, I)
end

@testset "JDiag Gabrieldernbach vs Python" begin
    for _ in 1:20
        # Get random input matrices.
        test_input = AJD.random_normal_commuting_matrices(10, 6)

        # Get the result of the Julia implementation.
        result = diagonalize(test_input, algorithm="jdiag_gabrieldernbach")

        # Get the result of Python implementation. Only the calculated filter is needed.
        _, V = py"jade"(test_input)

        @test isapprox(V, result.F)
    end
end


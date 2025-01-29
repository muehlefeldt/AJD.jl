# Theoretical Background

This module implements Approximate Joint Diagonalization (AJD) of multiple matrices.

## Problem Statement


The problem of Approximate Joint Diagonalaziation is to find a Matrix $V$ that "maximally diagonalizes" as set of matrices
$\{C_1, C_2, \ldots, C_n\}$. This problem of Joint Diagonalization occurs in various data analysis problems, such as
blind source separation and independent component analysis. The matrices $C_i$ are typically covariance matrices, which
means they are real and symmetric positive semi-definite.

In the special case where all matrices commute pairwise ($C_iC_j = C_jC_i$ for all $i,j$), an exact joint diagonalization exists. However, in practical applications like blind source separation and independent component analysis, this is rarely the case.

## Mathematical Formulation

The problem can be formally stated as an optimization problem:

```math
V = \underset{V \neq 0}{\operatorname{argmin}} \sum_{i=1}^n \text{off}(V^T C_i V)
```

where $\text{off}(M)$ represents the sum of squares of the off-diagonal elements of matrix $M$:

```math
\text{off}(M) = \sum_{i \neq j} |m_{ij}|^2
```
We require the condition that $V\neq 0$ to avoid the trivial solution of $V=0$.

## Implementation

This module provides two distinct algorithms to solve the AJD problem:

1. JADE (Jacobi Angles for Simultaneous Diagonalization)
2. FFDIAG (Fast Frobenius Diagonalization)

### References

Please refere to [References](references.md) for further literature and sources.

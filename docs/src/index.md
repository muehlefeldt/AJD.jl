```@meta
CurrentModule = AJD
```

# AJD.jl Documentation

Documentation for the [AJD.jl](https://github.com/muehlefeldt/AJD.jl) Julia Package to perform Approximate Joint Diagonalization (AJD) on multiple matrices.

## Contents
* [Getting Started Guide](getting-started.md)
* [Detailed documentation](detailed-docs.md)
* [Theory](theoretical-background.md)
* [References](references.md)

## Main AJD.jl functions

Below the main function of the AJD.jl package `diagonalize()` is introduced below. For the entire documentation, including non exported functions, of the package refer to the [detailed docs page](detailed-docs.md).

```@docs
diagonalize(
        A::Vector{<:AbstractMatrix{<:Number}};
        algorithm::String = "jdiag_gabrieldernbach",
        max_iter::Int = 1000,
        threshold::AbstractFloat = eps(),
        plot_matrix::Bool = false,
        plot_convergence::Bool = false
        )
```
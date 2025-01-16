```@meta
CurrentModule = AJD
```

# AJD.jl Documentation

Documentation for the [AJD.jl](https://github.com/gericke-n/AJD.jl) Package to perform Approximate Joint Diagonalization (AJD) on multiple matrices.

Below the main functions of the package are introduced. For the entire documentation of non exported functions of the package  refer to the [detailed docs page](./detailed-docs.md).


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
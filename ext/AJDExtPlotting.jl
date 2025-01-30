module AJDExtPlotting

using AJD
using Statistics: mean
using Plots: Plot, plot, heatmap, theme

function AJD.diagonalize(
    #::Type{T},
    A::Vector{<:AbstractMatrix{<:Number}},
    only_plot::Symbol;
    #x::Plots = None,
    algorithm::String = "jdiag_gabrieldernbach",
    max_iter::Int = 1000,
    threshold::AbstractFloat = eps(),
) #where {T<:Plot}
    #print(x)
    AJD.check_input(A, max_iter, threshold)

    if only_plot == :plot
        F, B, error_array, n_iter = AJD.get_diagonalization(
            A,
            algorithm = algorithm,
            max_iter = max_iter,
            threshold = threshold,
            only_plot = only_plot,
        )

        if n_iter >= max_iter
            @warn "Max iteration was reached. Consider increasing max_iter: diagonalize(M, max_iter=...)."
        end

        p = get_plot(F, B, error_array, algorithm)
    else
        throw(ArgumentError("Please use symbol only_plot=:plot to generate plots."))
    end
    return p
end

"""
    get_plot(
        filter::AbstractArray,
        diag_matrices::AbstractArray, 
        error_array::AbstractArray, 
        name::String)

In case of plot is user selected this generates heatmap plot of the filter matrix and the mean of diagonlaised matrices.
Also a lineplot of the error history of the algorithm calculation is created.
A combined plot is returned.
"""
function get_plot(
    filter::AbstractArray,
    diag_matrices::AbstractArray,
    error_array::AbstractArray,
    name::String,
)

    # Select PLots.jl theme.
    theme(:vibrant)

    # Get the plots and returned one combined plot.
    filter_plot = heatmap(real.(filter), yflip = true, title = "Filter Matrix")
    mean_diag_plot = heatmap(
        real.(mean((diag_matrices), dims = 3)[:, :, 1]),
        yflip = true,
        title = "Mean Diagonalized Matrices",
    )
    error_plot = plot(error_array, w = 3, title = "Error Convergence", label = name)
    return plot(
        filter_plot,
        mean_diag_plot,
        error_plot,
        layout = (3, 1),
        legend = false,
        size = (800, 1200),
    )
end

#export diagonalize

end # module
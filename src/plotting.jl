using Plots

function get_plot(
        filter::AbstractArray,
        diag_matrices::AbstractArray, 
        error_array::AbstractArray, 
        name::String)

    # Select PLots.jl theme.
    theme(:dark)

    # Get the plots and returned one combined plot.
    filter_plot = heatmap(real.(filter), yflip=true, title="Filter Matrix")
    mean_diag_plot = heatmap(real.(mean((diag_matrices), dims=3)[:, :, 1]), yflip=true, title="Mean Diagonalized Matrices")
    error_plot = plot(error_array, w=3, title="Error Convergence", label=name)
    return plot(filter_plot, mean_diag_plot, error_plot, layout=(3, 1), legend=false, size=(400, 1200))
end

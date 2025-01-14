using Plots

"""
    plot_matrix_heatmap(filter::AbstractMatrix, diag_matrices)

Plot a heatmap of the calculated filter matrix and the mean of the diagonlized matrices.
"""
function plot_matrix_heatmap(filter::AbstractMatrix, diag_matrices) 
    theme(:dark)
    # Subplot of the filter matrix.
    filter_plot = heatmap(filter, yflip=true, title="Filter Matrix")

    # Subplot of the mean of all the diagonalised matrices.
    mean_diag_plot = heatmap(mean(diag_matrices, dims=3)[:, :, 1], yflip=true, title="Mean Diagonalized Matrices")

    # Combine subplots and 
    combined_plot = plot(filter_plot, mean_diag_plot, layout=(2, 1), legend=false, size=(400, 800))
    display(combined_plot)
end

"""
    plot_convergence_lineplot(error_array::AbstractArray, name::String)

Plot the convergence error as recorded during the algorithm execution.
"""
function plot_convergence_lineplot(error_array::AbstractArray, name::String)
    theme(:dark)
    line_plot = plot(error_array, w=3, title="Error Convergence", label=name)
    display(line_plot)
end
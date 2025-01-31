"""
    generate_testdata(signal_sources::AbstractArray{<:Function},
    mixing_matrix::AbstractMatrix{<:Number}; <keyword_arguments>)

* `signal_sources`: Array of anonymous functions for generating time series data of the uncorrelated signals ``s_j`` of `BSS` e.g. [ s1 = x-> 1.4*sin(2x), s2 = 2.2sin(x)]
* `mixing_matrix`: mixing matrix by which the signals ``s_j`` are multiplied to get the measurements/observations ``x_i``

First Calculates from a given array of functions, which simulate the uncorrelated signals ``s_j(t)`` and a mixing matrix A, the observations ``x_i``:

``x_i(t) = \\sum_{t = 1}^{T} a_{i,j} s_j(t) ``.

Then a number of time delayed correlation matrices specified by `no_of_corr` is calculated.
# Arguments

* `delay`::Number = 1: time delay between signals
* `sample_time`::Number = 10: length of single time series (same for all observations)
* `no_of_samples`::Int = 100: number of observations made during `sample_time`
* `no_of_cor`::Int = 10: number of observations made over the entire measurement
# Example for `signal_sources` and `mixing_matrix`
```julia
signal_sources = [x->1.6sin(2pi*5x+5)+2sin(2pi*20x+27)+0.5sin(2pi*100x)+1,x->1.2(2pi*11x)+sin(2pi*2x)+0.7sin(2pi*111x+10)]
mixing_matrix = [0.32 -0.43; -1.31 0.34]
```
Mostly deprecated by generate_testdata for discrete values.
"""
#this function is mostly deprecated by the generate_testdata for measurements
#but it is nice for testing against papers which use simulated signals like the one
#in documetation generate_testdata.md
function generate_testdata(signal_sources::AbstractArray{<:Function}, 
    mixing_matrix::AbstractMatrix{<:Number};
    delay::Number = 1, sample_time::Number = 10,
    no_of_samples::Int = 100, no_of_cor::Int = 10)

    rows,columns = size(mixing_matrix)
    if columns != length(signal_sources)
        throw(ArgumentError("Signal source array and mixing matrix have different dimensions (columns of matrix don't match signals in signal_sources)."))
    end

    #initialize the matrix to be diagonalized
    C = Matrix{}[]

    for k in 0:no_of_cor-1

        x = zeros(rows,no_of_samples)
        x_delay = zeros(rows,no_of_samples)
        for row in 1:rows # axes won't work on Array of Functions
            for source in 1:length(signal_sources)
                #needs to be done since broadcasting on vector of anonymous functions doesn't seem to work
                #tried invoke and .|> but couldn't get it to work which is why iteration is necessary

                #observations at starting point t = 1+(k*T)
                x[row,:] = x[row,:] + mixing_matrix[row,source]*signal_sources[source].(range(k*sample_time+1,(k+1)*sample_time,length = no_of_samples))

                #time delayed observations
                x_delay[row,:] = x_delay[row,:] + mixing_matrix[row,source]*signal_sources[source].(range((k*sample_time+1+delay),(k+1)*sample_time+delay, length = no_of_samples))
            end
        end
        push!(C,generate_correlation_matrix(x,x_delay))
    end

    #convert C to appropriate type
    if all(isa.(C, Matrix{Float64})) == true
        C = convert(Vector{Matrix{Float64}},C)
    else
        C = convert(Vector{Matrix{ComplexF64}},C)
    end

    return C
end

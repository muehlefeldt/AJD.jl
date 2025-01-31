module AJDExtTimedTestdata

using AJD
using Random: Xoshiro, rand
using WAV: wavread
using Statistics: cor

directory = dirname(@__DIR__)
function AJD.get_test_data(type::Symbol, abs_path::String; n_dims = 10, n_matrices::Int = 10, delay = 1000,
    no_of_segments = n_matrices + 1)
 
    if type == :approx_diag
        if isempty(abs_path)
            file = "/test_data/channels3_room69_mix.wav"
            abs_path = directory * file
        end 
        data, _ = wavread(abs_path)
        return generate_testdata(
            data',
            delay = delay,
            no_of_segments = 6,
            show_warning = false,
        )
    elseif type == :approx_diag_large
        if isempty(abs_path)
            file = "/test_data/channels3_room69_mix.wav"
            abs_path = directory * file
        end
        data, _ = wavread(abs_path)
        return generate_testdata(
            data',
            delay = delay,
            no_of_segments =no_of_segments,
            show_warning = false,
        )
    
    else
        throw(ArgumentError("No valid type of test data selected."))
    end
end
"""
    generate_correlation_matrix(
        signal_one_data::AbstractArray,
        signal_two_data::AbstractArray,
    )

Calculates correlation matrix between observations ``x_i(t)`` and ``x_i(t+τ)``.

Inputs:
* `signal_one_data`: Array of dimension ``n \\times m``
* `signal_two_data`: Array of dimension ``n \\times m``

Should m be 1 this will give back a scalar value (pearson correlation coeff).
"""
function generate_correlation_matrix(
    signal_one_data::AbstractArray,
    signal_two_data::AbstractArray,
)

    if size(signal_one_data) != size(signal_two_data)
        throw(ArgumentError("Signals have different sizes!"))
    end
    
    C = cor(signal_one_data, signal_two_data, dims = 2)
    if length(C) == 1
        @warn "This is a scalar value due to inputs being size nx1. Should be matrix! Increase dimension!"
    end
    return C
end

"""
    generate_random_meas(
        no_of_meas::Int,
        no_of_samples::Int;
        seed = Xoshiro(),
        signal_type::DataType = Float64,
    )

*`no_of_meas`: used to declare how many measurements ``x_{i}`` are inside of the testdata
*`no_of_samples`: used to declare how many samples each measurement ``x_i`` has
*`seed`: seed for the `rand` function, might need to import the module Random
*`signal_type`: default is Float64 but ComplexF64 is possible as well

Similar to `random_commuting_matrices` function. Gives the opportunity to generate testdata randomly and pluck it into the discrete version of generate_testdata.
"""
function generate_random_measurements(
    no_of_meas::Int,
    no_of_samples::Int;
    seed = Xoshiro(),
    signal_type::DataType = Float64,
)
    measurements = Matrix{signal_type}(undef, no_of_meas, no_of_samples)
    for meas = 1:no_of_meas
        measurements[meas, :] = rand(seed, signal_type, no_of_samples)
    end
    return measurements
end

"""
    generate_testdata(measurements::AbstractArray;
    delay::Number = 10, no_of_segments::Int = 10)
* `measurements`: Matrix of rowwise measurements [``x_1``; ``x_2``;...; ``x_n``]
* `delay`: Time/index shift between observations to be correlated
* `no_of_segments`: Puts `signal_sources` into even segments to be correlated. If the number leads to uneven correlation will throw a warning if show_warning is true
* `show_warning`: If true will show a warning in case segments are uneven. Will lead to one less correlation matrix

Generate Correlation Matrices for discrete observations ``x_i``.

# Known Issue

If your data has a segment with variance close to 0 (e.g. due to all of the values being the same) the correlation matrix will have NaN values inside. Setting the number of segments to a lower value might help.
"""
function generate_testdata(measurements::AbstractArray;
    delay::Number = 10, 
    no_of_segments::Int = 10, 
    show_warning::Bool = true)

    x = copy(measurements)
    rows,columns = size(measurements)

    # if signal has length 100 and delay would be 99 there wouldn't be any data after observation 100 to correlate
    if columns/delay < 2
        throw(ArgumentError("Delay too big. Length of signals divided by delay less than 2. Delay shift would lead to array entry for non existent data."))
    end

    segmentation = columns/no_of_segments

    #in case segmentation leads to uneven segments i.e. 200 points and
    #3 segments -> last observation wouldn't be considered since index will be 66

    if isinteger(segmentation)
        segmentation = Int(segmentation)
    else
        # for testing should be false! could be also done with Logging module and
        # disable_logging(Logging.Warn) but that would supress all warnings for the test case
        # which might not be desirable!
        if show_warning == true
            @warn "Number of Segments leads to segments of different sizes! Will skip last segment.\n"
        end
        no_of_segments -= 1
        segmentation = floor(Int64,segmentation)
    end
    if delay > segmentation
        throw(ArgumentError("Delay bigger than segmentation is currently not supported"))
    end
    # C needs to be declared before the for loop otherwise won't be part of local scope
    # if statement therefore declaring C inside for loop won't work
    # and declaring C as Matrix[] or Array[] will push first entry to be
    # Array[...] or Matrix[...] with all other entries being of type Matrix or Array

    # initialize the matrix set with the type of signal_sources
    # this might lead to an error if the type of signals differs from eachother however
    # if only waveform files are used for testdata generation this shouldn't fail.
    
    C = Matrix{}[]

    for k in 1:no_of_segments-1
        #won't work for offset arrays however since this is only for testdata generation
        #of wav files which are read in shouldn't matter
        x_t = x[:,(k-1)*segmentation+1:k*segmentation]
        x_delay = x[:,(k-1)*segmentation+1+delay:k*segmentation+delay]
        push!(C,generate_correlation_matrix(x_t,x_delay))
    end
    if isnan.(sum(C)) != zeros(rows,rows)
        throw(ArgumentError("Number of segments leads to NaN inside of correlation matrix. See Documentation for further Info."))
    end
    # convert C to appropriate type
    if all(isa.(C, Matrix{Float64})) == true
        C = convert(Vector{Matrix{Float64}},C)
    else
        C = convert(Vector{Matrix{ComplexF64}},C)
    end
    return C
end



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
export generate_testdata,generate_random_measurements, generate_correlation_matrix
end
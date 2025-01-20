# Generating Test Data
The function `AJD.generate_testdata` is currently not exported but will give you the opportunity to create your own testsets to be diagonalized, consisting of an implementation of time delayed correlation matrices. The implementation of the time delayed correlation matrices used are explained [here (p.16)](https://publishup.uni-potsdam.de/opus4-ubp/frontdoor/deliver/index/docId/501/file/ziehe.pdf).

There are currently two ways to generate your own testdata, either using discrete signals or continous signals.

## Continous signals
The function 
```julia 
AJD.generate_testdata(signal_sources::AbstractArray{<:Function}, 
mixing_matrix::AbstractMatrix{<:Number}; delay::Number = 1, sample_time::Number = 10, no_of_samples::Int = 100, no_of_cor::Int = 10)
``` 
for continous signals is a bit complicated to be honest. The input `signal_sources` is an array of anonymous function like:
```julia
# used signals in signal_sources are similar to signals found in https://doi.org/10.21595/jve.2021.21961 p.1709
signal_sources = [t->1.6sin(2pi*5t+5)+2sin(2pi*20t+27)+0.5sin(2pi*100t)+1 , t->1.2(2pi*11t)+sin(2pi*2t)+0.7sin(2pi*111t+10)]
```
and relates to the unmixed signal used in BSS problems. For more information on BSS problems see: [Ziehe (p.5 ff.)](https://publishup.uni-potsdam.de/opus4-ubp/frontdoor/deliver/index/docId/501/file/ziehe.pdf).

Source Signals are time dependent functions but do not necessarily have to be harmonic functions (at least for testing) but can also include ``t^{t}``. 

Be careful though, NaN values can occur if the signal converges to zero over time!

The input `mixing_matrix` can be a random matrix with the only constraint being that it has to have the same column size as signals being recorded e.g. if ``s_j`` with j = 1,2 A has to have two columns.
```julia 
# mixing matrix is the same as in https://doi.org/10.21595/jve.2021.21961 p.1709
mixing_matrix =  [0.32 -0.43; -1.31 0.34]
```
Multiplying the Mixing Matrix and the signal sources leads to the observations ``x_i(t)`` used for the calculated time correlation matrix. (See documentation of function for further explanation of mathematics behind.)

The arguments `sample_time` is used as the ending time of the signal in seconds and  `no_of_samples` is the frequency/sampletime of the signal.

`no_of_cor` is used as an indicator how many time delayed correlation matrices are to be calculated.

The argument `delay` is used to indicate the time shift between two observations ``x(t)`` and ``x(t+\tau)``.

### Example continous signal

```julia 
using AJD
signal_sources = [x->1.6sin(2pi*5x+5)+2sin(2pi*20x+27)+0.5sin(2pi*100x)+1,x->1.2(2pi*11x)+sin(2pi*2x)+0.7sin(2pi*111x+10)]
mixing_matrix = [0.32 -0.43; -1.31 0.34]
testset_data = AJD.generate_testdata(signal_sources,mixing_matrix)
diagonalize(testset_data)
```
## Discrete signals
For realising discrete test data the files `channels2_room69_mix.wav`, `channels3_room69_mix.wav` and  `channels3_room69_mix_shortened.wav` are included in the repository. 

The data was generated using the github code found [here](https://github.com/fakufaku/bss_speech_dataset) - data consists of different soundfiles, which including a differing number of source signals.
The function 
```julia
generate_testdata(signal_sources::AbstractArray; delay::Number = 10, no_of_segments::Int = 10)
``` 
behaves differently for discrete signals than the continous counterpart.

`signal_sources` is a matrix with the observations ``x_i(t)`` ordered rowwise [``x_1(t)``; ``x_2(t)``; ``x_3(t)``].

`delay` is the shift index by which the time delay is emulated. `delay`= 1000 means x[`curr_column` + 1000].

`no_of_segments` divides the `signal_sources` into even segments to be correlated. Future implementation will also include the keyword `points_per_segment` which will behave similar to `sample_time` of analog counterpart, to be more intuitive.

If the segments are not equally divided the method will throw an error. Make sure the length of the observations is dividable by the number of segments or trim your data to be able to!

### Example discrete signal

```julia 
using WAV
data,fs = wavread("channels3_room69_mix.wav")
#currently needed since observations are inside columns
data = data' 
using AJD
testset_data = AJD.generate_testdata(data, delay = 1000,no_of_segments=6)
diagonalize(testset_data)

using WAV
data,fs = wavread("channels2_room69_mix.wav")
#currently needed since observations are inside columns
data = data' 
using AJD
#data has size 2x320336 which is why slicing is used
testset_data = AJD.generate_testdata(data[:,1:320300], delay = 1000,no_of_segments=100)
diagonalize(testset_data)
```

### Generating testset from [https://github.com/fakufaku/bss_speech_dataset](https://github.com/fakufaku/bss_speech_dataset)
If you want to generate your own test data or additional datasets refer to the github repository linked above and clone it. You'll need python!

For the code of the repo to work the following adaptations have to be made:
- the generate_samples.py has to be changed on line 120. Change the deprecated `np.float`to `float`

The generation will take a while and (*probably*) uses around 10GB of storage space. You can always cancel the generation of testdata during the process and get some testsets generated until then. After generating your testset you might want to edit the soundfiles since some of the testsets have 

### Known Issues

The discrete case won't work if one of the segments is a vector of only zeros, since the correlation of this vector will become NaN due to the variance being zero! The function will throw an error in that case!

Make sure to decrease number of segments or if added increase the number of points per segment.
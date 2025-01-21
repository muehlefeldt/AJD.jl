# Generating (time correlated) matrices
The function `AJD.generate_testdata` is currently not exported but will give you the opportunity to create your own testsets to be diagonalized, consisting of an implementation of time delayed correlation matrices. The implementation of the time delayed correlation matrices used are explained [here (p.16)](https://publishup.uni-potsdam.de/opus4-ubp/frontdoor/deliver/index/docId/501/file/ziehe.pdf).

If you just want to start with generating data to be diagonalized currently the recommended way is either to use the code:
```julia
data = AJD.get_test_data(:approx_diag_large, 10, 10)
filter_ = diagonalize(data)
```
which takes a waveform file part of the repository (`channels3_room69_mix.wav`) and reads the audio data from the file, which consists of three different audio signals mixed together and diagonalizes the set of matrices or the code

```julia
data = AJD.get_test_data(:exact_diag, 10, 10)
filter_ = diagonalize(data)
```

with data being a random set of normal commuting matrices with dimension `n_dim` ``\times`` `n_dim` ``\times`` `n_matrices`.

The first way is essentially the same as writing:

```julia
using WAV
data,_ = wavread("test_data/channels3_room69_mix.wav")
#data needs to be transposed since wavread concatenates signals columnwise
data = data' 
using AJD
#no_of_segments is always one more than n_matrices
testset_data = AJD.generate_testdata(data, delay = 1000,no_of_segments=11,show_warning = false)
filter_ = diagonalize(testset_data)
```
however the code above is a bit more versatile since delay can be adjusted and the warning message for wrong segmentation can be enabled.

Another way to generate a testset of correlation matrices from random vectors is provided with:

```julia
using AJD
#you can adjust this!
no_of_signals,no_of_samples = (10,1000)
# also includes a seed and signal_type (either Float64 or ComplexF64)
signals = AJD.generate_random_signals(no_of_signals,no_of_samples)
testset_data = AJD.generate_testdata(signals, delay = 100,no_of_segments=11,show_warning = false)
filter_ = diagonalize(testset_data)
```
which generates random time correlated matrices.

For more advanced testdata generation the following might be interesting.

There are currently two other ways to generate your own testdata, either using discrete or continous signals. The method for discrete testdata is used in the first example code. Continous signals are not the recommended way to implement testdata but is theoretically possible and was inspired by [Zhang,Xin2021](https://doi.org/10.21595/jve.2021.21961 p.1709).

## Continous signals
The function 
```julia 
AJD.generate_testdata(signal_sources::AbstractArray{<:Function}, 
mixing_matrix::AbstractMatrix{<:Number}; delay::Number = 1, sample_time::Number = 10, 
no_of_samples::Int = 100, no_of_cor::Int = 10)
``` 
for continous signals is a bit convoluted to be honest. The input `signal_sources` is an array of anonymous function like:
```julia
# similiar signals to signals found in https://doi.org/10.21595/jve.2021.21961 p.1709
signal_sources = [t->1.6sin(2pi*5t+5)+2sin(2pi*20t+27)+0.5sin(2pi*100t)+1 , t->1.2(2pi*11t)+sin(2pi*2t)+0.7sin(2pi*111t+10)]
```
and relates to the unmixed signal used in BSS problems. For more information on BSS problems see: [Ziehe (p.5 ff.)](https://publishup.uni-potsdam.de/opus4-ubp/frontdoor/deliver/index/docId/501/file/ziehe.pdf).

Source Signals are time dependent functions, but for **testing purposes only** do not necessarily only have to consist of harmonic functions and could include ``t^{t}``. 

Be careful though with choosing functions since if the variance of the signal becomes zero for long times the correlation matrix will include NaN values!

The input `mixing_matrix` can be a random matrix, the only constraint being that it has to have the same column size as signals in the source signals vector. Otherwise the calculation of the measurements/observations ``x_i = \sum_{j=1}^m A_{ij}s_j(t)`` won't be possible. (For further information about calculation of the measurements see [Ziehe (p.6 ff.)](https://publishup.uni-potsdam.de/opus4-ubp/frontdoor/deliver/index/docId/501/file/ziehe.pdf))
```julia 
# mixing matrix is the same as in https://doi.org/10.21595/jve.2021.21961 p.1709
mixing_matrix =  [0.32 -0.43; -1.31 0.34]
```
The measurements are used for calculating the time delayed correlation matrices ``C^k``
The argument `sample_time` is used as the ending time of the signal in seconds and `no_of_samples` is the frequency/sampletime of the signal.

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

The data was generated using the github code found [here](https://github.com/fakufaku/bss_speech_dataset).
The function 
```julia
generate_testdata(signal_sources::AbstractArray; delay::Number = 10, no_of_segments::Int = 10)
``` 
behaves differently for discrete signals than the continous counterpart.

`signal_sources` is a matrix with the observations ``x_i(t)`` ordered rowwise [``x_1(t)``; ``x_2(t)``; ``x_3(t)``].

`delay` is the shift index by which the time delay is emulated. `delay`= 1000 means x[`curr_column` + 1000].

`no_of_segments` divides the `signal_sources` into even segments to be correlated. 

Future implementation could also include the keyword `points_per_segment`, which will behave similar to `sample_time` of the analog counterpart, to be more intuitive. However that could mean (especially if the size of data is not known beforehand) that calculations could take a while.

If the segments are not equally divided the method will throw an error. Make sure the length of the observations is dividable by the number of segments or trim your data to be able to!

### Example discrete signal

```julia 
using WAV
data,fs = wavread("test_data/channels3_room69_mix.wav")
#currently needed since observations are inside columns
data = data' 
using AJD
testset_data = AJD.generate_testdata(data, delay = 1000,no_of_segments=6)
diagonalize(testset_data)

using WAV
data,fs = wavread("test_data/channels2_room69_mix.wav")
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
# Generating (time correlated) matrices
When generating testdata currently several options are available.
Firstly the easiest solution is the function:
```julia
using AJD
data = AJD.get_test_data(:exact_diag, n_dims = 10, n_matrices = 10)
filter_ = diagonalize(data, algorithm = FFDiag())
```
which will create a set of 10 matrices, which are normal and commuting with dimension ``10 \times 10`` and calculate the filter containing ``V``. Instead of `FFDiag()` algorithm `JDiagEdourdPineau()` or all other algorithms can be used (see getting started).

If you want to generate time correlation matrices from a testset the extension `AJD.ExtTimedTestdata` is currently available using:
``` julia
using AJD
using WAV
data = AJD.get_test_data(:approx_diag,"")
filter_ = diagonalize(data, algorithm = FFDiag())
```
which will create a time delayed correlation matrix with 3 measurement signals from the file `channels3_room69_mix.wav`, which can be found in the `test_data` folder and diagonalizes the set of matrices from the correlation matrix. The matrix in `data` has dimension ``5 \\times 10``.

For larger testdata the `:approx_diag_large` symbol can be used.

Implementation of the time delayed correlation matrices used are explained [here (p.16)](https://publishup.uni-potsdam.de/opus4-ubp/frontdoor/deliver/index/docId/501/file/ziehe.pdf).

## Note

The testdata functionality as well is implemented as an extension to AJD.jl. This may require you to add the `WAV` package to your Julia parent environment. This is done to limit the dependency footprint of AJD.jl.


Another way to generate a testset of correlation matrices from random vectors is provided with:

```julia
using AJD
using WAV
#you can adjust this!
no_of_signals,no_of_samples = (10,1000)
#get extension to generate signal
ext = Base.get_extension(AJD,:AJDExtTimedTestdata)
meas = ext.generate_random_measurements(no_of_signals,no_of_samples)
testset_data = ext.generate_testdata(meas, delay = 10,no_of_segments=11,show_warning = false)
filter_ = diagonalize(testset_data)
```
which generates random time correlated matrices. See documentation of `generate_testdata` for further information about the keywords. Currently the delay has to be smaller or equal to the elements in a segment e.g. if `no_of_segments` is 10 and `data` has thousand elements maximum `delay` is 100.

For more advanced testdata generation the following might be interesting.

## Further information on generate_testdata
For realising discrete test data the files `channels2_room69_mix.wav`, `channels3_room69_mix.wav` and  `channels3_room69_mix_shortened.wav` are included in the repository. You can load them into the `get_test_data` function by specifying the absolute path to the file and `:exact_diag_large` as the symbol. Might however not work when working from the temporary env.

### Generating testset from [https://github.com/fakufaku/bss_speech_dataset](`https://github.com/fakufaku/bss_speech_dataset`)
The data was generated using the github code found [here](https://github.com/fakufaku/bss_speech_dataset).

If you want to generate your own test data or additional datasets refer to the github repository linked above and clone it. You'll need python!

For the code of the repo to work the following adaptations have to be made:
- the generate_samples.py has to be changed on line 120. Change the deprecated `np.float`to `float`

The generation will take a while and (*probably*) uses around 10GB of storage space. You can always cancel the generation of testdata during the process and get some testsets generated until then. After generating your testset you might want to edit the soundfiles since some of the testsets have 

### Known Issues

The calculation of time correlation matrices won't work if one of the segments is a vector of only zeros, since the correlation of this vector will become NaN due to the variance being zero! The function will throw an error in that case!

Make sure to decrease number of segments to counteract.


"List of all algorithms. To loop over in tests."
const ALL_ALGORITHMS = [
    JDiagGabrielDernbach(),
    JDiagEdourdPineau(),
    JDiagCardoso(),
    FFDiag()
]

COMPLEX_ALGORITHMS = filter(supportscomplex, ALL_ALGORITHMS)

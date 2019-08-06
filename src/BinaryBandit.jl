module BinaryBandit

# precision constants: The rule of thumb is this: if you aggregate 10^N values, you will lose N digits of precision (http://www.ilikebigbits.com/2017_06_01_float_or_double.html); in our case N=3
# to be used relatively, i.e. if ( a - b ) > BB_numerical_precision_XX * min( abs( a ) + abs( b ) , floatmax( FloatXX ) (see also https://www.floating-point-gui.de/errors/comparison/)
const BB_numerical_precision_64 = 1.0e-13 # instead of 16
const BB_numerical_precision_32 = 1.0e-4 # instead of 7

include( "BB_DP_2.jl" )
include( "BB_evaluation_DP_2.jl" )

end # module

# BinaryBandit
A Julia package for optimization and evaluation of the multi-armed bandit problem with binary (success/failure) responses.

# Installation

Run the following commands in the Julia REPL:

```julia
using Pkg
Pkg.clone("https://github.com/PeterJacko/BinaryBandit.git")
```

# Basic Usage

Currently, only the Bayes-optimal (aka Bayesian decision theoretic) design is implemented. This is the design obtained by solving the problem by dynamic programming (using the backward recursion algorithm).
Currently, only the 2-armed problem is implemented.
Therefore, currently, all the function names start with "DP_2_", indicating that they refer to "dynamic programming" and 2 arms.

For frequentist evaluation of the Bayes-optimal design, the function is
```julia
function DP_2_NS( number_of_allocations :: Int64 , success_probability_arm_1 :: Float64 , success_probability_arm_2 :: Float64 , float_version :: Int64 = Int64( 64 ) , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
```
returning the mean and the variance of the number of successes (NS). For example, use
```julia
horizon = 60
success_probability_arm_1 = 0.3
success_probability_arm_2 = 0.5
BinaryBandit.DP_2_NS( horizon , success_probability_arm_1 , success_probability_arm_2 )
```
which will return `(27.667781619675154, 23.650456467947016)`. These three input arguments are required. The other arguments are optional: float_version is by default set to 64 bits (currently this is set despite changing this value); the parameters of the prior Beta distribution on each arm are by default set to ( 1 , 1 ) meaning that the prior is the uniform distribution.

For Bayesian evaluation or for online optimization of the Bayes-optimal design, the (recommended) function is
```julia
function DP_2_action_lin( number_of_allocations :: Int64 , float_version :: Int64 = Int64( 64 ) , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
```
returning the immediate action and the Bayes-expected number of successes (as Float64). For example, use
```julia
horizon = 60
BinaryBandit.DP_2_action_lin( horizon , 32 )
```
which will return `(3, 38.56234359741211)`. Note that the returned actions can take values 1 (allocation to arm 1 is strictly optimal), 2 (allocation to arm 2 is strictly optimal), or 3 (the difference in expected values between actions 1 and 2 is below a numerical precision threshold, so it can be concluded that allocation to either arm is near-optimal). The first input argument is required. The other arguments are optional: float_version is by default set to 64 bits (but 32 bits memory-wise allows for solving problems with larger horizon at a cost of a negligible inaccuracy); the parameters of the prior Beta distribution on each arm are by default set to ( 1 , 1 ) meaning that the prior is the uniform distribution.
For comparison of accuracy,
```julia
BinaryBandit.DP_2_action_lin( horizon )
```
will return `(3, 38.562343246635564)`, while 
```julia
BinaryBandit.DP_2_action_lin( horizon , 16 )
```
will return `(3, 38.5625)`. Note that although the 16-bit version memory-wise allows for solving problems with larger horizon, the cost is two-fold: the inaccuracy is notable (growing to around 1% for horizon 500) and the runtime is around five times larger.

For offline optimization of the Bayes-optimal design, the (recommended) function is
```julia
function DP_2_policy_bin_lin( number_of_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
```
returning the policy (i.e., actions for all states with linear indices in binary encoding stored in quadruples using hexadecimal numbers) and the Bayes-expected number of successes. For example, use
```julia
horizon = 60
BinaryBandit.DP_2_policy_bin_lin( horizon )
```
which will return
```julia
(UInt8[0x01, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55  …  0x9a, 0x95, 0x56, 0xa9, 0x5a, 0x9a, 0x95, 0xe9, 0xe9, 0x6b], 38.562343246635564)
```
The particular action can be read from the policy using the function
```julia
function DP_2_lin_index( number_of_allocations , number_of_successes_arm_1 , number_of_failures_arm_1 , number_of_successes_arm_2 , number_of_remaining_allocations )
```
which converts a 4D state to linear index. 

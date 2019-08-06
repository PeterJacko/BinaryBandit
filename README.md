# BinaryBandit
A Julia package for optimization and evaluation of the multi-armed bandit problem with binary (success/failure) responses. It is being developed by the G.O.A.L. (Group on Optimal Adaptive Learning), see https://www.lancaster.ac.uk/staff/jacko/goal/ for papers presenting results obtained with this package.

# Installation

Run the following commands in the Julia REPL:

```julia
using Pkg
Pkg.clone("https://github.com/PeterJacko/BinaryBandit.git")
```

# Basic Usage

Currently, only the Bayes-optimal (aka Bayesian decision theoretic) design is implemented. This is the design obtained by solving the problem by dynamic programming (using the backward recursion algorithm).
Currently, only the 2-armed problem is implemented.
Therefore, currently, all the function names start with `DP_2_`, indicating that they refer to *dynamic programming* and *2* arms.

For **frequentist evaluation** of the Bayes-optimal design, the function is
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

For **Bayesian evaluation** or for **online optimization** of the Bayes-optimal design, the (recommended) function is
```julia
function DP_2_action_lin( number_of_allocations :: Int64 , float_version :: Int64 = Int64( 64 ) , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
```
returning the immediate action and the Bayes-expected number of successes (as Float64). For example, use
```julia
horizon = 60
BinaryBandit.DP_2_action_lin( horizon , 32 )
```
which will return `(3, 38.56234359741211)`. Note that the returned actions can take values 1 (allocation to arm 1 is strictly optimal), 2 (allocation to arm 2 is strictly optimal), or 3 (the difference in expected values between actions 1 and 2 is below a numerical precision threshold, so it can be concluded that allocation to either arm is near-optimal). The first input argument is required. The other arguments are optional: float_version is by default set to 64 bits (but 32 bits memory-wise allows for solving problems with larger horizon at a cost of a negligible inaccuracy); the parameters of the prior Beta distribution on each arm are by default set to ( 1 , 1 ) meaning that the prior is the uniform distribution.
For comparison of accuracy, the 64-bit version
```julia
BinaryBandit.DP_2_action_lin( horizon )
```
will return `(3, 38.562343246635564)`, while the 16-bit version
```julia
BinaryBandit.DP_2_action_lin( horizon , 16 )
```
will return `(3, 38.5625)`. Note that although the 16-bit version memory-wise allows for solving problems with larger horizon, the cost is two-fold: the inaccuracy is notable (growing to around 1% for horizon 500) and the runtime is around five times larger.

For **offline optimization** of the Bayes-optimal design, the (recommended) function is
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
A particular action can be read from the policy vector using the function
```julia
function DP_2_lin_index( number_of_allocations , number_of_successes_arm_1 , number_of_failures_arm_1 , number_of_successes_arm_2 , number_of_remaining_allocations )
```
which converts a 4D state to linear index.

Note that the runtime grows quickly with the horizon, taking a few minutes for horizon around 1000, a few hours for horizon around 2000 and a few days for horizon around 4000. In terms of memory requirements, 32 GB RAM is able to store the whole policy for approximately horizon of 1500.

# Advanced Usage

See `test\runtests.jl` for a full list of functions, which include different implementations of those described above. In general, those starting with `DP_2_action` return just the immediate action, while those starting with `DP_2_policy` return the whole policy. Suffixes `_lin` and `_bin` indicate linear and binary encoding, respectively, of the policy and/or the value function. Suffix `_with_administration` is present in functions which assume that at the end of the trial, there will be a given additional number of allocations made in the subjects population using the arm that looks better at the end of the trial. A more general version of this are functions with suffix `_with_finale` (not included in `runtests.jl`), in which the final-period value function can be defined.

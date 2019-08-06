# Design: DP

# Float64 version of value_to_go
function DP_2_finale_mean( number_of_allocations :: Int64 , success_probability_arm_1 :: Float64 , success_probability_arm_2 :: Float64 , value_to_go_evaluation :: Array{ Float64 , 1 } , value_to_go :: Array{ Float64 , 1 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Uses linear indexing for value_to_go
  # Output is the expected value of the finale array (the sum of all finale elements multiplied by their probabilities under success_probability_arm_1 & 2 )

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0
            value_to_go_lin_index = 0
            for number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
                @inbounds begin

                value_to_go_lin_index += 1
                belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

                value_to_go_if_action_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + number_of_failures_arm_1 + 1 ] ) + ( 1.0 - belief_of_success_arm_1 ) * value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) ]
                value_to_go_if_action_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) , 2 ) ] ) + ( 1.0 - belief_of_success_arm_2 ) * value_to_go[ value_to_go_lin_index + number_of_failures_arm_1 + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) ]

                if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[ value_to_go_lin_index ] = value_to_go_if_action_1

                    value_to_go_evaluation[ value_to_go_lin_index ] = success_probability_arm_1 * value_to_go_evaluation[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + number_of_failures_arm_1 + 1 ] + ( 1.0 - success_probability_arm_1 ) * value_to_go_evaluation[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) ]

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[ value_to_go_lin_index ] = value_to_go_if_action_2

                    value_to_go_evaluation[ value_to_go_lin_index ] = success_probability_arm_2 * ( value_to_go_evaluation[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) , 2 ) ] ) + ( 1.0 - success_probability_arm_2 ) * value_to_go_evaluation[ value_to_go_lin_index + number_of_failures_arm_1 + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) ]

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    value_to_go[ value_to_go_lin_index ] = ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2

                    value_to_go_evaluation[ value_to_go_lin_index ] = ( success_probability_arm_1 * value_to_go_evaluation[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + number_of_failures_arm_1 + 1 ] + ( 1.0 - success_probability_arm_1 ) * value_to_go_evaluation[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) ] + success_probability_arm_2 * ( value_to_go_evaluation[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) , 2 ) ] ) + ( 1.0 - success_probability_arm_2 ) * value_to_go_evaluation[ value_to_go_lin_index + number_of_failures_arm_1 + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) ] ) / 2

                end

            end # @inbounds
            end
        end

        return value_to_go_evaluation[ 1 ]

end

function DP_2_NS( number_of_allocations :: Int64 , success_probability_arm_1 :: Float64 , success_probability_arm_2 :: Float64 , float_version :: Int64 = Int64( 64 ) , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Uses linear indexing for value_to_go
  # Output is the mean and the variance of the number of successes CHECK (as Float64 irrespectively of the value of "precision")

#        if float_version == 16
#
#            value_to_go_16 :: Array{ Float16 , 1 } = zeros( Float16 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
#            DP_2_finale_mean( number_of_allocations , value_to_go_16 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )
#
#        elseif float_version == 32
#
#            value_to_go_32 :: Array{ Float32 , 1 } = zeros( Float32 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
#            DP_2_finale_mean( number_of_allocations , value_to_go_32 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )
#
#        else

            value_to_go :: Array{ Float64 , 1 } = zeros( Float64 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )

            value_to_go_evaluation :: Array{ Float64 , 1 } = zeros( Float64 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
            value_to_go_lin_index = 0
            for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                value_to_go_lin_index += 1
                @inbounds begin
                    value_to_go_evaluation[ value_to_go_lin_index ] = number_of_successes_arm_1 + number_of_successes_arm_2
                end # @inbounds
            end

            NS_mean = DP_2_finale_mean( number_of_allocations , success_probability_arm_1 , success_probability_arm_2 , value_to_go_evaluation , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

            value_to_go = zeros( Float64 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )

            value_to_go_lin_index = 0
            for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                value_to_go_lin_index += 1
                @inbounds begin
                    value_to_go_evaluation[ value_to_go_lin_index ] = ( number_of_successes_arm_1 + number_of_successes_arm_2 - NS_mean ) ^ 2
                end # @inbounds
            end

            NS_variance = DP_2_finale_mean( number_of_allocations , success_probability_arm_1 , success_probability_arm_2 , value_to_go_evaluation , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

            return NS_mean , NS_variance

#        end

end

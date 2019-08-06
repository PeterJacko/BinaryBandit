function DP_2_lin_index( number_of_allocations , number_of_successes_arm_1 , number_of_failures_arm_1 , number_of_successes_arm_2 , number_of_remaining_allocations )
# This converts a 4D state to linear index
# number_of_successes_arm_1 , number_of_failures_arm_1 , number_of_successes_arm_2 \ge 0
# number_of_allocations , number_of_remaining_allocations \ge 1
# number_of_successes_arm_1 + number_of_failures_arm_1 + number_of_successes_arm_2 + number_of_remaining_allocations \le number_of_allocations

    return div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - ( number_of_allocations - number_of_remaining_allocations + 1 ) * ( number_of_allocations - number_of_remaining_allocations + 2 ) * ( number_of_allocations - number_of_remaining_allocations + 3 ) * ( number_of_allocations - number_of_remaining_allocations + 4 ) , 24 ) + div( ( number_of_allocations - number_of_remaining_allocations + 1 ) * ( number_of_allocations - number_of_remaining_allocations + 2 ) * ( number_of_allocations - number_of_remaining_allocations + 3 ) - ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 1 ) * ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 2 ) * ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 3 ) , 6 ) + div( ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 1 ) * ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 + 2 ) - ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 - number_of_failures_arm_1 + 1 ) * ( number_of_allocations - number_of_remaining_allocations - number_of_successes_arm_2 - number_of_failures_arm_1 + 2 ) , 2 ) + number_of_successes_arm_1 + 1

end

###############################################################################################################################

# Float64 version of value_to_go
function DP_2_action_lin_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float64 , 1 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Uses linear indexing for value_to_go
  # Output is the immediate action and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        value_to_go_if_action_1 = 0.0 # needed after the for loop
        value_to_go_if_action_2 = 0.0 # needed after the for loop
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

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[ value_to_go_lin_index ] = value_to_go_if_action_2

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    value_to_go[ value_to_go_lin_index ] = ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2

                end

            end # @inbounds
            end
        end
        if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 1 ) # action 1 (i.e., arm 1)

        elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 2 ) # action 2 (i.e., arm 2)

        else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

            action = Int8( 3 ) # randomise between actions 1 and 2

        end

        return action , value_to_go[ 1 ]

end

# Float32 version of value_to_go
function DP_2_action_lin_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float32 , 1 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Uses linear indexing for value_to_go
  # Output is the immediate action and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        value_to_go_if_action_1 = 0.0 # needed after the for loop
        value_to_go_if_action_2 = 0.0 # needed after the for loop
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

                    value_to_go[ value_to_go_lin_index ] = Float32( value_to_go_if_action_1 )

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[ value_to_go_lin_index ] = Float32( value_to_go_if_action_2 )

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    value_to_go[ value_to_go_lin_index ] = Float32( ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2 )

                end

            end # @inbounds
            end
        end
        if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 1 ) # action 1 (i.e., arm 1)

        elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 2 ) # action 2 (i.e., arm 2)

        else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

            action = Int8( 3 ) # randomise between actions 1 and 2

        end

        return action , Float64( value_to_go[ 1 ] )

end

# Float16 version of value_to_go
function DP_2_action_lin_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float16 , 1 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Uses linear indexing for value_to_go
  # Output is the immediate action and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        value_to_go_if_action_1 = 0.0 # needed after the for loop
        value_to_go_if_action_2 = 0.0 # needed after the for loop
        for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0
            value_to_go_lin_index = 0
            for number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
            @inbounds begin

                value_to_go_lin_index += 1
                belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

                # THESE TWO LINES ARE VERY SLOW
                value_to_go_if_action_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + number_of_failures_arm_1 + 1 ] ) + ( 1.0 - belief_of_success_arm_1 ) * value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) + ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) ]
                value_to_go_if_action_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[ value_to_go_lin_index + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) , 2 ) ] ) + ( 1.0 - belief_of_success_arm_2 ) * value_to_go[ value_to_go_lin_index + number_of_failures_arm_1 + div( ( number_of_observed_responses + 2 ) * ( number_of_observed_responses + 3 ) - ( number_of_observed_responses - number_of_successes_arm_2 + 2 ) * ( number_of_observed_responses - number_of_successes_arm_2 + 3 ) , 2 ) ]

                if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    # THE RIGHT-HAND SIDE IS VERY SLOW
                    value_to_go[ value_to_go_lin_index ] = Float16( value_to_go_if_action_1 )

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    # THE RIGHT-HAND SIDE IS VERY SLOW
                    value_to_go[ value_to_go_lin_index ] = Float16( value_to_go_if_action_2 )

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    # THE RIGHT-HAND SIDE IS VERY SLOW
                    value_to_go[ value_to_go_lin_index ] = Float16( ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2 )

                end

            end # @inbounds
            end
        end
        if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 1 ) # action 1 (i.e., arm 1)

        elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 2 ) # action 2 (i.e., arm 2)

        else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

            action = Int8( 3 ) # randomise between actions 1 and 2

        end

        return action , Float64( value_to_go[ 1 ] )

end

function DP_2_action_lin( number_of_allocations :: Int64 , float_version :: Int64 = Int64( 64 ) , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) ) #
  # This function implements the DP design for 2 arms
  # Uses linear indexing for value_to_go
  # Output is the immediate action and the Bayes-expected number of successes (as Float64 irrespectively of the value of "precision")

        if float_version == 16

            value_to_go_16 :: Array{ Float16 , 1 } = zeros( Float16 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
            DP_2_action_lin_with_finale( number_of_allocations , value_to_go_16 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        elseif float_version == 32

            value_to_go_32 :: Array{ Float32 , 1 } = zeros( Float32 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
            DP_2_action_lin_with_finale( number_of_allocations , value_to_go_32 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        else

            value_to_go :: Array{ Float64 , 1 } = zeros( Float64 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
            DP_2_action_lin_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        end

end

function DP_2_action_lin_with_administration( number_of_allocations :: Int64 , number_of_administered_allocations :: Int64 , float_version :: Int64 = Int64( 64 ) , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms; the better arm (comparing the beliefs) is administered after the trial number_of_administered_allocations times
  # Uses linear indexing for value_to_go
  # Output is the immediate action and the Bayes-expected number of successes (as Float64 irrespectively of the value of "precision")

        if number_of_administered_allocations == 0

            DP_2_action_lin( number_of_allocations , float_version , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        else

            if float_version == 16

                value_to_go_16 :: Array{ Float16 , 1 } = zeros( Float16 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
                value_to_go_lin_index = 0
                for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                    belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                    belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                    value_to_go_lin_index += 1
                    @inbounds begin
                        value_to_go_16[ value_to_go_lin_index ] = Float16( max( belief_of_success_arm_1 , belief_of_success_arm_2 ) )
                    end # @inbounds
                end
                value_to_go_16 .*= number_of_administered_allocations
                DP_2_action_lin_with_finale( number_of_allocations , value_to_go_16 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

            elseif float_version == 32

                value_to_go_32 :: Array{ Float32 , 1 } = zeros( Float32 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
                value_to_go_lin_index = 0
                for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                    belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                    belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                    value_to_go_lin_index += 1
                    @inbounds begin
                        value_to_go_32[ value_to_go_lin_index ] = Float32( max( belief_of_success_arm_1 , belief_of_success_arm_2 ) )
                    end # @inbounds
                end
                value_to_go_32 .*= number_of_administered_allocations
                DP_2_action_lin_with_finale( number_of_allocations , value_to_go_32 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

            else #all else is redirected to float_version == 64

                value_to_go :: Array{ Float64 , 1 } = zeros( Float64 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
                value_to_go_lin_index = 0
                for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                    belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                    belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                    value_to_go_lin_index += 1
                    @inbounds begin
                        value_to_go[ value_to_go_lin_index ] = max( belief_of_success_arm_1 , belief_of_success_arm_2 )
                    end # @inbounds
                end
                value_to_go .*= number_of_administered_allocations
                DP_2_action_lin_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

            end

        end

end

###############################################################################################################################

# Float64 version of value_to_go
function DP_2_action_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float64 , 3 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] of size ( number_of_allocations )^3 is the value_to_go at the finale (after all allocations are made)
  # Output is the immediate action and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        value_to_go_if_action_1 = 0.0 # needed after the for loop
        value_to_go_if_action_2 = 0.0 # needed after the for loop
        for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0 , number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
            @inbounds begin

                belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

                value_to_go_if_action_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 + 1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] ) + ( 1.0 - belief_of_success_arm_1 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 + 1 ,1+ number_of_successes_arm_2 ]
                value_to_go_if_action_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 + 1 ] ) + ( 1.0 - belief_of_success_arm_2 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ]

                if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_1

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_2

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2

                end

            end # @inbounds
        end
        if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 1 ) # action 1 (i.e., arm 1)

        elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 2 ) # action 2 (i.e., arm 2)

        else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

            action = Int8( 3 ) # randomise between actions 1 and 2

        end

        return action , value_to_go[1 ,1 ,1 ]

end

# Float32 version of value_to_go (calculations are done using Float64 variables, only the resulting value converted to Float32 for memory saving)
function DP_2_action_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float32 , 3 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] of size ( number_of_allocations )^3 is the value_to_go at the finale (after all allocations are made)
  # Output is the immediate action and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        value_to_go_if_action_1 = 0.0 # needed after the for loop
        value_to_go_if_action_2 = 0.0 # needed after the for loop
        for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0 , number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
            @inbounds begin

                belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

                value_to_go_if_action_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 + 1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] ) + ( 1.0 - belief_of_success_arm_1 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 + 1 ,1+ number_of_successes_arm_2 ]
                value_to_go_if_action_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 + 1 ] ) + ( 1.0 - belief_of_success_arm_2 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ]

                if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    # THE RIGHT-HAND SIDE IS VERY SLOW
                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = Float32( value_to_go_if_action_1 )

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    # THE RIGHT-HAND SIDE IS VERY SLOW
                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = Float32( value_to_go_if_action_2 )

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    # THE RIGHT-HAND SIDE IS VERY SLOW
                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = Float32( ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2 )

                end

            end # @inbounds
        end
        if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 1 ) # action 1 (i.e., arm 1)

        elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 2 ) # action 2 (i.e., arm 2)

        else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

            action = Int8( 3 ) # randomise between actions 1 and 2

        end

        return action , Float64( value_to_go[1 ,1 ,1 ] )

end

# Float16 version of value_to_go (calculations are done using Float64 variables, only the resulting value converted to Float16 for memory saving)
function DP_2_action_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float16 , 3 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] of size ( number_of_allocations )^3 is the value_to_go at the finale (after all allocations are made)
  # Output is the immediate action and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        value_to_go_if_action_1 = 0.0 # needed after the for loop
        value_to_go_if_action_2 = 0.0 # needed after the for loop
        for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0 , number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
            @inbounds begin

                belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

                # THESE TWO LINES ARE VERY SLOW
                value_to_go_if_action_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 + 1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] ) + ( 1.0 - belief_of_success_arm_1 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 + 1 ,1+ number_of_successes_arm_2 ]
                value_to_go_if_action_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 + 1 ] ) + ( 1.0 - belief_of_success_arm_2 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ]

                if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    # THE RIGHT-HAND SIDE IS VERY SLOW
                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = Float16( value_to_go_if_action_1 )

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    # THE RIGHT-HAND SIDE IS VERY SLOW
                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = Float16( value_to_go_if_action_2 )

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    # THE RIGHT-HAND SIDE IS VERY SLOW
                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = Float16( ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2 )

                end

            end # @inbounds
        end
        if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 1 ) # action 1 (i.e., arm 1)

        elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

            action = Int8( 2 ) # action 2 (i.e., arm 2)

        else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

            action = Int8( 3 ) # randomise between actions 1 and 2

        end

        return action , Float64( value_to_go[1 ,1 ,1 ] )

end

function DP_2_action( number_of_allocations :: Int64 , float_version :: Int64 = Int64( 64 ) , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) ) #
  # This function implements the DP design for 2 arms
  # Output is the immediate action and the Bayes-expected number of successes (as Float64 irrespectively of the value of "precision")

        if float_version == 16

            value_to_go_16 :: Array{ Float16 , 3 } = zeros( Float16 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
            DP_2_action_with_finale( number_of_allocations , value_to_go_16 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        elseif float_version == 32

            value_to_go_32 :: Array{ Float32 , 3 } = zeros( Float32 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
            DP_2_action_with_finale( number_of_allocations , value_to_go_32 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        else

            value_to_go :: Array{ Float64 , 3 } = zeros( Float64 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
            DP_2_action_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        end

end

function DP_2_action_with_administration( number_of_allocations :: Int64 , number_of_administered_allocations :: Int64 , float_version :: Int64 = Int64( 64 ) , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms; the better arm (comparing the beliefs) is administered after the trial number_of_administered_allocations times
  # Output is the immediate action and the Bayes-expected number of successes (as Float64 irrespectively of the value of "precision")

        if number_of_administered_allocations == 0

            DP_2_action( number_of_allocations , float_version , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        else

            if float_version == 16

                value_to_go_16 :: Array{ Float16 , 3 } = zeros( Float16 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
                for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                    belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                    belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                    @inbounds begin
                        value_to_go_16[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = Float16( max( belief_of_success_arm_1 , belief_of_success_arm_2 ) )
                    end # @inbounds
                end
                value_to_go_16 .*= number_of_administered_allocations
                DP_2_action_with_finale( number_of_allocations , value_to_go_16 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

            elseif float_version == 32

                value_to_go_32 :: Array{ Float32 , 3 } = zeros( Float32 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
                for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                    belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                    belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                    @inbounds begin
                        value_to_go_32[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = Float32( max( belief_of_success_arm_1 , belief_of_success_arm_2 ) )
                    end # @inbounds
                end
                value_to_go_32 .*= number_of_administered_allocations
                DP_2_action_with_finale( number_of_allocations , value_to_go_32 , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

            else #all else is redirected to float_version == 64

                value_to_go :: Array{ Float64 , 3 } = zeros( Float64 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
                for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                    belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                    belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                    @inbounds begin
                        value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = max( belief_of_success_arm_1 , belief_of_success_arm_2 )
                    end # @inbounds
                end
                value_to_go .*= number_of_administered_allocations
                DP_2_action_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

            end

        end

end

###############################################################################################################################

function DP_2_policy_bin_lin_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float64 , 1 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Uses linear indexing for value_to_go
  # Output is the policy (i.e., actions for all states) with linear indices in binary encoding and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        action_bin = zeros( UInt8 , div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 * 4 , 96 ) + 1 ) # this is the correct size of the action in binary encoding
        mod_lin = mod( div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 * 4 , 24 ) + 1 , 4 ) # modulo of size of action_lin divided by 4
        action_bin_temp :: UInt8 = 0x00 # encoded up to 4 actions to be written at the next position of action_bin
        action_bin_temp_index = mod( 4 - mod_lin , 4 ) # number of actions currently encoded in action_bin_temp; initialise (up to 3 first values may be unused) so that the final value in action_lib refers to the initial state
        lin_index = 0 # last filled position of vector action_bin
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
                    action_bin_temp = ( action_bin_temp * 0x04 ) + 0x01 # action 1 is optimal, action 2 is not optimal

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[ value_to_go_lin_index ] = value_to_go_if_action_2
                    action_bin_temp = ( action_bin_temp * 0x04 ) + 0x02 # action 1 is not optimal, action 2 is optimal

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    value_to_go[ value_to_go_lin_index ] = ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2
                    action_bin_temp = ( action_bin_temp * 0x04 ) + 0x03 # randomise between actions 1 and 2 as both actions are optimal

                end

                action_bin_temp_index += 1
                if action_bin_temp_index == 4
                    lin_index += 1
                    action_bin[ lin_index ] = action_bin_temp
                    action_bin_temp = 0x00
                    action_bin_temp_index = 0
                end

            end # @inbounds
            end
        end

        return action_bin , value_to_go[ 1 ]

end

function DP_2_policy_bin_lin( number_of_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) ) #
  # This function implements the DP design for 2 arms
  # Output is the policy (i.e., actions for all states) with linear indices in binary encoding and the Bayes-expected number of successes

        value_to_go :: Array{ Float64 , 1 } = zeros( Float64 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
        DP_2_policy_bin_lin_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

end

function DP_2_policy_bin_lin_with_administration( number_of_allocations :: Int64 , number_of_administered_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms; the better arm (comparing the beliefs) is administered after the trial number_of_administered_allocations times
  # Uses linear indexing for value_to_go
  # Output is the policy (i.e., actions for all states) with linear indices in binary encoding and the Bayes-expected number of successes

        if number_of_administered_allocations == 0

            DP_2_policy_bin_lin( number_of_allocations , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        else

                value_to_go :: Array{ Float64 , 1 } = zeros( Float64 , div( ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 , 6 ) + 1 )
                value_to_go_lin_index = 0
                for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                    belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                    belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                    value_to_go_lin_index += 1
                    @inbounds begin
                        value_to_go[ value_to_go_lin_index ] = max( belief_of_success_arm_1 , belief_of_success_arm_2 )
                    end # @inbounds
                end
                value_to_go .*= number_of_administered_allocations
                DP_2_policy_bin_lin_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        end

end

###############################################################################################################################

function DP_2_policy_bin_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float64 , 3 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] of size ( number_of_allocations )^3 is the value_to_go at the finale (after all allocations are made)
  # Output is the policy (i.e., actions for all states) with linear indices in binary encoding and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        action_bin = zeros( UInt8 , div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 * 4 , 96 ) + 1 ) # this is the correct size of the action in binary encoding
        mod_lin = mod( div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 * 4 , 24 ) + 1 , 4 ) # modulo of size of action_lin divided by 4
        action_bin_temp :: UInt8 = 0x00 # encoded up to 4 actions to be written at the next position of action_bin
        action_bin_temp_index = mod( 4 - mod_lin , 4 ) # number of actions currently encoded in action_bin_temp; initialise (up to 3 first values may be unused) so that the final value in action_lib refers to the initial state
        lin_index = 0 # last filled position of vector action_bin
        for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0 , number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
            @inbounds begin

                belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

                value_to_go_if_action_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 + 1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] ) + ( 1.0 - belief_of_success_arm_1 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 + 1 ,1+ number_of_successes_arm_2 ]
                value_to_go_if_action_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 + 1 ] ) + ( 1.0 - belief_of_success_arm_2 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ]

                if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_1
                    action_bin_temp = ( action_bin_temp * 0x04 ) + 0x01 # action 1 is optimal, action 2 is not optimal

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_2
                    action_bin_temp = ( action_bin_temp * 0x04 ) + 0x02 # action 1 is not optimal, action 2 is optimal

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2
                    action_bin_temp = ( action_bin_temp * 0x04 ) + 0x03 # randomise between actions 1 and 2 as both actions are optimal

                end

                action_bin_temp_index += 1
                if action_bin_temp_index == 4
                    lin_index += 1
                    action_bin[ lin_index ] = action_bin_temp
                    action_bin_temp = 0x00
                    action_bin_temp_index = 0
                end

            end # @inbounds
        end

        return action_bin , value_to_go[1 ,1 ,1 ]

end

function DP_2_policy_bin( number_of_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Output is the policy (i.e., actions for all states) with linear indices and the Bayes-expected number of successes

        value_to_go :: Array{ Float64 , 3 } = zeros( Float64 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
        DP_2_policy_bin_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

end

function DP_2_policy_bin_with_administration( number_of_allocations :: Int64 , number_of_administered_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms; the better arm (comparing the beliefs) is administered after the trial number_of_administered_allocations times
  # Output is the policy (i.e., actions for all states) with linear indices and the Bayes-expected number of successes

        if number_of_administered_allocations == 0

            DP_2_policy_bin( number_of_allocations , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        else

            # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations )
            value_to_go :: Array{ Float64 , 3 } = zeros( Float64 , number_of_allocations + 1, number_of_allocations + 1, number_of_allocations + 1 )
            # value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] of size ( number_of_allocations )^3 is the value_to_go at the finale (after all allocations are made)
            for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                @inbounds begin
                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = max( belief_of_success_arm_1 , belief_of_success_arm_2 )
                end # @inbounds
            end
            value_to_go .*= number_of_administered_allocations
            DP_2_policy_bin_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        end

end

###############################################################################################################################

function DP_2_policy_lin_with_finale_new( number_of_allocations :: Int64 , value_to_go :: Array{ Float64 , 3 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] of size ( number_of_allocations )^3 is the value_to_go at the finale (after all allocations are made)
  # Output is the policy (i.e., actions for all states) with linear indices and the Bayes-expected number of successes

    @inbounds begin

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        action_lin = zeros( Int8 , div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 * 4 , 24 ) + 1 ) # this is the correct size of the action
        lin_index = 0
        for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0 , number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )

                lin_index += 1
                belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) * ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 )
                belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) * ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                # both believes are multiplied by the product of denominators, i.e., ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) * ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
#                belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
#                belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

                value_to_go_if_action_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 + 1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] ) + ( ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) * ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) - belief_of_success_arm_1 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 + 1 ,1+ number_of_successes_arm_2 ]
                value_to_go_if_action_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 + 1 ] ) + ( ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) * ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) - belief_of_success_arm_2 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ]
                # both values are multiplied by the product of denominators, i.e., ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) * ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )

                if value_to_go_if_action_1 > ( BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 ) + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_1 / ( ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) * ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                    action_lin[ lin_index ] = 1 # action 1 (i.e., arm 1)

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_2 / ( ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) * ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                    action_lin[ lin_index ] = 2 # action 2 (i.e., arm 2)

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = ( value_to_go_if_action_1 + value_to_go_if_action_2 ) * 0.5 / ( ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) * ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                    action_lin[ lin_index ] = 3 # randomise between actions 1 and 2

                end

        end

    end # @inbounds

    return action_lin , value_to_go[1 ,1 ,1 ]

end

function DP_2_policy_lin_new( number_of_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Output is the policy (i.e., actions for all states) with linear indices and the Bayes-expected number of successes

        value_to_go :: Array{ Float64 , 3 } = zeros( Float64 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
        DP_2_policy_lin_with_finale_new( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

end

###############################################################################################################################

function DP_2_policy_lin_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float64 , 3 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] of size ( number_of_allocations )^3 is the value_to_go at the finale (after all allocations are made)
  # Output is the policy (i.e., actions for all states) with linear indices and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        action_lin = zeros( Int8 , div( number_of_allocations * ( number_of_allocations + 1 ) * ( number_of_allocations + 2 ) * ( number_of_allocations + 3 ) - 1 * 2 * 3 * 4 , 24 ) + 1 ) # this is the correct size of the action
        lin_index = 0
        for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0 , number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
            @inbounds begin

                lin_index += 1
                belief_of_success_arm_1 = ( ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 ) )
                belief_of_success_arm_2 = ( ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 ) )

                value_to_go_if_action_1 = belief_of_success_arm_1 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 + 1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] ) + ( 1.0 - belief_of_success_arm_1 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 + 1 ,1+ number_of_successes_arm_2 ]
                value_to_go_if_action_2 = belief_of_success_arm_2 * ( 1.0 + value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 + 1 ] ) + ( 1.0 - belief_of_success_arm_2 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ]

                if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_1
                    action_lin[ lin_index ] = 1 # action 1 (i.e., arm 1)

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_2
                    action_lin[ lin_index ] = 2 # action 2 (i.e., arm 2)

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2
                    action_lin[ lin_index ] = 3 # randomise between actions 1 and 2

                end

            end # @inbounds
        end

        return action_lin , value_to_go[1 ,1 ,1 ]

end

function DP_2_policy_lin( number_of_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Output is the policy (i.e., actions for all states) with linear indices and the Bayes-expected number of successes

        value_to_go :: Array{ Float64 , 3 } = zeros( Float64 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
        DP_2_policy_lin_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

end

function DP_2_policy_lin_with_administration( number_of_allocations :: Int64 , number_of_administered_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms; the better arm (comparing the beliefs) is administered after the trial number_of_administered_allocations times
  # Output is the policy (i.e., actions for all states) with linear indices and the Bayes-expected number of successes

        if number_of_administered_allocations == 0

            DP_2_policy_lin( number_of_allocations , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        else

            # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations )
            value_to_go :: Array{ Float64 , 3 } = zeros( Float64 , number_of_allocations + 1, number_of_allocations + 1, number_of_allocations + 1 )
            # value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] of size ( number_of_allocations )^3 is the value_to_go at the finale (after all allocations are made)
            for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                @inbounds begin
                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = max( belief_of_success_arm_1 , belief_of_success_arm_2 )
                end # @inbounds
            end
            value_to_go .*= number_of_administered_allocations
            DP_2_policy_lin_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        end

end

###############################################################################################################################

function DP_2_policy_with_finale( number_of_allocations :: Int64 , value_to_go :: Array{ Float64 , 3 } , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Output is the policy (i.e., actions for all states) with action[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 , number_of_allocations - number_of_observed_responses ] and the Bayes-expected number of successes

        # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations ) not needed as it is given by value_to_go

        # backwards recursion: t-th step
        action = zeros( Int8 , number_of_allocations , number_of_allocations , number_of_allocations , number_of_allocations )
        for number_of_observed_responses = ( number_of_allocations - 1 ) : -1 : 0 , number_of_successes_arm_2 = 0 : number_of_observed_responses , number_of_failures_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 ) , number_of_successes_arm_1 = 0 : ( number_of_observed_responses - number_of_successes_arm_2 - number_of_failures_arm_1 )
            @inbounds begin

                belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_observed_responses - number_of_successes_arm_1 - number_of_failures_arm_1 )

                value_to_go_if_action_1 = belief_of_success_arm_1 * ( 1 + value_to_go[1+ number_of_successes_arm_1 + 1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] ) + ( 1 - belief_of_success_arm_1 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 + 1 ,1+ number_of_successes_arm_2 ]
                value_to_go_if_action_2 = belief_of_success_arm_2 * ( 1 + value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 + 1 ] ) + ( 1 - belief_of_success_arm_2 ) * value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ]

                if ( value_to_go_if_action_1 - value_to_go_if_action_2 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_1
                    action[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 , number_of_allocations - number_of_observed_responses ] = 1 # action 1 (i.e., arm 1)

                elseif ( value_to_go_if_action_2 - value_to_go_if_action_1 ) > BB_numerical_precision_64 * ( value_to_go_if_action_1 + value_to_go_if_action_2 )

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = value_to_go_if_action_2
                    action[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 , number_of_allocations - number_of_observed_responses ] = 2 # action 2 (i.e., arm 2)

                else #if value_to_go_if_action_1 approx== value_to_go_if_action_2

                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = ( value_to_go_if_action_1 + value_to_go_if_action_2 ) / 2
                    action[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 , number_of_allocations - number_of_observed_responses ] = 3 # randomise between actions 1 and 2

                end

            end # @inbounds
        end

        return action , value_to_go[1 ,1 ,1 ]

end

function DP_2_policy( number_of_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms
  # Output is the policy (i.e., actions for all states) with action[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 , number_of_allocations - number_of_observed_responses ] and the Bayes-expected number of successes

        value_to_go :: Array{ Float64 , 3 } = zeros( Float64 , number_of_allocations + 1 , number_of_allocations + 1 , number_of_allocations + 1 )
        DP_2_policy_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

end

function DP_2_policy_with_administration( number_of_allocations :: Int64 , number_of_administered_allocations :: Int64 , prior_success_arm_1 :: Int64 = Int64( 1 ) , prior_failure_arm_1 :: Int64 = Int64( 1 ) , prior_success_arm_2 :: Int64 = Int64( 1 ) , prior_failure_arm_2 :: Int64 = Int64( 1 ) )
  # This function implements the DP design for 2 arms; the better arm (comparing the beliefs) is administered after the trial number_of_administered_allocations times
  # Output is the policy (i.e., actions for all states) with action[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 , number_of_allocations - number_of_observed_responses ] and the Bayes-expected number of successes

        if number_of_administered_allocations == 0

            DP_2_policy( number_of_allocations , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        else

            # backwards recursion: the finale (i.e., number_of_observed_responses = number_of_allocations )
            value_to_go :: Array{ Float64 , 3 } = zeros( Float64 , number_of_allocations + 1, number_of_allocations + 1, number_of_allocations + 1 )
            # value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] of size ( number_of_allocations )^3 is the value_to_go at the finale (after all allocations are made)
            for number_of_successes_arm_2 = 0 : number_of_allocations , number_of_failures_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 , number_of_successes_arm_1 = 0 : number_of_allocations - number_of_successes_arm_2 - number_of_failures_arm_1
                belief_of_success_arm_1 = ( prior_success_arm_1 + number_of_successes_arm_1 ) / ( prior_success_arm_1 + prior_failure_arm_1 + number_of_successes_arm_1 + number_of_failures_arm_1 )
                belief_of_success_arm_2 = ( prior_success_arm_2 + number_of_successes_arm_2 ) / ( prior_success_arm_2 + prior_failure_arm_2 + number_of_allocations - number_of_successes_arm_1 - number_of_failures_arm_1 )
                @inbounds begin
                    value_to_go[1+ number_of_successes_arm_1 ,1+ number_of_failures_arm_1 ,1+ number_of_successes_arm_2 ] = max( belief_of_success_arm_1 , belief_of_success_arm_2 )
                end # @inbounds
            end
            value_to_go .*= number_of_administered_allocations
            DP_2_policy_with_finale( number_of_allocations , value_to_go , prior_success_arm_1 , prior_failure_arm_1 , prior_success_arm_2 , prior_failure_arm_2 )

        end

end

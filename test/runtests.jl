using BinaryBandit

horizon = 60

BinaryBandit.DP_2_action_lin( horizon , 64 )
BinaryBandit.DP_2_action_lin( horizon , 32 )
BinaryBandit.DP_2_action_lin( horizon , 16 )

BinaryBandit.DP_2_action_lin_with_administration( horizon , 0 , 64 )
BinaryBandit.DP_2_action_lin_with_administration( horizon , 0 , 32 )
BinaryBandit.DP_2_action_lin_with_administration( horizon , 0 , 16 )

BinaryBandit.DP_2_action( horizon , 64 )
BinaryBandit.DP_2_action( horizon , 32 )
BinaryBandit.DP_2_action( horizon , 16 )

BinaryBandit.DP_2_action_with_administration( horizon , 0 , 64 )
BinaryBandit.DP_2_action_with_administration( horizon , 0 , 32 )
BinaryBandit.DP_2_action_with_administration( horizon , 0 , 16 )

BinaryBandit.DP_2_policy_bin_lin( horizon )
BinaryBandit.DP_2_policy_bin_lin_with_administration( horizon , 0 )

BinaryBandit.DP_2_policy_bin( horizon )
BinaryBandit.DP_2_policy_bin_with_administration( horizon , 0 )

BinaryBandit.DP_2_policy_lin_new( horizon )

BinaryBandit.DP_2_policy_lin( horizon )
BinaryBandit.DP_2_policy_lin_with_administration( horizon , 0 )

BinaryBandit.DP_2_policy( horizon )
BinaryBandit.DP_2_policy_with_administration( horizon , 0 )


success_probability_arm_1 = 0.3
success_probability_arm_2 = 0.5

BinaryBandit.DP_2_NS( horizon , success_probability_arm_1 , success_probability_arm_2 )

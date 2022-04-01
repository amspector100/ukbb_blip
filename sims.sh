#!/usr/bin/env
# Common arguments for all simulations
NREPS=1 # make 128
NPROCESSES=1 # more cores = faster simulations
COMMON_ARGS="
        --num_causal [0,2,4,6,8,10]
        --hg2 [0.005,0.0005]
        --reps $NREPS
        --num_processes $NPROCESSES
"

# small p
ARGS1="${COMMON_ARGS}
	--chrome 12
	--start 133000001
"
python3.9 blip_sims.py $ARGS1

# medium p
ARGS2="${COMMON_ARGS}
	--chrome 10
	--start 134000001
"
python3.9 blip_sims.py $ARGS2

# large p
ARGS3="${COMMON_ARGS}
	--chrome 1
	--start 237000001
"
python3.9 blip_sims.py $ARGS3

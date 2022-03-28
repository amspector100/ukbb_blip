import os
import sys
import time
import datetime
import json
import numpy as np
import pandas as pd
from multiprocessing import Pool
from functools import partial

### Multiprocessing helper
def _one_arg_function(list_of_inputs, args, func, kwargs):
	"""
	Globally-defined helper function for pickling in multiprocessing.
	:param list of inputs: List of inputs to a function
	:param args: Names/args for those inputs
	:param func: A function
	:param kwargs: Other kwargs to pass to the function. 
	"""
	new_kwargs = {}
	for i, inp in enumerate(list_of_inputs):
		new_kwargs[args[i]] = inp
	return func(**new_kwargs, **kwargs)


def apply_pool(func, constant_inputs={}, num_processes=1, **kwargs):
	"""
	Spawns num_processes processes to apply func to many different arguments.
	This wraps the multiprocessing.pool object plus the functools partial function. 
	
	Parameters
	----------
	func : function
		An arbitrary function
	constant_inputs : dictionary
		A dictionary of arguments to func which do not change in each
		of the processes spawned, defaults to {}.
	num_processes : int
		The maximum number of processes spawned, defaults to 1.
	kwargs : dict
		Each key should correspond to an argument to func and should
		map to a list of different arguments.
	Returns
	-------
	outputs : list
		List of outputs for each input, in the order of the inputs.
	Examples
	--------
	If we are varying inputs 'a' and 'b', we might have
	``apply_pool(
		func=my_func, a=[1,3,5], b=[2,4,6]
	)``
	which would return ``[my_func(a=1, b=2), my_func(a=3,b=4), my_func(a=5,b=6)]``.
	"""

	# Construct input sequence
	args = sorted(kwargs.keys())
	num_inputs = len(kwargs[args[0]])
	for arg in args:
		if len(kwargs[arg]) != num_inputs:
			raise ValueError(f"Number of inputs differs for {args[0]} and {arg}")
	inputs = [[] for _ in range(num_inputs)]
	for arg in args:
		for j in range(num_inputs):
			inputs[j].append(kwargs[arg][j])

	# Construct partial function
	partial_func = partial(
		_one_arg_function, args=args, func=func, kwargs=constant_inputs,
	)

	# Don't use the pool object if num_processes=1
	num_processes = min(num_processes, len(inputs))
	if num_processes == 1:
		all_outputs = []
		for inp in inputs:
			all_outputs.append(partial_func(inp))
	else:
		with Pool(num_processes) as thepool:
			all_outputs = thepool.map(partial_func, inputs)

	return all_outputs

def rejset_power(rej_sets, beta):
	power = 0
	nfd = 0
	for s in rej_sets:
		if np.any(beta[list(s)] != 0):
			power += 1 / len(s)
		else:
			nfd += 1
	fdp = nfd / max(1, len(rej_sets))
	return power, fdp
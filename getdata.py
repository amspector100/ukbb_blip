import sys
import os
import warnings
import wget
import requests

TRAITS = [
	'body_HEIGHTz', 
	'disease_CARDIOVASCULAR', 
	'biochemistry_HDLcholesterol',
	'biochemistry_LDLdirect'
]

REP_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fnature13835/MediaObjects/41586_2015_BFnature13835_MOESM8_ESM.xls"
REP_FNAME = "41586_2015_BFnature13835_MOESM8_ESM.xls"

def create_dir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)

def pull_main_data():
	"""
	Downloads data necessary for main analysis.
	"""
	file_directory = os.path.dirname(os.path.abspath(__file__))

	# Step 1: download bolt_337K data
	bolt_dir = f"{file_directory}/data/polyfun_results/"
	create_dir(bolt_dir)
	os.chdir(bolt_dir)
	# for each trait download summary statistics
	base_url = f"https://storage.googleapis.com/broad-alkesgroup-public/polyfun_results/"
	for trait in TRAITS:
		fname = f"bolt_337K_unrelStringentBrit_MAF0.001_v3.{trait}.bgen.stats.gz"
		if os.path.exists(fname):
			print(f"Summary statistics for trait={trait} are already downloaded.")
		else:
			print(f"Downloading summary statistics for trait={trait}.")
			wget.download(base_url + fname)

	# Step 2: download SuSiE outputs
	alpha_dir = bolt_dir + "alphas/"
	create_dir(alpha_dir)
	os.chdir(alpha_dir)
	# for each trait download SuSiE output
	for trait in TRAITS:
		fname = f"{trait}.nonfunct.alphas.parquet"
		if os.path.exists(fname):
			print(f"SuSiE model for trait={trait} is already downloaded.")
		else:
			print(f"Downloading SuSiE model for trait={trait}.")
			wget.download(base_url + fname)

	# # Step 3: download data for replication analysis
	# farh_dir = f"{file_directory}/data/farh2015/"
	# create_dir(farh_dir)
	# os.chdir(farh_dir)
	# if os.path.exists(REP_FNAME):
	# 	print("Results from Farh (2015) are already downloaded.")
	# else:
	# 	wget.download(REP_URL)


if __name__ == '__main__':
	pull_main_data()
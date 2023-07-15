import sys
import os
import warnings
import wget
import requests
import gdown

TRAITS = [
	'body_HEIGHTz', 
	'disease_CARDIOVASCULAR', 
	'biochemistry_HDLcholesterol',
	'biochemistry_LDLdirect'
]

# urls for data download
SUSIE_GDRIVE_IDS = {
	"body_HEIGHTz":"1Gkkp_0ExgSyoS_Ip3irWbnz6HIKXcfEz",
	"disease_CARDIOVASCULAR":"1fQYAYZvByvivIM9Y9cBnfX3OJF5i_5Lf",
	"biochemistry_HDLcholesterol":"1WZv0r3Q94WD0HHZ2qey1AtnjbSH8dNpJ",
	"biochemistry_LDLdirect":"1t55ncZw3gtu6K525wJKuahnxE5xnohJn",
}

BOLT_GDRIVE_IDS = {
	"body_HEIGHTz":"1hQNGzAWbrnZxAf6TfyoQ_M0_307-T3jl",
	"disease_CARDIOVASCULAR":"1_mjijbhqgmFHKLZiWQ5FVRj36oINsrf2",
	"biochemistry_HDLcholesterol":"1geE4hwXqBlDETg2ylLNlVx4Kuexnq6-_",
	"biochemistry_LDLdirect":"16egXMkSIscvrO1c5tyzaI1dzzOFfp0X2",
}

LD_DRIVE_IDS = {
	1:{"chr1_237000001_240000001.gz":"1LxLhYZspPFeKDSjXfsfp9_FGL82Dm0rZ",
		"chr1_237000001_240000001.npz":"1JEA2JkjrhKzakiinnaTBpWeqXz74GdK0"},
	10:{"chr10_134000001_137000001.gz":"1WBi0gPz0Jy0xd9ldSXia9ITKBAwaYJB4",
		"chr10_134000001_137000001.npz":"1_MpSco9ejWNF0YzxhpHieLx614kbrSmM"},
	12:{"chr12_133000001_136000001.gz":"1-4a5P6pUBlSdhrCW0Gz2IRS5eb5albay",
		"chr12_133000001_136000001.npz":"1XaOqp5bMa3WWzuVHUFepprNNlmbZ_Swt"},
}

# list of LDs which are now available for free
FREE_LDS = [
    (1, 237000001),
    (10, 134000001),
    (12, 133000001)
]

def data_availability_msg(chrome, start):
    return f"""
        Unfortunately, new since the posting of our paper, the 
        raw data now costs money to download due to its enormous size.
        However, we have made the three raw LD matrices used in our
        simulations available for free, namely:

        {[f'chromosome={x}, start position={y}' for (x,y) in FREE_LDS]}

        where each tuple above lists the chromosome number and starting position of the loci.

        To download the ld matrix for chrome={chrome}, start={start}, 
        visit https://registry.opendata.aws/ukbb-ld/.
    """

# Original URL at which the farh2015 data was hosted
# REP_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fnature13835/MediaObjects/41586_2015_BFnature13835_MOESM8_ESM.xls"
# REP_FNAME = "41586_2015_BFnature13835_MOESM8_ESM.xls"

def create_dir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)

def download_free_ld_data(chrome=1):
	file_directory = os.path.dirname(os.path.abspath(__file__))
	raw_ld_dir = f"{file_directory}/data/ld/"
	create_dir(raw_ld_dir)
	os.chdir(raw_ld_dir)
	for fname in LD_DRIVE_IDS[chrome].keys():
		if os.path.exists(fname):
			print(f"{fname} is already downloaded.")
		else:
			print(f"Downloading {fname}.")
			gdrive_id = LD_DRIVE_IDS[chrome][fname]
			gdown.download(id=snp_gdrive_id, output=fname, quiet=False)
	print("Finished downloading data; now preprocessing.")
	os.chdir(file_directory)

def pull_main_data(download_snps=True, download_susie=True, download_bolt=False):
	"""
	Downloads data necessary for main analysis.
	"""
	file_directory = os.path.dirname(os.path.abspath(__file__))

	# Option 1: download names of ~19.6 million SNPs 
	# included in the analysis
	# (caching these values saves preprocessing)
	create_dir("main_cache/")
	if download_snps:
		os.chdir(f"{file_directory}/main_cache/")
		fname = 'snps.json'
		snp_gdrive_id = "1-qvarB5uzI67IIrToFqdD8f_zD9p4-5l"
		if os.path.exists(fname):
			print("The list of SNPs is already downloaded.")
		else:
			print("Downloading list of SNPs.")
			gdown.download(id=snp_gdrive_id, output=fname, quiet=False)

	# Option 2: download bolt_337K data
	bolt_dir = f"{file_directory}/data/polyfun_results/"
	create_dir(bolt_dir)
	os.chdir(bolt_dir)
	# for each trait download summary statistics
	if download_bolt:
		for trait in TRAITS:
			fname = f"bolt_337K_unrelStringentBrit_MAF0.001_v3.{trait}.bgen.stats.gz"
			if os.path.exists(fname):
				print(f"Summary statistics for trait={trait} are already downloaded.")
			else:
				print(f"Downloading summary statistics for trait={trait}.")
				gdown.download(id=BOLT_GDRIVE_IDS[trait], output=fname, quiet=False)

	# Option 3: download SuSiE outputs
	alpha_dir = bolt_dir + "alphas/"
	create_dir(alpha_dir)
	os.chdir(alpha_dir)
	# for each trait download SuSiE output
	if download_susie:
		for trait in TRAITS:
			fname = f"{trait}.nonfunct.alphas.parquet"
			if os.path.exists(fname):
				print(f"SuSiE model for trait={trait} is already downloaded.")
			else:
				print(f"Downloading SuSiE model for trait={trait}.")
				print(SUSIE_GDRIVE_IDS[trait])
				gdown.download(id=SUSIE_GDRIVE_IDS[trait], output=fname, quiet=False)

	# # Step 3: download data for replication analysis
	# farh_dir = f"{file_directory}/data/farh2015/"
	# create_dir(farh_dir)
	# os.chdir(farh_dir)
	# if os.path.exists(REP_FNAME):
	# 	print("Results from Farh (2015) are already downloaded.")
	# else:
	# 	wget.download(REP_URL)

	# Reset
	os.chdir(file_directory)


if __name__ == '__main__':
	pull_main_data(download_bolt=False)
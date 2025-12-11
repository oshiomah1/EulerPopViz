import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pkg_resources #pip install setuptools
from geovar import *


# Define the base data path
data_path = "/quobyte/bmhenngrp/from-lssc0/"

# Filepath to the VCF file (THIS IS CURRENTLY CHR 21 but it will be expannded to ALL chromosomes when caapap v2 is ready)
#vcf_file = "{}/data/genomes/CAAPA_freeze2_PHASED_common_rare/snps/chr21.20240612_Freeze2.common_rare_phased.snps.vcf.gz".format(data_path)
# vcf_file = "{}data/genomes/2025_CAAPA_Freeze2B/combined/Freeze2.BiallelicSNPs_shapeit5_Phased_common_rare_allchr.vcf.gz".format(data_path)
# # Filepath to the population panel file ( 2 COLUMN tab delimited file with sample id and pop name)
# population_panel = "{}/projects/CAAPA2_functional_annotation/geovar/data/cleaned_population_panel.txt".format(data_path)

# # Reading the population dataframe
# pop_df = read_pop_panel(population_panel)
# pop_df

# # Writing out VCF to a Frequency Table (*Majestiq will use the ouptut from this as side project)
# af_df = vcf_to_freq_table(vcf_file, pop_df=pop_df, outfile="/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/freq_mat_file_allchr.csv", minor_allele=True)

# # Print the beginning of the allele frequency table
# af_df.head()

# # Writing out VCF to a Frequency Table (*Majestiq will use the ouptut from this as side project)
# af_df = vcf_to_freq_table(vcf_file, pop_df=pop_df, outfile="/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/freq_mat_file_v2.csv", minor_allele=True)

# # Print the beginning of the allele frequency table
# af_df.head()


# # File paths
# freq_file = "/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/freq_mat_file_v2.csv"
 
# # Load files
# freq_df = pd.read_csv(freq_file, sep="\t")  # Assuming CSV is comma-separated

# # Creating the GeoVar Object
# geovar_test = GeoVar()


# # Adding in the frequency file (all of it)
# geovar_test.add_freq_mat(freq_mat_file="/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/freq_mat_file.csv")

# #once this is done create a zipped version of the freq_mat_file in the results folder to use downstream


# # Generate a geovar binning with the binning we used in our paper
# geovar_test.geovar_binning()

# # Printing details about the GeoVar object
# print(geovar_test)

# from tqdm import tqdm
# import gzip

# #load in James Kitchens Function from github 
# def geovar_codes_streaming_fixed(geovar_obj, freq_mat_file):
#     """Version of GeoVar code generation algorithm that streams through file to avoid memory overflow.
#     Args:
#         freq_mat_file (:obj:`string`): filepath to
#         frequency table file (see example notebook for formatting).
#     """
#     assert geovar_obj.bins is not None
#     geovar_codes = []
#     # Setting up the testing bins
#     test_bins = np.array([x[1] for x in geovar_obj.bins])
#     with gzip.open(freq_mat_file,'r') as f:
#         header = f.readline()
#         # Take the population labels currently
#         pops = np.array(header.split()[6:])
#         geovar_obj.pops = pops
#         for line in tqdm(f):
#             # Split after the 6th column ...
#             maf_vector = np.array(line.split()[6:]).astype(np.float64)
#             cur_geovar = np.digitize(maf_vector, test_bins, right=True)
#             cur_geovar_code = "".join([str(i) for i in cur_geovar])
#             geovar_codes.append(cur_geovar_code)
#     # Setting the variables here
#     geovar_obj.geovar_codes = np.array(geovar_codes)
#     geovar_obj.n_variants = geovar_obj.geovar_codes.size
#     geovar_obj.n_populations = geovar_obj.pops.size

# geovar = GeoVar(bins=[(0,0), (0,0.01), (0.01,0.05), (0.05,0.1), (0.1,1.0)])
# geovar_codes_streaming_fixed(geovar_obj=geovar, freq_mat_file="/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/freq_mat_file.csv.gz") 
# counts_output = geovar.count_geovar_codes()
# geovar_code_counts = pd.DataFrame({"codes":counts_output[0],"counts":counts_output[1]})
# geovar_code_counts.to_csv("/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/data/caapa_chr21.txt", sep=' ', header=False, index=False)

geovar_test = GeoVar()


# Adding in the frequency file (all of it)
 
geovar_test.add_freq_mat(freq_mat_file="/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/freq_mat_file_allchr.csv")


# Generate a geovar binning with the binning we used in our paper
geovar_test.geovar_binning()

# Printing details about the GeoVar object
print(geovar_test)

from tqdm import tqdm
import gzip

#load in James Kitchens Function from github 
def geovar_codes_streaming_fixed(geovar_obj, freq_mat_file):
    """Version of GeoVar code generation algorithm that streams through file to avoid memory overflow.
    Args:
        freq_mat_file (:obj:`string`): filepath to
        frequency table file (see example notebook for formatting).
    """
    assert geovar_obj.bins is not None
    geovar_codes = []
    # Setting up the testing bins
    test_bins = np.array([x[1] for x in geovar_obj.bins])
    with gzip.open(freq_mat_file,'r') as f:
        header = f.readline()
        # Take the population labels currently
        pops = np.array(header.split()[6:])
        geovar_obj.pops = pops
        for line in tqdm(f):
            # Split after the 6th column ...
            maf_vector = np.array(line.split()[6:]).astype(np.float64)
            cur_geovar = np.digitize(maf_vector, test_bins, right=True)
            cur_geovar_code = "".join([str(i) for i in cur_geovar])
            geovar_codes.append(cur_geovar_code)
    # Setting the variables here
    geovar_obj.geovar_codes = np.array(geovar_codes)
    geovar_obj.n_variants = geovar_obj.geovar_codes.size
    geovar_obj.n_populations = geovar_obj.pops.size

geovar = GeoVar(bins=[(0,0), (0,0.01), (0.01,0.05), (0.05,0.1), (0.1,1.0)])
geovar_codes_streaming_fixed(geovar_obj=geovar, freq_mat_file="/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/freq_mat_file_allchr.csv.gz") 
counts_output = geovar.count_geovar_codes()
geovar_code_counts = pd.DataFrame({"codes":counts_output[0],"counts":counts_output[1]})
geovar_code_counts.to_csv("/quobyte/bmhenngrp/from-lssc0/projects/CAAPA2_functional_annotation/geovar/results/allchr.txt", sep=' ', header=False, index=False)
 
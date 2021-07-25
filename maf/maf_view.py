import os
import pandas as pd


# input_dir = r'D:/junmo/wd/utuc/maf/rmhd_maf/mutect2/sample1/'
input_dir = r'E:/'

# test_sample = r'rmHd_20S-31099-A4-14_mutect.maf'
# test_sample = r'rmHd_20S-14292-A1-7_mutect.maf'
# test_sample = r'rmHd_20S-31099-A4-3_mutect.maf'
test_sample = r'test2.vep.maf'

test_sample_path = input_dir + test_sample


maf_df = pd.read_csv(test_sample_path, sep='\t')

print(list(maf_df.columns))

# print(maf_df)

#  ['t_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count']


print(maf_df[['t_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count']])

# maf_df = maf_df[['t_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count']]

maf_df.dropna(axis=0, inplace=True)

print(maf_df)




import os
from glob import glob
import pandas as pd


# bed file sorting

input_dir = r'/home/jun9485/data/WES/T-only_bed/'
input_format = '*.bed'

input_path_lst = glob(input_dir + input_format)

# print(input_path_lst)

header = ['CHROM', 'START', 'END']

for input_path in input_path_lst:
    bed_df = pd.read_csv(input_path, sep='\t', names=header)
    # print(bed_df.head())
    bed_df.sort_values(by=['CHROM', 'START'], ascending=[True, True], inplace=True)
    # print(bed_df.head())
    print(input_path)
    save_name = input_path.split(r'/')[-1]

    output_path = input_dir + save_name

    bed_df.to_csv(output_path, index=False, header=None, sep="\t")

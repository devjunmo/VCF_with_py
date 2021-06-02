
from glob import glob
import os
import pandas as pd
import sys


input_dir = r'/myData/WES/data/vcf/hard/WES1_210420/Td30_Odx/pass_only_Teratoma_specifics/subsets_GT_trim/'
output_dir = input_dir + r'bed/'
input_format = r'*.csv'
prefix = 'SNP_'

if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)

input_lst = glob(input_dir + input_format)


for i in range(len(input_lst)):
    f_name = input_lst[i].split(r'/')[-1].split(r'.')[0].split(r'_')[-1]
    csv_file = input_lst[i]
    csv_df = pd.read_csv(csv_file, low_memory=False)
    # print(csv_df)
    bed_df = pd.DataFrame(columns=['CHROM', 'START', 'END', 'span'])
    bed_df['CHROM'] = csv_df.pop('CHROM')
    bed_df['START'] = csv_df.pop('POS')
    bed_df['START'] = bed_df['START'] - 1
    bed_df['span'] = csv_df.pop('ALT')
    bed_df['span'] = bed_df['span'].map(len)
    bed_df['END'] = bed_df['START'] + bed_df['span']
    bed_df.pop('span')
    # print(bed_df.head())

    output_path = output_dir + prefix + f_name + '.bed'

    bed_df.to_csv(output_path, index=False, header=None, sep="\t")


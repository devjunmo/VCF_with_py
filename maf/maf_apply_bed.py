
"""
maf file에 bed file을 inner join
bed file의 컬럼명은 대응하는 maf 파일의 컬럼명과 동일하게 설정
"""

import pandas as pd
import os
import natsort
from glob import glob


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


INPUT_DIR = r'E:/UTUC_data/gdc_hg38/maf/2nd/DP_AF_filtered_maf/true_maf'
OUTPUT_DIR_NAME = r'True_positive_maf'

OUTPUT_DIR = os.path.join(INPUT_DIR, OUTPUT_DIR_NAME)

maf_format = '*.maf'
bed_format = '*.txt'

maf_join_col = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
bed_join_col = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']

maf_path_lst = natsort.natsorted(glob(os.path.join(INPUT_DIR, maf_format)))
bed_path_lst = natsort.natsorted(glob(os.path.join(INPUT_DIR, bed_format)))


if len(maf_path_lst) != len(bed_path_lst):
    print(f'{bcolors.FAIL}len(maf_path_lst) != len(bed_path_lst){bcolors.ENDC}')
    exit(1)


def inner_join_maf_bed(_maf_path, _bed_path):
    maf_df = pd.read_csv(_maf_path, sep='\t')
    bed_df = pd.read_csv(_bed_path, sep='\t')

    return pd.merge(maf_df, bed_df, how='inner', \
                    left_on=maf_join_col, right_on=bed_join_col)


for i in range(len(maf_path_lst)):
    
    # true_maf = inner_join_maf_bed(maf_path_lst[i], bed_path_lst[i])

    f_name = os.path.splitext(os.path.basename(maf_path_lst[i]))
    print(f_name)

    break
    # true_maf.to_csv(output_path, index=False, sep='\t')



from glob import glob
import os
from typing import IO
import pandas as pd


# allele count 기준으로 filtering 하기
# normal 2이하 / tumor 5이상 / t-depth 30이상


# colname 생략 없이 출력
pd.set_option('display.max_seq_items', None)
# row 생략 없이 출력
# pd.set_option('display.max_rows', None)
# col 생략 없이 출력
pd.set_option('display.max_columns', None)

input_dir = r'/home/jun9485/data/utuc/annotation/maf/'
input_format = r'*.maf'

output_dir_name = r'filtered_maf/'

output_dir = input_dir + output_dir_name

if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)


input_maf_lst = glob(input_dir + input_format)


for i in range(len(input_maf_lst)):
    maf_file = input_maf_lst[i]

    maf_df = pd.read_csv(maf_file, sep='\t')

    # print(maf_df)

    print(maf_df.columns)

    print(maf_df['t_depth'])
    print(maf_df['t_ref_count'])
    print(maf_df['t_alt_count'])
    print(maf_df['n_ref_count'])
    print(maf_df['n_alt_count']) # 2 미만


    break
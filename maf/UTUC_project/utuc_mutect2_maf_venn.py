import venn
import os
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd


input_dir = r'E:/UTUC_data/WES/maf/mutect2/sample1/'
input_format = r'*.maf'

pair_info = r'E:/UTUC_data/utuc_NT_pair.csv'


pair_df = pd.read_csv(pair_info)
pair_df.set_index('Tumor', inplace=True)
pair_dict = pair_df.to_dict('index')

print(pair_dict)


input_lst = glob(input_dir + input_format)

print(input_lst)


for i in range(len(input_lst)):
    input_maf = input_lst[i]
    print(input_maf)
    t_name = input_maf.split('\\')[-1].split(r'.')[0].split(r'_')[1] # 20S-14292-A1-7
    tumor_grade = pair_dict[t_name]['Tumor_Grade'] # low

    maf_df = pd.read_csv(input_maf, sep='\t', low_memory=False)

    print(maf_df)
    
    
    break
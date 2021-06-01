from glob import glob
import os
import subprocess as sp
import pandas as pd
import numpy as np


# 재활용 x // 단순 스크립트 //


# ref df 수정 (DP30 이상으로..)

INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/'
ref_df_path = r'/myData/IPS_project_sample_lst/WES_Variant_Info_Origin_Teratoma_pair_210513.xlsx'
ref_df = pd.read_excel(ref_df_path)
# print(ref_df)
# print(ref_df["# SNP_common"])


teratoma_specific_dir = INPUT_DIR + r'Teratoma_specifics_pass_only/subsets_DP30/'
origin_specific_dir = INPUT_DIR + r'Origin_specifics/subsets_DP30/'



snp_t_path_lst = glob(teratoma_specific_dir + '*_SNP_*')

# snp_o_path_lst = glob(origin_specific_dir + '*_SNP_*')



# T-only(snp)
for i in range(len(snp_t_path_lst)):
    subset_t = pd.read_csv(snp_t_path_lst[i], low_memory=False)
    # subset_o = pd.read_csv(snp_o_path_lst[i], low_memory=False)

    # print(subset.head())
    sample_name_t = snp_t_path_lst[i].split(r'/')[-1].split(r'.')[0].split('_')[-1]; print(sample_name_t)
    # sample_name_o = snp_o_path_lst[i].split(r'/')[-1].split(r'.')[0].split('_')[-1]; print(sample_name_o)


    # print(sample_name)
    ref_idx_t = ref_df[ref_df['teratoma'] == sample_name_t].index
    # ref_idx_o = ref_df[ref_df['origin cell'] == sample_name_o].index
    
    t_count = subset_t.shape[0]; print(t_count)
    # o_count = subset_o.shape[0]; print(o_count)

    # break


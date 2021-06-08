from glob import glob
import os
import subprocess as sp
import pandas as pd
import numpy as np
from enum import Enum
from class_input_from_csv import MakePairInputList
from class_bcftools_isec import Mk_vcf_intersection

import matplotlib.pyplot as plt


# 이거 안되는 코드 #


# extract_T_only.py로 
# 1) Tonly와 2) T-specific이긴 한데 Origin이 탈락한 부분을 
# 얻고 / 비율구하고 / 걸린 필터 비율 구하기

###################### [HYPER PARAMETERS] ###########################

class TeratomaOrigin(Enum): # csv 파일의 컬럼명 순서로 구성
    Teratoma = 0
    orgin = 1

INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/'
# INPUT_DIR = r'/myData/WES/data/vcf/raw/'
# INPUT_DIR = r'/myData/WES/data/vcf/cnn/WES1_210420/'
# INPUT_DIR = r'/myData/WES/data/vcf/re_hard/'

SNP_INPUT_FORMAT = r'hardFiltered_SNP*.vcf.gz'
INDEL_INPUT_FORMAT = r'hardFiltered_INDEL*.vcf.gz'
# SNP_INPUT_FORMAT = r'SNP*.vcf.gz'
# INDEL_INPUT_FORMAT = r'INDEL*.vcf.gz'
# SNP_INPUT_FORMAT = r'SNP*.vcf.gz'
# INDEL_INPUT_FORMAT = r'INDEL*.vcf.gz'


# input파일의 prefix
# /input/dir/<사 이 부 분 입 력>/name*
PREFIX_SNP_T = r'Teratoma_specifics/SNP_'
PREFIX_SNP_O = 'hardFiltered_SNP_'

PREFIX_INDEL_T = r'Teratoma_specifics/INDEL_'
PREFIX_INDEL_O = 'hardFiltered_INDEL_'


pair_path = r'/myData/WES/src/Origin_Teratoma_pairs.csv'

enum_data = TeratomaOrigin

filter_comp = 'PASS' # raw vcf일땐 '.' // PASS인것만 사용하겠다

is_only_PASS = False

root_output_dir_name_snp = r'snp_isec_Tsp-pass_and_raw_IPS/'
root_output_dir_name_indel = r'indel_isec_Tsp-pass_and_raw_IPS/'

# is_qsub = False

run_mode = "A" 
#####################################################################
 

# step1. isec

if run_mode == "I":
    snp_path = INPUT_DIR + SNP_INPUT_FORMAT
    indel_path = INPUT_DIR + INDEL_INPUT_FORMAT

    obj = MakePairInputList(pair_path, snp_path, indel_path, enum_data)
    obj.trim_pair_df()
    pair_names_df = obj.pair_info_df



    isec_obj = Mk_vcf_intersection(pair_names_df, INPUT_DIR, \
                                    PREFIX_SNP_T, PREFIX_SNP_O, \
                                    PREFIX_INDEL_T, PREFIX_INDEL_O, \
                                    root_output_dir_name_snp, root_output_dir_name_indel, \
                                    _is_only_PASS = is_only_PASS, _filter_comp = filter_comp)

    isec_obj.run_isec()


# step2. gather_specific_isec_data.py 활용
# T_only_data / bath_T_and_I_pos_IPS_view 디렉토리에 각각 생성


# step3. bcftools_mk_subset.py를 활용하여 csv형태로 만들어줌


# step4. csv파일들을 가져와 데이터 프레임으로 만든 후 원하는 데이터를 뽑아내기

elif run_mode == "A":

    ref_df_path = r'/myData/IPS_project_sample_lst/WES_Variant_Info_Origin_Teratoma_pair_210513.xlsx'
    ref_df = pd.read_excel(ref_df_path)
    # print(ref_df)
    # print(ref_df["# SNP_common"])

    t_only_dir = r'/myData/WES/data/vcf/hard/WES1_210420/T_only_data/subsets/'
    droped_ips_dir = r'/myData/WES/data/vcf/hard/WES1_210420/both_T_and_I_pos_IPS_view/subsets/'

    snp_t_only_path_lst = glob(t_only_dir + '*_SNP_*')
    # print(t_only_path_lst)
    snp_dropped_ips_path_lst = glob(droped_ips_dir + '*_SNP_*')
    # print(droped_ips_path_lst)

    # T-only(snp)
    for i in range(len(snp_t_only_path_lst)):
        subset = pd.read_csv(snp_t_only_path_lst[i], low_memory=False)
        # print(subset.head())
        sample_name = snp_t_only_path_lst[i].split(r'/')[-1].split(r'.')[0].split('_')[-1]
        # print(sample_name)
        ref_idx = ref_df[ref_df['teratoma'] == sample_name].index
        t_only_count = subset.shape[0]
        # print(t_only_count)
        ref_df.loc[ref_idx, 'SNP_T-only'] = t_only_count
        ref_df.loc[ref_idx, 'SNP_T-only %'] = round((t_only_count / (ref_df.loc[ref_idx, '# SNP_teratoma-specific'])) * 100, 2)
        # break

    

    # def plot_generator(): # 값list, 빈도list를 
    #     yield 0
    
    # g = plot_generator()

    # dropped ips(snp)
    for i in range(len(snp_dropped_ips_path_lst)):
        subset = pd.read_csv(snp_dropped_ips_path_lst[i], low_memory=False)
        sample_name = snp_dropped_ips_path_lst[i].split(r'/')[-1].split(r'.')[0].split('_')[-1]
        print(sample_name)
        ref_idx = ref_df[ref_df['origin cell'] == sample_name].index
        dropped_ips_count = subset.shape[0]
        # print(dropped_ips_count)
        ref_df.loc[ref_idx, 'SNP_dropped_ips'] = dropped_ips_count

        # break


    print(ref_df)
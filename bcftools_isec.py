import glob
import os
import subprocess as sp
import pandas as pd
import numpy as np
from enum import Enum
from class_input_from_csv import MakePairInputList


# os.system('clear')

# 1. 페어링된 샘플 리스트 불러오기
# 2. \n단위로 읽어서 리스트에 넣어주기
# 3. d1 d2 변수에 할당해주기

# 4. bcftools isec -p {output_dir} d1 d2

###################### [HYPER PARAMETERS] ###########################

## 파라미터 조정 후 밑에 인풋 갯수만큼 함수만 추가해주기 ##

# 자주쓴다면.. 아에 하이퍼 파라미터를 리스트로 한꺼번에 받아서
# 파라미터만 보고 쓰게 수정이 가능하긴 함

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
PREFIX_SNP = 'hardFiltered_SNP_'
PREFIX_INDEL = 'hardFiltered_INDEL_'
# PREFIX_SNP = 'SNP_'
# PREFIX_INDEL = 'INDEL_'
# PREFIX_SNP = 'SNP_cnn_'
# PREFIX_INDEL = 'INDEL_cnn_'

pair_path = r'/myData/WES/src/Origin_Teratoma_pairs.csv'

enum_data = TeratomaOrigin

filter_comp = 'PASS' # raw vcf일땐 '.' // PASS인것만 사용하겠다

is_only_PASS = False

root_output_dir_name_snp = r'snp_isec_include_filtered/'
root_output_dir_name_indel = r'indel_isec_include_filtered/'

is_qsub = False

#####################################################################

snp_path = INPUT_DIR + SNP_INPUT_FORMAT
indel_path = INPUT_DIR + INDEL_INPUT_FORMAT

obj = MakePairInputList(pair_path, snp_path, indel_path, enum_data)
obj.trim_pair_df()
pair_names_df = obj.pair_info_df

# print(pair_names_df)

# exit(0)

for rows in pair_names_df.itertuples():

    snp_teratoma_data_path = INPUT_DIR + PREFIX_SNP + rows[1] + '.vcf.gz'
    snp_origin_data = INPUT_DIR + PREFIX_SNP + rows[2] + '.vcf.gz'

    indel_teratoma_data_path = INPUT_DIR + PREFIX_INDEL + rows[1] + '.vcf.gz'
    indel_origin_data = INPUT_DIR + PREFIX_INDEL + rows[2] + '.vcf.gz'

    output_dir_name = rows[1] + '_' + rows[2]

    snp_output_dir = INPUT_DIR + root_output_dir_name_snp + output_dir_name
    indel_output_dir = INPUT_DIR + root_output_dir_name_indel + output_dir_name

    if os.path.isdir(snp_output_dir) is False:
        if os.path.isdir(INPUT_DIR + root_output_dir_name_snp) is False:
            os.mkdir(INPUT_DIR + root_output_dir_name_snp)
        os.mkdir(snp_output_dir)
    if os.path.isdir(indel_output_dir) is False:
        if os.path.isdir(INPUT_DIR + root_output_dir_name_indel) is False:
            os.mkdir(INPUT_DIR + root_output_dir_name_indel)
        os.mkdir(indel_output_dir)
    

    if is_only_PASS:
        sp.call(f"bcftools isec -f {filter_comp} -p {snp_output_dir}/ {snp_teratoma_data_path} {snp_origin_data}", shell=True)
        sp.call(f"bcftools isec -f {filter_comp} -p {indel_output_dir}/ {indel_teratoma_data_path} {indel_origin_data}", shell=True)
    
    else:
        sp.call(f"bcftools isec -p {snp_output_dir}/ {snp_teratoma_data_path} {snp_origin_data}", shell=True)
        sp.call(f"bcftools isec -p {indel_output_dir}/ {indel_teratoma_data_path} {indel_origin_data}", shell=True)

    

#################################################


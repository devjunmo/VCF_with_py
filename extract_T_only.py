import glob
import os
import subprocess as sp
import pandas as pd
import numpy as np
from enum import Enum
from class_input_from_csv import MakePairInputList
from class_bcftools_isec import Mk_vcf_intersection


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

isec_obj = Mk_vcf_intersection(pair_names_df, INPUT_DIR, PREFIX_SNP, PREFIX_INDEL, \
                                root_output_dir_name_snp, root_output_dir_name_indel, \
                                _is_only_PASS = is_only_PASS, _filter_comp = filter_comp)

isec_obj.run_isec()
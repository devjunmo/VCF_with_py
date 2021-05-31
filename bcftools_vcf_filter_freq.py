
import pandas as pd
import os
from class_input_from_csv import MakePairInputList
from enum import Enum
import subprocess as sp


class TeratomaOrigin(Enum): # csv 파일의 컬럼명 순서로 구성
    Teratoma = 0
    orgin = 1

# INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/'
# INPUT_DIR = r'/myData/WES/data/vcf/raw/'
# INPUT_DIR = r'/myData/WES/data/vcf/cnn/WES1_210420/'
INPUT_DIR = r'/myData/WES/data/vcf/re_hard/'

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

pair_file_name = r'IPS_Terratoma_pairs.csv'

enum_data = TeratomaOrigin

filter_comp = 'PASS' # raw vcf일땐 '.' // PASS인것만 사용하겠다

is_qsub = False

#####################################################################

pair_path = INPUT_DIR + pair_file_name
snp_path = INPUT_DIR + SNP_INPUT_FORMAT
indel_path = INPUT_DIR + INDEL_INPUT_FORMAT

obj = MakePairInputList(pair_path, snp_path, indel_path, enum_data)
obj.trim_pair_df()
pair_names_df = obj.pair_info_df

for rows in pair_names_df.itertuples():
    snp_teratoma_data_path = INPUT_DIR + PREFIX_SNP + rows[1] + '.vcf.gz'
    snp_origin_data = INPUT_DIR + PREFIX_SNP + rows[2] + '.vcf.gz'

    indel_teratoma_data_path = INPUT_DIR + PREFIX_INDEL + rows[1] + '.vcf.gz'
    indel_origin_data = INPUT_DIR + PREFIX_INDEL + rows[2] + '.vcf.gz'

    sp.call(f"bcftools query -f '%CHROM %POS %FILTER\n' {snp_teratoma_data_path} > {INPUT_DIR}testT.sub", shell=True)
    sp.call(f"bcftools query -f '%CHROM %POS %FILTER\n' {snp_origin_data} > {INPUT_DIR}testO.sub", shell=True)
    
    df_text = pd.read_csv(f'{INPUT_DIR}testT.sub', sep=' ', header=None, \
        names=['CHROM', 'POS', 'FILTER'])
    df_text2 = pd.read_csv(f'{INPUT_DIR}testO.sub', sep=' ', header=None, \
        names=['CHROM', 'POS', 'FILTER'])

    print('--------------------------------------------------')
    print(rows[1])
    print('--------------------------------------------------')
    print(df_text.groupby(['FILTER']).count())
    print('--------------------------------------------------')
    print(rows[2])
    print('--------------------------------------------------')
    print(df_text2.groupby(['FILTER']).count())




# 샘플별로 subset만들고, 
# 지정한 디렉토리에 저장

## isec -> gather후 subset 만드는 순서로 적용


import pandas as pd
import os
from class_input_from_csv import MakePairInputList
from enum import Enum
import subprocess as sp
import glob



# runmode #

### with_pair_list: 페어리스트에 등록된 요소만 돌리기, IPS, Teratoma 둘 다 일때
### T_only: 테라토마만 있을때, 디렉토리 내 모든 파일 돌리기 

# RUN_MODE = 'with_pair_list' 
RUN_MODE = 'specific_only' 


class TeratomaOrigin(Enum): # csv 파일의 컬럼명 순서로 구성
    # Teratoma = 0
    orgin = 1

# INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/'
# INPUT_DIR = r'/myData/WES/data/vcf/raw/'
# INPUT_DIR = r'/myData/WES/data/vcf/cnn/WES1_210420/'
# INPUT_DIR = r'/myData/WES/data/vcf/re_hard/'
# INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/Teratoma_specifics/'
# INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/Origin_specifics/'
# INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/commons/'
# INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/incFilter_Teratoma_specifics/'
# INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/incFilter_Origin_specifics/'
# INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/incFilter_commons_T/'
INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/incFilter_commons_O/'


# SNP_INPUT_FORMAT = r'hardFiltered_SNP*.vcf.gz'
# INDEL_INPUT_FORMAT = r'hardFiltered_INDEL*.vcf.gz'
SNP_INPUT_FORMAT = r'SNP_*.vcf.gz'
INDEL_INPUT_FORMAT = r'INDEL_*.vcf.gz'
# SNP_INPUT_FORMAT = r'SNP*.vcf.gz'
# INDEL_INPUT_FORMAT = r'INDEL*.vcf.gz'


# input파일의 prefix
# PREFIX_SNP = 'hardFiltered_SNP_'
# PREFIX_INDEL = 'hardFiltered_INDEL_'
PREFIX_SNP = 'SNP_'
PREFIX_INDEL = 'INDEL_'
# PREFIX_SNP = 'SNP_cnn_'
# PREFIX_INDEL = 'INDEL_cnn_'


PREFIX_SAVE_SNP = 'subset_SNP_'
PREFIX_SAVE_INDEL = 'subset_INDEL_'
# PREFIX_SAVE_SNP = 'raw_subset_SNP_'
# PREFIX_SAVE_INDEL = 'raw_subset_INDEL_'


SAVE_DIR_NAME = 'subsets/'
# SAVE_DIR_NAME = 'raw_subsets/'

# T-specific일때 IPS파일이 없어서 path list가 모두 삭제되는 현상 발생 해결 목적
pair_path = r'/myData/WES/src/Origin_Teratoma_pairs.csv'

enum_data = TeratomaOrigin


# filter_comp = 'PASS' # raw vcf일땐 '.' // PASS인것만 사용하겠다

is_qsub = False

# subset_comp = "'%CHROM\t%POS\t%DP\t%MQ\t%QD\t%SOR\n'" # 이방식 안됨. 어쩔수 없음 할때마다 코드 내부의 sp.call부분 손봐야함
# header = ['CHROM', 'POS', 'DP', 'MQ', 'QD', 'SOR', 'GT']
header = ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO', 'GT', 'AD', 'DP', 'GQ', 'PL']

#####################################################################

snp_path = INPUT_DIR + SNP_INPUT_FORMAT
indel_path = INPUT_DIR + INDEL_INPUT_FORMAT

SAVE_DIR = INPUT_DIR + SAVE_DIR_NAME

if os.path.isdir(SAVE_DIR) is False:
    os.mkdir(SAVE_DIR)


def drop_and_reset(__df, __idx):
    __df.drop(__idx, inplace=True)
    __df.reset_index(inplace=True, drop=True)


def rm_GT_more_than_one(_df):
    idx = _df[(_df["GT"] != "0/1") & (_df["GT"] != "1/1")].index
    drop_and_reset(_df, idx)


def rm_DP_less_than_30(_df):
    idx = _df[_df["DP"] < 30].index
    drop_and_reset(_df, idx)


# print('???\n', pair_names_df)

# exit(0)
if RUN_MODE == 'with_pair_list':
    obj = MakePairInputList(pair_path, snp_path, indel_path, enum_data)
    obj.trim_pair_df()
    pair_names_df = obj.pair_info_df

    # print(pair_names_df)
    for rows in pair_names_df.itertuples():
        snp_teratoma_data_path = INPUT_DIR + PREFIX_SNP + rows[1] + '.vcf.gz'
        snp_origin_data = INPUT_DIR + PREFIX_SNP + rows[2] + '.vcf.gz'

        indel_teratoma_data_path = INPUT_DIR + PREFIX_INDEL + rows[1] + '.vcf.gz'
        indel_origin_data = INPUT_DIR + PREFIX_INDEL + rows[2] + '.vcf.gz'

        sp.call(f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT\t%AD\t%DP\t%GQ\t%PL]\n' {snp_teratoma_data_path} > {INPUT_DIR}snpT.tmp", shell=True)
        sp.call(f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT\t%AD\t%DP\t%GQ\t%PL]\n' {snp_origin_data} > {INPUT_DIR}snpO.tmp", shell=True)
        sp.call(f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT\t%AD\t%DP\t%GQ\t%PL]\n' {indel_teratoma_data_path} > {INPUT_DIR}indT.tmp", shell=True)
        sp.call(f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT\t%AD\t%DP\t%GQ\t%PL]\n' {indel_origin_data} > {INPUT_DIR}indO.tmp", shell=True)

        # sp.call(f"rm {INPUT_DIR}*.tmp", shell=True)
        
        df_text1 = pd.read_csv(f'{INPUT_DIR}snpT.tmp', sep='\t', header=None, \
            names=header, low_memory=False)
        df_text2 = pd.read_csv(f'{INPUT_DIR}snpO.tmp', sep='\t', header=None, \
            names=header, low_memory=False)
        df_text3 = pd.read_csv(f'{INPUT_DIR}indT.tmp', sep='\t', header=None, \
            names=header, low_memory=False)
        df_text4 = pd.read_csv(f'{INPUT_DIR}indO.tmp', sep='\t', header=None, \
            names=header, low_memory=False)

        df_lst = [df_text1, df_text2, df_text3, df_text4]

        for vcf_df in df_lst:
            rm_GT_more_than_one(vcf_df)
            rm_DP_less_than_30(vcf_df)
            exit(0)

        # VCF handling code 삽입 filter_comp 변수 활용!!

        df_text1.to_csv(f'{SAVE_DIR}{PREFIX_SAVE_SNP + rows[1]}.csv', header=True, index=False)
        df_text2.to_csv(f'{SAVE_DIR}{PREFIX_SAVE_SNP + rows[2]}.csv', header=True, index=False)
        df_text3.to_csv(f'{SAVE_DIR}{PREFIX_SAVE_INDEL + rows[1]}.csv', header=True, index=False)
        df_text4.to_csv(f'{SAVE_DIR}{PREFIX_SAVE_INDEL + rows[2]}.csv', header=True, index=False)


elif RUN_MODE == 'specific_only':
    obj = MakePairInputList(pair_path, snp_path, indel_path, enum_data)
    obj.trim_pair_df()
    pair_names_df = obj.pair_info_df
    # exit(0)

    # print(pair_names_df)
    for rows in pair_names_df.itertuples():
        snp_teratoma_data_path = INPUT_DIR + PREFIX_SNP + rows[1] + '.vcf.gz'
        indel_teratoma_data_path = INPUT_DIR + PREFIX_INDEL + rows[1] + '.vcf.gz'

        sp.call(f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT\t%AD\t%DP\t%GQ\t%PL]\n' {snp_teratoma_data_path} > {INPUT_DIR}snpT.tmp", shell=True)
        sp.call(f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO[\t%GT\t%AD\t%DP\t%GQ\t%PL]\n' {indel_teratoma_data_path} > {INPUT_DIR}indT.tmp", shell=True)

        # sp.call(f"rm {INPUT_DIR}*.tmp", shell=True)
        
        df_text1 = pd.read_csv(f'{INPUT_DIR}snpT.tmp', sep='\t', header=None, \
            names=header, low_memory=False)

        df_text3 = pd.read_csv(f'{INPUT_DIR}indT.tmp', sep='\t', header=None, \
            names=header, low_memory=False)

        df_lst = [df_text1, df_text3]

        for vcf_df in df_lst:
            rm_GT_more_than_one(vcf_df)
            rm_DP_less_than_30(vcf_df)

        # VCF handling code 삽입 filter_comp 변수 활용!!

        df_text1.to_csv(f'{SAVE_DIR}{PREFIX_SAVE_SNP + rows[1]}.csv', header=True, index=False)
        df_text3.to_csv(f'{SAVE_DIR}{PREFIX_SAVE_INDEL + rows[1]}.csv', header=True, index=False)



elif RUN_MODE == 'T_only':

    snp_list = glob.glob(snp_path)
    indel_list = glob.glob(indel_path)

    if len(snp_list) != len(indel_list):
        print('error: snp목록과 indel목록의 갯수가 일치하지 않음')
        exit(1)
    
    for i in range(len(snp_list)):
        sample_name = snp_list[i].split(r'/')[-1].split(r'.')[0].split(r'_')[-1]
        sp.call(f"bcftools query -f '%CHROM\t%POS\t%DP\t%MQ\t%QD\t%SOR[\t%GT]\n' {snp_list[i]} > {INPUT_DIR}snpT.tmp", shell=True)
        sp.call(f"bcftools query -f '%CHROM\t%POS\t%DP\t%MQ\t%QD\t%SOR[\t%GT]\n' {indel_list[i]} > {INPUT_DIR}indT.tmp", shell=True)

        df_text1 = pd.read_csv(f'{INPUT_DIR}snpT.tmp', sep='\t', header=None, \
            names=header, index_col=None)
        df_text2 = pd.read_csv(f'{INPUT_DIR}indT.tmp', sep='\t', header=None, \
            names=header, index_col=None)

        # df_text1 = df_text1[df_text1['FILTER'] == filter_comp]
        # df_text2 = df_text2[df_text2['FILTER'] == filter_comp]

        rm_GT_more_than_one(df_text1)
        rm_GT_more_than_one(df_text2)
        
        df_text1.to_csv(f'{SAVE_DIR}{PREFIX_SAVE_SNP + sample_name}.csv', header=True, index=False)
        df_text2.to_csv(f'{SAVE_DIR}{PREFIX_SAVE_INDEL + sample_name}.csv', header=True, index=False)

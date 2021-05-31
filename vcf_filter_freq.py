# filter되어 없어진 값들의 필터된 빈도를 출력하는 코드

import pandas as pd
import os
from class_input_from_csv import MakePairInputList
from enum import Enum
import vcfpy


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

# print(pair_names_df)

# exit(0)


def set_property_from_col(_record, _columns):
    tmp_lst = []
    for i in range(len(columns)):
        if _columns[i] == "CHROM":
            tmp_lst.append(_record.CHROM)
        elif _columns[i] == "POS":
            tmp_lst.append(_record.POS)
        elif _columns[i] == "FILTER":
            tmp_lst.append(_record.FILTER)
    return tmp_lst



def record2df(_recordObj, _dfHeaderList):
    vcf_df = pd.DataFrame(columns=_dfHeaderList)

    count = 0
    for record in _recordObj:
        line = set_property_from_col(record, _dfHeaderList)
        vcf_df.loc[str(count)] = line
        count += 1
        if count % 10000 == 0:
            print(f'{count}row processed...')
    return vcf_df



columns=['CHROM', 'POS', 'FILTER']

for rows in pair_names_df.itertuples():

    snp_teratoma_data_path = INPUT_DIR + PREFIX_SNP + rows[1] + '.vcf.gz'
    snp_origin_data = INPUT_DIR + PREFIX_SNP + rows[2] + '.vcf.gz'

    indel_teratoma_data_path = INPUT_DIR + PREFIX_INDEL + rows[1] + '.vcf.gz'
    indel_origin_data = INPUT_DIR + PREFIX_INDEL + rows[2] + '.vcf.gz'


    reader_snp_te = vcfpy.Reader.from_path(snp_teratoma_data_path)
    reader_snp_ip = vcfpy.Reader.from_path(snp_origin_data)
    reader_indel_te = vcfpy.Reader.from_path(indel_teratoma_data_path)
    reader_indel_ip = vcfpy.Reader.from_path(indel_origin_data)

    print(record2df(reader_snp_te, columns).groupby(['FILTER']).count())
    print(record2df(reader_snp_ip, columns).groupby(['FILTER']).count())
    print(record2df(reader_indel_te, columns).groupby(['FILTER']).count())
    print(record2df(reader_indel_ip, columns).groupby(['FILTER']).count())





    








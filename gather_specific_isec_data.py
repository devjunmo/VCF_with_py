
# 특정 디렉토리에 isec으로 부터 나온 specific or common data를 디렉토리 명을 이용하여 이름 붙여 저장하기
# 만들어진 데이터 그룹은 bcftools_mk_subset.py로 서브셋화 시켜서 R로 플로팅

import os
from enum import Enum
import shutil
import subprocess as sp

######################################################################################################

# ~ _specifics 디렉토리 만들어 줘야 함. 향후 자동 생성으로 코드 수정할것 

# 완성 후 vcf -> vcf.gz으로 변환해 줘야 함. sp.call로 자동화 가능.

class TeratomaOrigin(Enum): # csv 파일의 컬럼명 순서로 구성
    Teratoma = 0
    orgin = 1


# RUN_TYPE = 'SNP'
RUN_TYPE = 'INDEL'


# SNP_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/snp_isec/'
# INDEL_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/indel_isec/'
# SNP_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/snp_isec_include_filtered/'
# INDEL_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/indel_isec_include_filtered/'
# SNP_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/snp_isec_Tsp-pass_and_raw_IPS/'
# INDEL_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/indel_isec_Tsp-pass_and_raw_IPS/'
# SNP_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/snp_isec_pass_only/'
# INDEL_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/indel_isec_pass_only/'
# OUTPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/Teratoma_specifics/'
SNP_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/Td30_Odx/snp_isec_pass_only/'
INDEL_INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/Td30_Odx/indel_isec_pass_only/'


# vcf_name = '0000.vcf' # 0000.vcf : teratoma // 0001.vcf : origin // 0002.vcf : T_common // 0003.vcf: O_common
# vcf_name = '0001.vcf'
# vcf_name = '0002.vcf'
vcf_name = '0003.vcf'

OUTPUT_ROOT = r'/myData/WES/data/vcf/hard/WES1_210420/Td30_Odx/'

# SPECIFIC_DIR = 'Teratoma_specifics/'
# SPECIFIC_DIR = 'Origin_specifics/'
# SPECIFIC_DIR = 'commons/'
# SPECIFIC_DIR = 'incFilter_Teratoma_specifics/'
# SPECIFIC_DIR = 'incFilter_Origin_specifics/'
# SPECIFIC_DIR = 'incFilter_Commons_T/'
# SPECIFIC_DIR = 'T_only_data/'
# SPECIFIC_DIR = 'pass_only_Teratoma_specifics/'
# SPECIFIC_DIR = 'pass_only_Origin_specifics/'
# SPECIFIC_DIR = 'pass_only_T_commons/'
SPECIFIC_DIR = 'pass_only_O_commons/'

# SPECIFIC_DIR = 'both_T_and_I_pos_IPS_view/'

OUTPUT_SPECIFIC_DIR = OUTPUT_ROOT + SPECIFIC_DIR

enum_data = TeratomaOrigin

# 얘 잊지마..
# target_sample_name = enum_data.Teratoma.name
target_sample_name = enum_data.orgin.name 

######################################################################################################

# 딕셔너리 사용해서 key = name / value = path로 저장하자

path_dict = dict()

if os.path.isdir(OUTPUT_SPECIFIC_DIR) is False:
    os.mkdir(OUTPUT_SPECIFIC_DIR)


def get_dir_name(_path, _target_sample_type):
    if _target_sample_type == enum_data.Teratoma.name:
        _sample_name = _path.split(r'/')[-2].split(r'_')[0]
    elif _target_sample_type == enum_data.orgin.name:
        _sample_name = _path.split(r'/')[-2].split(r'_')[1]
    return _sample_name


def mk_path_dict(root_dir, prefix, _vcf_name, _path_dict, _target_sample_type):
    files = os.listdir(root_dir)
    for _file in files:
        path = os.path.join(root_dir, _file)
        if os.path.isdir(path):
            mk_path_dict(path, prefix, _vcf_name, _path_dict, _target_sample_type) # 디렉토리면 재귀호출
        else:
            if _file == _vcf_name:
                _name = get_dir_name(path, _target_sample_type)
                path_dict[_name] = path



if RUN_TYPE == 'SNP':
    mk_path_dict(SNP_INPUT_DIR, '', vcf_name, path_dict, target_sample_name)

    for path in path_dict.items():
        f_name = path[0]
        f_path = path[1] 
        shutil.copy(f_path, OUTPUT_SPECIFIC_DIR + 'SNP_' + f_name + '.vcf')
    os.chdir(OUTPUT_SPECIFIC_DIR)
    sp.call(r"find . -maxdepth 1 -name 'SNP*.vcf' -exec bgzip {} {}.gz \;", shell=True)
    sp.call(r"find . -maxdepth 1 -name 'SNP*.vcf.gz' -exec tabix -p vcf {} \;", shell=True)


elif RUN_TYPE == 'INDEL':
    mk_path_dict(INDEL_INPUT_DIR, '', vcf_name, path_dict, target_sample_name)
    for path in path_dict.items():
        f_name = path[0]
        f_path = path[1] 
        shutil.copy(f_path, OUTPUT_SPECIFIC_DIR + 'INDEL_' + f_name + '.vcf')
    os.chdir(OUTPUT_SPECIFIC_DIR)
    sp.call(r"find . -maxdepth 1 -name 'INDEL*.vcf' -exec bgzip {} {}.gz \;", shell=True)
    sp.call(r"find . -maxdepth 1 -name 'INDEL*.vcf.gz' -exec tabix -p vcf {} \;", shell=True)



import os
from enum import Enum
import shutil
import subprocess as sp

######################################################################################################

# WGS, WES 비교용

class TeratomaOrigin(Enum): # 여기서 이거 필요없음. 한큐에 Teratoma, origin 다 됨
    Teratoma = 0
    orgin = 1


# RUN_TYPE = 'SNP'
RUN_TYPE = 'INDEL'



# SNP_INPUT_DIR = r'/data_244/VCF/WGS_WES_isec/snp_isec_pass_only/'
# INDEL_INPUT_DIR = r'/data_244/VCF/WGS_WES_isec/indel_isec_pass_only/'
# SNP_INPUT_DIR = r'/data_244/VCF/WGS_WES_isec_interval_apply/snp_isec_pass_only/'
# INDEL_INPUT_DIR = r'/data_244/VCF/WGS_WES_isec_interval_apply/indel_isec_pass_only/'
SNP_INPUT_DIR = r'/data_244/VCF/noDP_WGS_WES_isec_interval_apply/snp_isec_pass_only/'
INDEL_INPUT_DIR = r'/data_244/VCF/noDP_WGS_WES_isec_interval_apply/indel_isec_pass_only/'

# 0000.vcf : WES-spe // 0001.vcf : WGS-spe // 0002.vcf : WES-common // 0003.vcf: WGS-common
# vcf_name = '0000.vcf' 
# vcf_name = '0001.vcf'
vcf_name = '0002.vcf'
# vcf_name = '0003.vcf'

OUTPUT_ROOT = r'/data_244/VCF/gatherd_noDP_WGS_WES_interval_apply/'

# SPECIFIC_DIR = 'WES_specific/'
# SPECIFIC_DIR = 'WGS_specific/'
SPECIFIC_DIR = 'common/'

OUTPUT_SPECIFIC_DIR = OUTPUT_ROOT + SPECIFIC_DIR

enum_data = TeratomaOrigin

# 얘 필요없음. 액팅안함
# target_sample_name = enum_data.Teratoma.name
target_sample_name = enum_data.orgin.name 

######################################################################################################

# 딕셔너리 사용해서 key = name / value = path로 저장하자

path_dict = dict()

if os.path.isdir(OUTPUT_ROOT) is False:
    os.mkdir(OUTPUT_ROOT)

if os.path.isdir(OUTPUT_SPECIFIC_DIR) is False:
    os.mkdir(OUTPUT_SPECIFIC_DIR)


# 여기가 데이터 디렉토리에 따라 살짝식 다름.. 좋은 코드라고 할수 없긴 함

def get_dir_name(_path, _target_sample_type):
    _sample_name = _path.split(r'/')[-2]
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


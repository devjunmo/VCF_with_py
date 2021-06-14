from glob import glob
import os
import subprocess as sp
import pandas as pd
import numpy as np
from enum import Enum
from class_input_from_csv import MakePairInputList
from class_bcftools_isec import Mk_vcf_intersection



# path 딕셔너리 생성 
# key: hard~snp/indel~T/I   val: [WESpath, WGSpath]

wgs_vcf_dir = r'/data_244/WGS/HN00146173/gs/hardfiltered/'
wes_vcf_dir = r'/data_244/WES/hardfilterd_vcf/dp_apply/'

isec_output_dir = r'/data_244/VCF/WGS_WES_isec/'


input_format = 'hardFiltered*.vcf.gz'

path_dict = dict()

root_output_dir_name_snp = r'snp_isec_pass_only/'
root_output_dir_name_indel = r'indel_isec_pass_only/'

is_PASS_only = True
filter_comp = 'PASS'

output_root_dir = r''

input_wgs_path_lst = glob(wgs_vcf_dir + input_format)
# input_wgs_path_lst = [wgs_flag, input_wgs_path_lst]

input_wes_path_lst = glob(wes_vcf_dir + input_format)
# input_wes_path_lst = [wes_flag, input_wes_path_lst]


# dict val = [wes_path, wgs_path]

for wes_input in input_wes_path_lst:
    key_data = wes_input.split(r'/')[-1].split('.')[0]
    path_dict[key_data] = [wes_input]

for wgs_input in input_wgs_path_lst:
    key_data = wgs_input.split(r'/')[-1].split('.')[0]
    try:
        path_dict[key_data].append(wgs_input)
    except KeyError:
        # del path_dict[key_data]
        continue

rm_key_lst = []

for k, v in path_dict.items():
    if len(v) < 2:
        rm_key_lst.append(k)

for rm_key in rm_key_lst:
    del path_dict[rm_key]



for k, v in path_dict.items():
    output_dir_name = k.split(r'_')[-1]
    snp_output_dir = isec_output_dir + root_output_dir_name_snp + output_dir_name
    indel_output_dir = isec_output_dir + root_output_dir_name_indel + output_dir_name

    wes_path = v[0]
    wgs_path = v[1]

    var_type = k.split(r'_')[1]

    if os.path.isdir(snp_output_dir) is False:
        if os.path.isdir(isec_output_dir + root_output_dir_name_snp) is False:
            os.mkdir(isec_output_dir + root_output_dir_name_snp)
        os.mkdir(snp_output_dir)
    if os.path.isdir(indel_output_dir) is False:
        if os.path.isdir(isec_output_dir + root_output_dir_name_indel) is False:
            os.mkdir(isec_output_dir + root_output_dir_name_indel)
        os.mkdir(indel_output_dir)
    
    if is_PASS_only:
        if var_type == 'SNP':
            sp.call(f"bcftools isec -f {filter_comp} -p {snp_output_dir}/ {wes_path} {wgs_path}", shell=True)
        elif var_type == 'INDEL':
            sp.call(f"bcftools isec -f {filter_comp} -p {indel_output_dir}/ {wes_path} {wgs_path}", shell=True)

    # else:
    #     sp.call(f"bcftools isec -p {snp_output_dir}/ {snp_tumor_data_path} {snp_origin_data}", shell=True)
    #     sp.call(f"bcftools isec -p {indel_output_dir}/ {indel_teratoma_data_path} {indel_origin_data}", shell=True)




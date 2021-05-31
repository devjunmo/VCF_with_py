import glob
import os
import subprocess as sp
import natsort



############################################################

# INPUT_DIR = r'/myData/WES/data/vcf/hard/WES1_210420/'
INPUT_DIR = r'/myData/WES/data/vcf/cnn/WES1_210420/'
OUTPUT_DIR = INPUT_DIR

# 나중에 다른 프로젝트 추가되면 리펙토링

# pipeline type : g_hard / g_cnn
pipeline_type = 'g_cnn'

file_type = "ALL" # ST / SI / IndT / IndI / ALL

# hard filter
ST = r'hard*SNP_Teratoma*.vcf.gz'
SI = r'hard*SNP_hiPS*.vcf.gz'
IndT = r'hard*INDEL_Teratoma*.vcf.gz'
IndI = r'hard*INDEL_hiPS*.vcf.gz' 

# cnn filter
ST_c = r'hard*SNP_Teratoma*.vcf.gz'
SI_c = r'hard*SNP_hiPS*.vcf.gz'
IndT_c = r'hard*INDEL_Teratoma*.vcf.gz'
IndI_c = r'hard*INDEL_hiPS*.vcf.gz' 





############################################################


def mk_vcf_file_list(_input_dir, _input_form, group_name):
    _input_path_list = glob.glob(_input_dir + _input_form)
    _input_path_list = natsort.natsorted(_input_path_list)
    _output_path = _input_dir + group_name + '.list'
    
    f = open(f'{_output_path}', mode='w')
    for i in range(len(_input_path_list)):
        data = f'{_input_path_list[i]}\n'
        f.write(data)
    f.close



if pipeline_type == 'g_hard':
    if file_type == "ST":
        mk_vcf_file_list(INPUT_DIR, ST, "SNP_Teratoma")
    elif file_type == "SI":
        mk_vcf_file_list(INPUT_DIR, SI, "SNP_IPS")
    elif file_type == "IndT":
        mk_vcf_file_list(INPUT_DIR, IndT, "INDEL_Teratoma")
    elif file_type == "IndI":
        mk_vcf_file_list(INPUT_DIR, IndI, "INDEL_IPS")
    elif file_type == "ALL":
        mk_vcf_file_list(INPUT_DIR, ST, "SNP_Teratoma")
        mk_vcf_file_list(INPUT_DIR, SI, "SNP_IPS")
        mk_vcf_file_list(INPUT_DIR, IndT, "INDEL_Teratoma")
        mk_vcf_file_list(INPUT_DIR, IndI, "INDEL_IPS")


elif pipeline_type == 'g_cnn':
    if file_type == "ST":
        mk_vcf_file_list(INPUT_DIR, ST_c, "SNP_Teratoma_cnn")
    elif file_type == "SI":
        mk_vcf_file_list(INPUT_DIR, SI_c, "SNP_IPS_cnn")
    elif file_type == "IndT":
        mk_vcf_file_list(INPUT_DIR, IndT_c, "INDEL_Teratoma_cnn")
    elif file_type == "IndI":
        mk_vcf_file_list(INPUT_DIR, IndI_c, "INDEL_IPS_cnn")
    elif file_type == "ALL":
        mk_vcf_file_list(INPUT_DIR, ST_c, "SNP_Teratoma_cnn")
        mk_vcf_file_list(INPUT_DIR, SI_c, "SNP_IPS_cnn")
        mk_vcf_file_list(INPUT_DIR, IndT_c, "INDEL_Teratoma_cnn")
        mk_vcf_file_list(INPUT_DIR, IndI_c, "INDEL_IPS_cnn")








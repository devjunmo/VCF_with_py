import glob
import os
import subprocess as sp
import natsort
import pandas as pd
import numpy as np
from enum import Enum


# 디렉토리 순회해서 sed명령 줘서 txt파일 만들기
# txt파일들 하나씩 읽어와서 df화 시키기
# 그 후 새로운 파이썬 파일 만들어서 밴다이어그램 그리기
# 벤다이어그램 그리고 나서는 vcfpy를 사용하여 dp정보, mq정보 뽑아보기


# target_snp_root_dir = r'/myData/WES/data/vcf/hard/WES1_210420/snp_isec/' # 순회할 루트 디렉토리
# target_indel_root_dir = r'/myData/WES/data/vcf/hard/WES1_210420/indel_isec/' # 순회할 루트 디렉토리

# target_snp_root_dir = r'/myData/WES/data/vcf/raw/snp_isec/' # 순회할 루트 디렉토리
# target_indel_root_dir = r'/myData/WES/data/vcf/raw/indel_isec/' # 순회할 루트 디렉토리

# target_snp_root_dir = r'/myData/WES/data/vcf/cnn/WES1_210420/snp_isec/' # 순회할 루트 디렉토리
# target_indel_root_dir = r'/myData/WES/data/vcf/cnn/WES1_210420/indel_isec/' # 순회할 루트 디렉토리

target_snp_root_dir = r'/myData/WES/data/vcf/re_hard/snp_isec/' # 순회할 루트 디렉토리
target_indel_root_dir = r'/myData/WES/data/vcf/re_hard/indel_isec/' # 순회할 루트 디렉토리

# run_mode = 'count' # 'count', 'df'
run_mode = 'df'


# save_path = r'/myData/WES/data/vcf/hard/WES1_210420/'
# save_path = r'/myData/WES/data/vcf/raw/'
save_path = r'/myData/WES/data/vcf/re_hard/'


def mk_vcf_pos_count_files_in_dir(root_dir, prefix):
    files = os.listdir(root_dir)
    for file in files:
        path = os.path.join(root_dir, file)
        if os.path.isdir(path):
            mk_vcf_pos_count_files_in_dir(path, prefix) # 디렉토리면 재귀호출
        else:
            flag = path.split('.')[-1]
            if flag == "vcf":
                output_name = path.split('/')[-2] # 파일이면 작업수행 1. 이름따기, 2. sed먹이기 3. txt파일화 하기 
                output_dir = root_dir + '/'
                output_path = output_dir + output_name + '.count'
                sp.call(f"sed '/#/'d {path} | wc -l >> {output_path}", shell=True)


# 다시 순회하면서 기록한 텍스트 파일을 읽어서 리스트로 가져온 다음 데이터프레임에 붙인다
# index: sample / col: 차1 차2 교1 교2

def extract_count_from_file(_file_path):
    # 한줄씩 읽어서 input_path_list에 넣기 
    _count_lst = []
    f = open(_file_path, 'r')
    for i in f.readlines():
        _count_lst.append(i[:-1])
    return _count_lst

def count_file_to_df(root_dir, prefix, empty_df):
    files = os.listdir(root_dir)
    for file in files:
        path = os.path.join(root_dir, file)
        if os.path.isdir(path):
            count_file_to_df(path, prefix, empty_df) # 디렉토리면 재귀호출
        else:
            flag = path.split('.')[-1]
            if flag == "count":
                index_name = path.split('.')[-2].split('/')[-1]
                count_lst = extract_count_from_file(path)
                empty_df.loc[index_name] = count_lst


if run_mode == 'count':
    # count 파일 만들기
    mk_vcf_pos_count_files_in_dir(target_snp_root_dir, "")
    mk_vcf_pos_count_files_in_dir(target_indel_root_dir, "")

elif run_mode == 'df':
    # count파일로 df만들기
    snp_count_df = pd.DataFrame(columns=['snp_p_T', 'snp_p_I', 'snp_s_T', 'snp_s_I'])
    indel_count_df = pd.DataFrame(columns=['indel_p_T', 'indel_p_I', 'indel_s_T', 'indel_s_I'])

    count_file_to_df(target_snp_root_dir, "", snp_count_df)
    count_file_to_df(target_indel_root_dir, "", indel_count_df)

    print(snp_count_df)
    print(indel_count_df)

    snp_count_df.to_csv(save_path + 'WES_Teratoma_ips_snp_count.csv', header=True, index=True)
    indel_count_df.to_csv(save_path + 'WES_Teratoma_ips_indel_count.csv', header=True, index=True)






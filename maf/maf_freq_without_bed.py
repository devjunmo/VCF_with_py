from numpy.core.fromnumeric import shape, var
import pandas as pd
import numpy as np
from glob import glob
import subprocess as sp
import maf_functions 
import natsort
import os

# bed파일을 maf에 붙여주는 이유는 GT 처리때문임 

# SNP기준 코드

# maf파일에 SNP는 Start_Position = End_Position이고,
# maf's Start_Position = maf's End_Position = bed's End_Position

# 루프 첫머리에 sample name이 서로 맞는지 확인해주는 과정이 필요함.

maf_WGS_inputs_dir = r'/data_244/VCF/1000G_VCF_sample/WGS/maf/rm_hd_maf/'
maf_WGS_input_format = r'*.maf'

bed_WGS_inputs_dir = r'/data_244/VCF/gatherd_noDP_WGS_WES_interval_apply/WGS_specific/processed_bed/'
bed_WGS_input_format = r'SNP*.bed'

maf_WES_inputs_dir = r'/data_244/VCF/gatherd_WGS_WES_interval_apply/WES_specific/maf/rm_hd_maf/'
maf_WES_input_format = r'rmHd_SNP_*.maf'

bed_WES_inputs_dir = r'/data_244/VCF/gatherd_WGS_WES_interval_apply/WES_specific/processed_bed/'
bed_WES_input_format = r'SNP*.bed'

maf_Tsp_inputs_dir = r'/data_244/WES/Tsp_annotation_210608/DP_filter_apply/maf/rmHd_maf/'
maf_Tsp_input_format = 'rmHd_SNP_*.maf'

bed_Tsp_input_dir = r'/data_244/WES/Tsp_annotation_210608/DP_filter_apply/teratoma_specific_processed_bed_210609/'
bed_Tsp_input_format = r'SNP*.bed'


bed_header = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']

key = ['Chromosome', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']


var_classification_ref = r'/data_244/WES/maf_variant_classification.csv'

freq_out_root_dir = r'/data_244/VCF/no_dp_WGS_WES_var_class_interval_apply/'
wgs_specific_freq_dir = r'WGS_freq/'
wes_specific_freq_dir = r'WES_freq/'


# freq_output_dir = r'/data_244/VCF/WGS_WES_var_class/'
# freq_under_30_name = 'freq_under30.csv'
# freq_over_30_name = 'freq_over30.csv'
# freq_tonly_name = 'tonly_over30.csv'
# freq_not_tonly_name = 'not_tonly_over30.csv'


# tonly는 영상에서 가져온 데이터라서
t_only_bed_dir = r'/home/jun9485/data/WES/T-only_bed/'
t_only_format = r'SNP*.bed'
t_only_header = ['Chromosome', 'Start_Position', 'End_Position']
t_only_processing_key = ['Chromosome', 'End_Position']


maf_WGS_inputs = maf_WGS_inputs_dir + maf_WGS_input_format
bed_WGS_processed = bed_WGS_inputs_dir + bed_WGS_input_format

maf_WES_DP_inputs = maf_WES_inputs_dir + maf_WES_input_format
bed_WES_DP_processed = bed_WES_inputs_dir + bed_WES_input_format

maf_Tsp_inputs = maf_Tsp_inputs_dir + maf_Tsp_input_format
bed_Tsp_processed = bed_Tsp_input_dir + bed_Tsp_input_format



maf_WGS_input_lst = glob(maf_WGS_inputs)
maf_WGS_input_lst = natsort.natsorted(maf_WGS_input_lst)

maf_WES_input_lst = glob(maf_WES_DP_inputs)
maf_WES_input_lst = natsort.natsorted(maf_WES_input_lst)

maf_Tsp_input_lst = glob(maf_Tsp_inputs)
maf_Tsp_input_lst = natsort.natsorted(maf_Tsp_input_lst)


bed_WGS_input_lst = glob(bed_WGS_processed)
bed_WGS_input_lst = natsort.natsorted(bed_WGS_input_lst)
bed_WES_input_lst = glob(bed_WES_DP_processed)
bed_WES_input_lst = natsort.natsorted(bed_WES_input_lst)

bed_Tsp_input_lst = glob(bed_Tsp_processed)
bed_Tsp_input_lst = natsort.natsorted(bed_Tsp_input_lst)


if os.path.isdir(freq_out_root_dir) is False:
    os.mkdir(freq_out_root_dir)

if os.path.isdir(freq_out_root_dir + wgs_specific_freq_dir) is False:
    os.mkdir(freq_out_root_dir + wgs_specific_freq_dir)

if os.path.isdir(freq_out_root_dir + wes_specific_freq_dir) is False:
    os.mkdir(freq_out_root_dir + wes_specific_freq_dir)


# print(maf_WGS_input_lst) # '/data_244/VCF/gatherd_WGS_WES/WGS_specific/maf/rm_hd_maf/rmHd_SNP_Teratoma-7.maf'   '/data_244/VCF/gatherd_WGS_WES/WGS_specific/maf/rm_hd_maf/rmHd_SNP_hiPS66-B.maf'
# print(maf_WES_input_lst) # '/data_244/VCF/gatherd_WGS_WES/WES_specific/maf/rm_hd_maf/rmHd_SNP_Teratoma-7.maf'   '/data_244/VCF/gatherd_WGS_WES/WES_specific/maf/rm_hd_maf/rmHd_SNP_hiPS66-B.maf'
# print(maf_Tsp_input_lst) # '/data_244/WES/Tsp_annotation_210608/DP_filter_apply/maf/rmHd_maf/rmHd_SNP_Teratoma-6.maf'   '/data_244/WES/Tsp_annotation_210608/DP_filter_apply/maf/rmHd_maf/rmHd_SNP_Teratoma-26.maf'


# print(bed_WGS_input_lst) # '/data_244/VCF/gatherd_WGS_WES/WGS_specific/processed_bed/SNP_Teratoma-7.bed'    '/data_244/VCF/gatherd_WGS_WES/WGS_specific/processed_bed/SNP_hiPS66-B.bed'
# print(bed_WES_input_lst) # '/data_244/VCF/gatherd_WGS_WES/WES_specific/processed_bed/SNP_Teratoma-7.bed'    '/data_244/VCF/gatherd_WGS_WES/WES_specific/processed_bed/SNP_hiPS66-B.bed'
# print(bed_Tsp_input_lst) # '/data_244/WES/Tsp_annotation_210608/DP_filter_apply/teratoma_specific_processed_bed_210609/SNP_Teratoma-6.bed'  '/data_244/WES/Tsp_annotation_210608/DP_filter_apply/teratoma_specific_processed_bed_210609/SNP_Teratoma-26.bed'


def get_common_sample_lst(_input_lst1, _input_lst2):
    _input_lst1 = [_input.split(r'/')[-1] for _input in _input_lst1]
    _input_lst2 = [_input.split(r'/')[-1] for _input in _input_lst2]
    _input_lst1, _input_lst_2 = set(_input_lst1), set(_input_lst2)
    _commons = list(_input_lst1 & _input_lst_2)
    return _commons


common_tsp_WGS_maf_names = get_common_sample_lst(maf_WGS_input_lst, maf_Tsp_input_lst)
common_tsp_WGS_bed_names = get_common_sample_lst(bed_WGS_input_lst, bed_Tsp_input_lst)

common_tsp_WGS_maf_names = natsort.natsorted(common_tsp_WGS_maf_names)
common_tsp_WGS_bed_names = natsort.natsorted(common_tsp_WGS_bed_names)

# print(common_tsp_WGS_maf_names)
# print(common_tsp_WGS_bed_names)


common_tsp_maf_path_lst = []
common_WGS_maf_path_lst = []

common_tsp_bed_path_lst = []
common_WGS_bed_path_lst = []


for maf_name in common_tsp_WGS_maf_names:
    common_tsp_maf_path_lst.append(maf_Tsp_inputs_dir + maf_name)
    common_WGS_maf_path_lst.append(maf_WGS_inputs_dir + maf_name)

for bed_name in common_tsp_WGS_bed_names:
    common_tsp_bed_path_lst.append(bed_Tsp_input_dir + bed_name)
    common_WGS_bed_path_lst.append(bed_WGS_inputs_dir + bed_name)


# dp_under_30_variant_df = pd.DataFrame()
# dp_more_30_variant_df = pd.DataFrame()
# t_only_variant_df = pd.DataFrame()
# not_tonly_variant_df = pd.DataFrame()


# var_class_ref_df = pd.read_csv(var_classification_ref, index_col='var_class')

# ref_grp_dict = maf_functions.var_ref_df_to_dict(var_class_ref_df, 'var_grp')

# processed_tonly_bed_data_list = []

# print(ref_dict)


# exit(0)



def mk_grp_freq_table(_freq_table, _ref_grp_dict, _save_variant_df):
    
    var_classes = list(_freq_table.index)
    var_freq_counts = list(_freq_table.values)

    for j in range(len(var_classes)):
        variant_key = var_classes[j]
        variant_count = var_freq_counts[j]

        _save_variant_df.fillna(0, inplace=True)

        try:
            _save_variant_df.loc[sample_name, _ref_grp_dict[variant_key]] = _save_variant_df.loc[sample_name, _ref_grp_dict[variant_key]] + variant_count
        
        except KeyError:
            _save_variant_df.loc[sample_name, _ref_grp_dict[variant_key]] = 0
            _save_variant_df.loc[sample_name, _ref_grp_dict[variant_key]] = _save_variant_df.loc[sample_name, _ref_grp_dict[variant_key]] + variant_count


# merge할때 dtype이 동일해야 함. 특히 Tonly는 sex크로모좀이 없을때도 있어서 chr dtype이 int로 될때가 있음. 
def casting_key_to_str(_df_lst, _key):
    for _df in _df_lst:
        for _col in _key:
            _df[_col] = _df[_col].astype(str)


def gt_processing(_bed_df, _maf_df, _key, _how='left', _inc_nan_col='Variant_Classification'):

    casting_key_to_str([_bed_df, _maf_df], _key)

    # print('in the gt processign function')
    # print(wgs_maf_df.loc[:, ["Chromosome", "Start_Position", "End_Position", "Variant_Classification"]])
    # print(_bed_df)
    # print(_maf_df)
    # print(_bed_df.loc[:, ["Chromosome", "Start_Position", "End_Position"]])
    # print(_maf_df.loc[:, ["Chromosome", "Start_Position", "End_Position", "Variant_Classification"]])
    # print('func out!')

    _joined_df = pd.merge(_bed_df, _maf_df, how=_how, on = _key)

    _processed_maf_df = _joined_df[-_joined_df[_inc_nan_col].isnull()] # bed파일에 붙였고, bed에 매치되는게 없는 maf포지션에 대해서는 NaN이 박힘
    
    return _processed_maf_df



# exit(0)

# WES specific data analysis

# for i in range(len(maf_WES_input_lst)):
#     sample_type = maf_WES_input_lst[i].split(r'/')[-1].split('.')[0].split('_')[1]
#     sample_name = maf_WES_input_lst[i].split(r'/')[-1].split('.')[0].split('_')[2]
#     sample_name = sample_type + '_' + sample_name

#     print('[ WES ]', sample_name)

#     # maf data 불러오기
    
#     wes_maf_df = pd.read_csv(maf_WES_input_lst[i], sep='\t', low_memory=False)
#     # print(wes_maf_df)

#     # bed data 불러오기

#     wes_bed_df = pd.read_csv(bed_WES_input_lst[i], sep='\t', names=bed_header, low_memory=False)
#     # print(wes_bed_df)
#     # print(wes_bed_df.dtypes)

    
#     # GT 제거 전처리

#     processed_WES_maf = gt_processing(wes_bed_df, wes_maf_df, key)
#     # print(processed_WES_maf)
#     # print(wes_bed_df.dtypes)

#     # variant check

#     wes_variant_freq = pd.value_counts(processed_WES_maf['Variant_Classification'])

#     print(wes_variant_freq)

#     wes_variant_freq.to_csv(freq_out_root_dir + wes_specific_freq_dir + sample_name + '.csv')
#     # break

# exit(0)



# WGS specific data analysis

for i in range(len(maf_WGS_input_lst)):
    sample_type = maf_WGS_input_lst[i].split(r'/')[-1].split('.')[0].split('_')[1]
    sample_name = maf_WGS_input_lst[i].split(r'/')[-1].split('.')[0].split('_')[2]
    sample_name = sample_type + '_' + sample_name

    print('[ WGS ]', sample_name)

    # maf data 불러오기
    
    wgs_maf_df = pd.read_csv(maf_WGS_input_lst[i], sep='\t', low_memory=False)
    # print(wgs_maf_df.loc[:, ["Chromosome", "Start_Position", "End_Position", "Variant_Classification"]])

    # bed data 불러오기

    # wgs_bed_df = pd.read_csv(bed_WGS_input_lst[i], sep='\t', names=bed_header, low_memory=False, \
    #             dtype={'Chromosome':object, 'Start_Position':object, 'End_Position':object})
    # print(wgs_bed_df)
    # print(wgs_bed_df.dtypes)
    

    # GT 제거 전처리

    # processed_WGS_maf = gt_processing(wgs_bed_df, wgs_maf_df, key)
    # print(wgs_bed_df)
    # print(wgs_bed_df.dtypes)
    # print(processed_WGS_maf.loc[:, ["Chromosome", "Start_Position", "End_Position"]])

    # variant check

    wgs_variant_freq = pd.value_counts(wgs_maf_df['Variant_Classification'])

    print(wgs_variant_freq)

    # wgs_variant_freq.to_csv(freq_out_root_dir + wgs_specific_freq_dir + sample_name + '.csv')

    # break

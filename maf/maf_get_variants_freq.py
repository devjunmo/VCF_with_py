from numpy.core.fromnumeric import shape, var
import pandas as pd
import numpy as np
from glob import glob
import subprocess as sp
import maf_functions 

# bed파일을 maf에 붙여주는 이유는 GT 처리때문임 

# SNP기준 코드
# maf파일에 SNP는 Start_Position = End_Position이고,
# maf's Start_Position = maf's End_Position = bed's End_Position

# 루프 첫머리에 sample name이 서로 맞는지 확인해주는 과정이 필요함.

maf_Tsp_NoDP_inputs = r'/home/jun9485/data/WES/Tsp_annotation_210608/no_DP_filter/maf/rmHd_maf/rmHd_SNP_*.maf'
bed_Tsp_NoDP_processed = r'/home/jun9485/data/WES/Tsp_annotation_210608/no_DP_filter/teratoma_specific_processed_bed_210609/SNP*.bed'

maf_Tsp_T_DP_inputs = r'/home/jun9485/data/WES/Tsp_annotation_210608/DP_filter_apply/maf/rmHd_maf/rmHd_SNP_*.maf'
bed_Tsp_T_DP_processed = r'/home/jun9485/data/WES/Tsp_annotation_210608/DP_filter_apply/teratoma_specific_processed_bed_210609/SNP*.bed'

bed_header = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']
key = ['Chromosome', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']

var_classification_ref = r'/home/jun9485/data/WES/maf_variant_classification.csv'

freq_output_dir = r'/home/jun9485/data/WES/Tsp_annotation_210608/maf_analysis/'
freq_under_30_name = 'freq_under30.csv'
freq_over_30_name = 'freq_over30.csv'
freq_tonly_name = 'tonly_over30.csv'
freq_not_tonly_name = 'not_tonly_over30.csv'

t_only_bed_dir = r'/home/jun9485/data/WES/T-only_bed/'
t_only_format = r'SNP*.bed'
t_only_header = ['Chromosome', 'Start_Position', 'End_Position']
t_only_processing_key = ['Chromosome', 'End_Position']


maf_Tsp_NoDP_input_lst = glob(maf_Tsp_NoDP_inputs)
maf_Tsp_T_DP_input_lst = glob(maf_Tsp_T_DP_inputs)

bed_Tsp_NoDP_input_lst = glob(bed_Tsp_NoDP_processed)
bed_Tsp_T_DP_input_lst = glob(bed_Tsp_T_DP_processed)
bed_Tonly_dpApply_input_lst = glob(t_only_bed_dir + t_only_format)

# print(bed_Tsp_NoDP_input_lst[7])
# print(bed_Tsp_T_DP_input_lst[7])
# print(bed_Tonly_dpApply_input_lst[7])

# exit()


dp_under_30_variant_df = pd.DataFrame()
dp_more_30_variant_df = pd.DataFrame()
t_only_variant_df = pd.DataFrame()
not_tonly_variant_df = pd.DataFrame()


var_class_ref_df = pd.read_csv(var_classification_ref, index_col='var_class')

ref_grp_dict = maf_functions.var_ref_df_to_dict(var_class_ref_df, 'var_grp')

processed_tonly_bed_data_list = []

# print(ref_dict)


# exit(0)


def sample_match_test():
    pass


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



for i in range(len(maf_Tsp_NoDP_input_lst)):
    sample_type = maf_Tsp_NoDP_input_lst[i].split(r'/')[-1].split('.')[0].split('_')[1]
    sample_name = maf_Tsp_NoDP_input_lst[i].split(r'/')[-1].split('.')[0].split('_')[2]
    sample_name = sample_type + '_' + sample_name

    print(sample_name)

    # maf data 불러오기
    
    noDP_maf_df = pd.read_csv(maf_Tsp_NoDP_input_lst[i], sep='\t', low_memory=False)
    dP_maf_df = pd.read_csv(maf_Tsp_T_DP_input_lst[i], sep='\t', low_memory=False)

    # bed data 불러오기

    noDP_bed_df = pd.read_csv(bed_Tsp_NoDP_input_lst[i], sep='\t', names=bed_header, low_memory=False)
    dP_bed_df = pd.read_csv(bed_Tsp_T_DP_input_lst[i], sep='\t', names=bed_header, low_memory=False)
    t_only_bed_df = pd.read_csv(bed_Tonly_dpApply_input_lst[i], sep='\t', names=t_only_header, low_memory=False)

    # print(dP_bed_df)
    # print(t_only_bed_df)


    # Tonly bedfile과 Tsp30 bed file을 left join하여 Tonly에 REF ALT 달아주기

    # t_only_bed_df['Chromosome'] = t_only_bed_df['Chromosome'].astype('str')
    # t_only_bed_df['Start_Position'] = t_only_bed_df['Start_Position'].astype('str')
    # t_only_bed_df['End_Position'] = t_only_bed_df['End_Position'].astype('str')
    
    casting_key_to_str([noDP_maf_df, dP_maf_df, noDP_bed_df, dP_bed_df], key)
    casting_key_to_str([t_only_bed_df], t_only_header)


    t_only_bed_df_processed = pd.merge(t_only_bed_df, dP_bed_df, how='left', on = t_only_processing_key)
    # print(t_only_bed_df_processed)
    # print(type(t_only_bed_df_processed))

    # exit(0)

    
    # GT 제거 전처리

    no_dp_join_df = pd.merge(noDP_bed_df, noDP_maf_df, how='left', \
                            on = key)

    dp_join_df = pd.merge(dP_bed_df, dP_maf_df, how='left', \
                            on = key)

    processed_maf_df_no_dp = no_dp_join_df[-no_dp_join_df['Variant_Classification'].isnull()] # bed파일에 붙였고, bed에 매치되는게 없는 maf포지션에 대해서는 NaN이 박힘
    processed_maf_df_apply_dp = dp_join_df[-dp_join_df['Variant_Classification'].isnull()]


    # (over / under) DP 30 maf file join

    joined_maf = pd.merge(processed_maf_df_no_dp, processed_maf_df_apply_dp, how='left', \
                            on = key)


    # under 30 변이 테이블 생성

    dp_under_30_mask = joined_maf['Variant_Classification_y'].isnull()
    dp_under_30_maf_df =  joined_maf[dp_under_30_mask]
    dp_under_30_variant_freq = pd.value_counts(dp_under_30_maf_df['Variant_Classification_x'])
    # print(dp_under_30_variant_freq)



    # over 30 변이 테이블 생성
    dp_more_30_maf_df = joined_maf[-dp_under_30_mask]
    dp_more_30_variant_freq = pd.value_counts(dp_more_30_maf_df['Variant_Classification_x'])
    # print(dp_more_30_variant_freq)

    # exit()



    # T-only 변이 테이블 생성 -> silent가 많다

    # get_processed_t_only_bed(t_only_bed_df, dP_bed_df, t_only_header)

    joined_maf_t_only = pd.merge(t_only_bed_df_processed, \
                                 processed_maf_df_apply_dp, how='left', on = key)

    # print(joined_maf_t_only)

    t_only_variant_freq = pd.value_counts(joined_maf_t_only['Variant_Classification'])
    # print(t_only_variant_freq)

    # exit()


    # T-only가 아닌 부분에 대한 변이 테이블 생성

    # print(processed_maf_df_apply_dp)

    t_only_bed_df_processed['flag_col'] = 'jun'


    joined_maf_not_t_only = pd.merge(t_only_bed_df_processed, \
                                    processed_maf_df_apply_dp, how='right', on = key)

    # print(list(joined_maf_not_t_only.columns))
    # print(joined_maf_not_t_only['flag_col'])
    # print(pd.value_counts(joined_maf_not_t_only['flag_col']))

    not_tonly_idx = joined_maf_not_t_only['flag_col'].isna()
    not_tonly_maf_df = joined_maf_not_t_only[not_tonly_idx]
    # print(not_tonly_maf_df.shape)
    not_t_only_variant_freq = pd.value_counts(not_tonly_maf_df['Variant_Classification'])
    print(not_t_only_variant_freq)

    # exit()



    mk_grp_freq_table(dp_under_30_variant_freq, ref_grp_dict, dp_under_30_variant_df)
    mk_grp_freq_table(dp_more_30_variant_freq, ref_grp_dict, dp_more_30_variant_df)
    mk_grp_freq_table(t_only_variant_freq, ref_grp_dict, t_only_variant_df)
    mk_grp_freq_table(not_t_only_variant_freq, ref_grp_dict, not_tonly_variant_df)


# exit()

# print(dp_under_30_variant_df)
# print(dp_under_30_variant_df.shape)
dp_under_30_variant_df.to_csv(freq_output_dir + freq_under_30_name)
dp_more_30_variant_df.to_csv(freq_output_dir + freq_over_30_name)
t_only_variant_df.to_csv(freq_output_dir + freq_tonly_name)
not_tonly_variant_df.to_csv(freq_output_dir + freq_not_tonly_name)

    # break


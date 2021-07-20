import pandas as pd
import numpy as np
from glob import glob
import subprocess as sp
import maf_functions 

# bed파일을 maf에 붙여주는 이유는 GT 처리때문임 

# SNP기준 코드
# maf파일에 SNP는 Start_Position = End_Position이고,
# maf's Start_Position = maf's End_Position = bed's End_Position

maf_Tsp_NoDP_inputs = r'/home/jun9485/data/WES/Tsp_annotation_210608/no_DP_filter/maf/rmHd_maf/rmHd_SNP_*.maf'
bed_Tsp_NoDP_processed = r'/home/jun9485/data/WES/Tsp_annotation_210608/no_DP_filter/teratoma_specific_processed_bed_210609/SNP*.bed'

maf_Tsp_T_DP_inputs = r'/home/jun9485/data/WES/Tsp_annotation_210608/DP_filter_apply/maf/rmHd_maf/rmHd_SNP_*.maf'
bed_Tsp_T_DP_processed = r'/home/jun9485/data/WES/Tsp_annotation_210608/DP_filter_apply/teratoma_specific_processed_bed_210609/SNP*.bed'

bed_header = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']
key = ['Chromosome', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']

var_classification_ref = r'/home/jun9485/data/WES/maf_variant_classification.csv'


maf_Tsp_NoDP_input_lst = glob(maf_Tsp_NoDP_inputs)
maf_Tsp_T_DP_input_lst = glob(maf_Tsp_T_DP_inputs)

bed_Tsp_NoDP_input_lst = glob(bed_Tsp_NoDP_processed)
bed_Tsp_T_DP_input_lst = glob(bed_Tsp_T_DP_processed)

var_class_ref_df = pd.read_csv(var_classification_ref, index_col='var_class')

ref_dict = maf_functions.var_ref_df_to_dict(var_class_ref_df, 'var_grp')

print(ref_dict)

print(ref_dict["Silent"])
# exit(0)

# def flip_TF_series(_tmp_lst):
#     for i in range(len(_tmp_lst)):
#         if tmp_lst[i] == True:
#             tmp_lst[i] = False
#         elif tmp_lst[i] == False:
#             tmp_lst[i] = True
    
#     return_series = pd.Series(_tmp_lst)
#     return return_series



for i in range(len(maf_Tsp_NoDP_input_lst)):
    # print(maf_Tsp_NoDP_input_lst[i].split(r'/')[-1].split('.')[0])
    # print(maf_Tsp_T_DP_input_lst[i].split(r'/')[-1].split('.')[0])
    # print(bed_Tsp_NoDP_input_lst[i].split(r'/')[-1].split('.')[0])
    # print(bed_Tsp_T_DP_input_lst[i].split(r'/')[-1].split('.')[0])
    # break
    print(maf_Tsp_NoDP_input_lst[i].split(r'/')[-1].split('.')[0])
    
    noDP_maf_df = pd.read_csv(maf_Tsp_NoDP_input_lst[i], sep='\t', low_memory=False)
    dP_maf_df = pd.read_csv(maf_Tsp_T_DP_input_lst[i], sep='\t', low_memory=False)

    noDP_bed_df = pd.read_csv(bed_Tsp_NoDP_input_lst[i], sep='\t', names=bed_header, low_memory=False)
    dP_bed_df = pd.read_csv(bed_Tsp_T_DP_input_lst[i], sep='\t', names=bed_header, low_memory=False)

    # print(noDP_bed_df['End_Position'].head())
    # print(noDP_maf_df['End_Position'].head())
    # print(noDP_bed_df['Start_Position'].head())
    # print(dP_bed_df['Start_Position'].head())

    

    # break


    no_dp_join_df = pd.merge(noDP_bed_df, noDP_maf_df, how='left', \
                            on = key)

    dp_join_df = pd.merge(dP_bed_df, dP_maf_df, how='left', \
                            on = key)

    # print(noDP_bed_df.shape)
    # print(no_dp_join_df.shape)
    # print(no_dp_join_df.isnull().values.any()) # True. True = NaN이 섞여있다. 
    #                                            # 여기서 NaN 포함 행은 processed specific이니까 기존 bed에 없던, 즉 bad GT 포지션임.

    # print(type(no_dp_join_df['Variant_Classification'].isnull())) # <class 'pandas.core.series.Series'>

    # nan_df_no_dp = no_dp_join_df[no_dp_join_df['Variant_Classification'].isnull()]
    processed_maf_df_no_dp = no_dp_join_df[-no_dp_join_df['Variant_Classification'].isnull()] # bed파일에 붙였고, bed에 매치되는게 없는 maf포지션에 대해서는 NaN이 박힘
    processed_maf_df_apply_dp = dp_join_df[-dp_join_df['Variant_Classification'].isnull()]
    # print(nan_df_no_dp['Variant_Classification'])
    # print(nan_df_no_dp.shape[0]) # 걸러진 GT의 갯수 

    # print('DP 필터 적용 전, GT 전처리 전:', processed_maf_df_no_dp.shape[0])
    # print('DP 필터 적용 후, GT 전처리 전:', dP_maf_df.shape[0])
    # print('DP 필터 적용 후, GT 전처리 후:', processed_maf_df_apply_dp.shape[0])
    # print('\n')

    joined_maf = pd.merge(processed_maf_df_no_dp, processed_maf_df_apply_dp, how='left', \
                            on = key)

    # # print(joined_maf['Variant_Classification_y'])
    # # print(joined_maf['Variant_Classification_y'].isnull())
    # tmp_lst = list(joined_maf['Variant_Classification_y'].isnull())
    # srs = flip_TF_series(tmp_lst) # y부분에 NaN이 아닌부분을 True로 놓아서 마이너스 기호로 빼줄것. 즉, DP적용이 못붙은 포지션만 남기려는 목적. < 생각해보니 굳이..? 그냥 마스크를 []안에 넣으면 될거같은데

    # # print(srs)

    # # print(joined_maf.shape)

    # DP_under_30_maf_df = joined_maf[-srs]

    # # print(interest_maf.shape)

    # # print(interest_maf['Variant_Classification_x'])

    # freq = pd.value_counts(interest_maf['Variant_Classification_x'])
    # print(freq)

    dp_under_30_mask = joined_maf['Variant_Classification_y'].isnull()

    dp_under_30_maf_df =  joined_maf[dp_under_30_mask]

    dp_under_30_variant_freq = pd.value_counts(dp_under_30_maf_df['Variant_Classification_x'])

    print(dp_under_30_variant_freq)

    print(type(dp_under_30_variant_freq)) # <class 'pandas.core.series.Series'>

    print(dp_under_30_variant_freq['Silent'])
    print(dp_under_30_variant_freq.index)
    print(dp_under_30_variant_freq.values)
    print(type(dp_under_30_variant_freq.index)) # <class 'pandas.core.indexes.base.Index'>
    print(type(dp_under_30_variant_freq.values)) # <class 'numpy.ndarray'>

    tmp_idx_lst = list(dp_under_30_variant_freq.index)
    tmp_val_lst = list(dp_under_30_variant_freq.values)

    print(type(tmp_idx_lst)) # list
    print(type(tmp_val_lst)) # list

    for j in range(len(tmp_idx_lst)):
        variant_key = tmp_idx_lst[j]
        variant_grp = tmp_val_lst[j]

        






    break
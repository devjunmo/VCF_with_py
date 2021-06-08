import pandas as pd
import numpy as np
from glob import glob
import subprocess as sp

# bed파일을 maf에 붙여주는 이유는 GT 처리때문임 

# SNP기준 코드
# maf파일에 SNP는 Start_Position = End_Position이고,
# maf's Start_Position = maf's End_Position = bed's End_Position

maf_Tsp_NoDP_inputs = r'/data_244/WES/no_DP_filter_samples/pass_only_Teratoma_specifics/maf/rm_top_hash_maf/rmHd_SNP_*.maf'
bed_Tsp_NoDP_processed = r'/data_244/WES/no_DP_filter_samples/pass_only_Teratoma_specifics/teratoma_specific_processed_bed_210606/SNP*.bed'

maf_Tsp_T_DP_inputs = r'/data_244/WES/T_DP_O_NDP_samples/pass_only_Teratoma_specifics/maf/rm_top_hash_maf/rmHd_SNP_*.maf'
bed_Tsp_T_DP_processed = r'/data_244/WES/T_DP_O_NDP_samples/pass_only_Teratoma_specifics/teratoma_specific_processed_bed_210606/SNP*.bed'

bed_header = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']
key = ['Chromosome', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2']


maf_Tsp_NoDP_input_lst = glob(maf_Tsp_NoDP_inputs)
maf_Tsp_T_DP_input_lst = glob(maf_Tsp_T_DP_inputs)

bed_Tsp_NoDP_input_lst = glob(bed_Tsp_NoDP_processed)
bed_Tsp_T_DP_input_lst = glob(bed_Tsp_T_DP_processed)

for i in range(len(maf_Tsp_NoDP_input_lst)):
    print(maf_Tsp_NoDP_input_lst[i].split(r'/')[-1].split('.')[0])
    
    noDP_maf_df = pd.read_csv(maf_Tsp_NoDP_input_lst[i], sep='\t', low_memory=False)
    dP_maf_df = pd.read_csv(maf_Tsp_T_DP_input_lst[i], sep='\t', low_memory=False)

    noDP_bed_df = pd.read_csv(bed_Tsp_NoDP_input_lst[i], sep='\t', names=bed_header, low_memory=False)
    dP_bed_df = pd.read_csv(bed_Tsp_T_DP_input_lst[i], sep='\t', names=bed_header, low_memory=False)

    # print(noDP_bed_df['End_Position'].head())
    # print(noDP_maf_df['End_Position'].head())


    no_dp_join_df = pd.merge(noDP_bed_df, noDP_maf_df, how='left', \
                            on = key)

    dp_join_df = pd.merge(dP_bed_df, dP_maf_df, how='left', \
                            on = key)

    # print(noDP_bed_df.shape)
    # print(no_dp_join_df.shape)
    # print(no_dp_join_df.isnull().values.any()) # True. True = NaN이 섞여있다. 
    #                                            # 여기서 NaN 포함 행은 processed specific이니까 기존 bed에 없던, 즉 bad GT 포지션임.

    nan_df_no_dp = no_dp_join_df[no_dp_join_df['Variant_Classification'].isnull()]
    # print(nan_df_no_dp['Variant_Classification'])
    print(nan_df_no_dp.shape) # 걸러진 GT의 갯수


    


    # break
from numpy.core.fromnumeric import var
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

freq_output_dir = r'/home/jun9485/data/WES/Tsp_annotation_210608/maf_analysis/'
freq_under_30_name = 'freq_under30.csv'
freq_over_30_name = 'freq_over30.csv'


maf_Tsp_NoDP_input_lst = glob(maf_Tsp_NoDP_inputs)
maf_Tsp_T_DP_input_lst = glob(maf_Tsp_T_DP_inputs)

bed_Tsp_NoDP_input_lst = glob(bed_Tsp_NoDP_processed)
bed_Tsp_T_DP_input_lst = glob(bed_Tsp_T_DP_processed)


dp_under_30_variant_df = pd.DataFrame()
dp_more_30_variant_df = pd.DataFrame()


var_class_ref_df = pd.read_csv(var_classification_ref, index_col='var_class')

ref_dict = maf_functions.var_ref_df_to_dict(var_class_ref_df, 'var_grp')

# print(ref_dict)


# exit(0)



for i in range(len(maf_Tsp_NoDP_input_lst)):
    sample_type = maf_Tsp_NoDP_input_lst[i].split(r'/')[-1].split('.')[0].split('_')[1]
    sample_name = maf_Tsp_NoDP_input_lst[i].split(r'/')[-1].split('.')[0].split('_')[2]
    sample_name = sample_type + '_' + sample_name

    print(maf_Tsp_NoDP_input_lst[i].split(r'/')[-1].split('.')[0])
    
    noDP_maf_df = pd.read_csv(maf_Tsp_NoDP_input_lst[i], sep='\t', low_memory=False)
    dP_maf_df = pd.read_csv(maf_Tsp_T_DP_input_lst[i], sep='\t', low_memory=False)

    noDP_bed_df = pd.read_csv(bed_Tsp_NoDP_input_lst[i], sep='\t', names=bed_header, low_memory=False)
    dP_bed_df = pd.read_csv(bed_Tsp_T_DP_input_lst[i], sep='\t', names=bed_header, low_memory=False)

    no_dp_join_df = pd.merge(noDP_bed_df, noDP_maf_df, how='left', \
                            on = key)

    dp_join_df = pd.merge(dP_bed_df, dP_maf_df, how='left', \
                            on = key)

    processed_maf_df_no_dp = no_dp_join_df[-no_dp_join_df['Variant_Classification'].isnull()] # bed파일에 붙였고, bed에 매치되는게 없는 maf포지션에 대해서는 NaN이 박힘
    processed_maf_df_apply_dp = dp_join_df[-dp_join_df['Variant_Classification'].isnull()]


    joined_maf = pd.merge(processed_maf_df_no_dp, processed_maf_df_apply_dp, how='left', \
                            on = key)

    dp_under_30_mask = joined_maf['Variant_Classification_y'].isnull()

    dp_under_30_maf_df =  joined_maf[dp_under_30_mask]

    dp_under_30_variant_freq = pd.value_counts(dp_under_30_maf_df['Variant_Classification_x'])

    print(dp_under_30_variant_freq)

    under30_var_classes = list(dp_under_30_variant_freq.index)
    under30_var_freq_counts = list(dp_under_30_variant_freq.values)

    for j in range(len(under30_var_classes)):
        variant_key = under30_var_classes[j]
        variant_count = under30_var_freq_counts[j]

        dp_under_30_variant_df.fillna(0, inplace=True)

        try:
            dp_under_30_variant_df.loc[ref_dict[variant_key], sample_name] = dp_under_30_variant_df.loc[ref_dict[variant_key], sample_name] + variant_count
        except KeyError:
            dp_under_30_variant_df.loc[ref_dict[variant_key], sample_name] = 0
            dp_under_30_variant_df.loc[ref_dict[variant_key], sample_name] = dp_under_30_variant_df.loc[ref_dict[variant_key], sample_name] + variant_count
            


        # dp_under_30_variant_df.loc[ref_dict[variant_key], sample_name] = dp_under_30_variant_df.loc[ref_dict[variant_key], sample_name] + variant_count
    
    # print(dp_under_30_variant_df)

    

print(dp_under_30_variant_df)
dp_under_30_variant_df.to_csv(dp_under_30_variant_df+freq_under_30_name, sep=',')
    # break


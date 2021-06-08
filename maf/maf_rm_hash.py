from glob import glob
import os
import subprocess as sp

maf_Tsp_NoDP_input_dir = r'/home/jun9485/data/WES/no_DP_filter_samples/pass_only_Teratoma_specifics/maf/rm_top_hash_maf/'
bed_Tsp_NoDP_processed = r'/home/jun9485/data/WES/no_DP_filter_samples/pass_only_Teratoma_specifics/teratoma_specific_processed_bed_210606/*'

maf_Tsp_T_DP_input_dir = r'/home/jun9485/data/WES/T_DP_O_NDP_samples/pass_only_Teratoma_specifics/maf/rm_top_hash_maf/'
bed_Tsp_T_DP_processed = r'/home/jun9485/data/WES/T_DP_O_NDP_samples/pass_only_Teratoma_specifics/teratoma_specific_processed_bed_210606/*'


# NO DP

# maf_Tsp_NoDP_input_lst = glob(maf_Tsp_NoDP_input_dir + '*')

# os.chdir(maf_Tsp_NoDP_input_dir)

# for maf_Tsp_NoDP in maf_Tsp_NoDP_input_lst:

#     f_name = maf_Tsp_NoDP.split(r'/')[-1].split(r'.')[0]
#     print(f_name)
#     sp.call(rf'grep -v "#" {f_name}.maf > rmHd_{f_name}.maf', shell=True)



# DP

maf_Tsp_T_DP_input_lst = glob(maf_Tsp_T_DP_input_dir + '*')

os.chdir(maf_Tsp_T_DP_input_dir)

for maf_Tsp_T_DP in maf_Tsp_T_DP_input_lst:

    f_name = maf_Tsp_T_DP.split(r'/')[-1].split(r'.')[0]
    print(f_name)
    sp.call(rf'grep -v "#" {f_name}.maf > rmHd_{f_name}.maf', shell=True)

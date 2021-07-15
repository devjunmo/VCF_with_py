from glob import glob
import os
import subprocess as sp

input_maf_dir = r'E:/UTUC_data/WES/maf/mutect1/sample1/'
output_maf_dir_name = r'rm_hd_maf/'

# maf_Tsp_NoDP_input_dir = r'/home/jun9485/data/WES/Tsp_annotation_210608/no_DP_filter/maf/'

# maf_Tsp_T_DP_input_dir = r'/home/jun9485/data/WES/Tsp_annotation_210608/DP_filter_apply/maf/'


input_maf_lst = glob(input_maf_dir + '*.maf')

old = os.getcwd()

os.chdir(input_maf_dir)

output_dir = input_maf_dir + output_maf_dir_name

if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)


for maf_file in input_maf_lst:
    f_name = maf_file.split(r'/')[-1].split(r'.')[0]
    print(f_name)
    sp.call(rf'grep -v "#" {f_name}.maf > rmHd_{f_name}.maf', shell=True)


sp.call(rf'mv rmHd_* {output_maf_dir_name}', shell = True)

os.chdir(old)


# NO DP

# maf_Tsp_NoDP_input_lst = glob(maf_Tsp_NoDP_input_dir + '*.maf')

# os.chdir(maf_Tsp_NoDP_input_dir)

# for maf_Tsp_NoDP in maf_Tsp_NoDP_input_lst:

#     f_name = maf_Tsp_NoDP.split(r'/')[-1].split(r'.')[0]
#     print(f_name)
#     sp.call(rf'grep -v "#" {f_name}.maf > rmHd_{f_name}.maf', shell=True)



# # DP

# maf_Tsp_T_DP_input_lst = glob(maf_Tsp_T_DP_input_dir + '*.maf')

# os.chdir(maf_Tsp_T_DP_input_dir)

# for maf_Tsp_T_DP in maf_Tsp_T_DP_input_lst:

#     f_name = maf_Tsp_T_DP.split(r'/')[-1].split(r'.')[0]
#     print(f_name)
#     sp.call(rf'grep -v "#" {f_name}.maf > rmHd_{f_name}.maf', shell=True)

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from jun_tools import jun_mtd as jm  # pip install YjmTool
import os


# 디렉토리의 마프파일 가져오기
# {Tumor_id: [N_id, row_count]} 딕셔너리 만들기


# input_dir = r"E:/stemcell/somatic_analysis/maf/vardict/DP_AF_filtered_maf/exonic_maf/origin"
# input_format = r"*.maf"
# input_lst = jm.get_input_path_list(input_dir, input_format, False)

# output_dir_name = r"comp_result"
# output_dir = jm.set_output_dir(input_dir, output_dir_name)


# mutation_profile_dict = dict()

# mut_prof_df = pd.DataFrame(columns=['Origin', 'mut_count'])


# def mutation_profile():
#     pass


# def logical_profile():
#     pass


# for i in range(len(input_lst)):
#     input_maf_df = pd.read_csv(input_lst[i], sep='\t')
#     # print(input_maf_df)
#     row_count = input_maf_df.shape[0]
#     t_id = input_maf_df['Tumor_Sample_Barcode'][0]
#     n_id = input_maf_df['Matched_Norm_Sample_Barcode'][0]

#     # mutation_profile[t_id] = n_id + ':' + row_count

#     mut_prof_df.loc[t_id] = [n_id, row_count]

# print(mut_prof_df)

# sns.set(font_scale=1)

# sns.barplot(data=mut_prof_df,
#             x=mut_prof_df.index,
#             y="mut_count",
#             hue="Origin")

# plt.legend(loc=2, prop={'size': 15})
# # plt.xticks(rotation=-90)
# plt.show()







root_dir = r"E:/stemcell/somatic_analysis/maf/vardict/DP_AF_filtered_maf/exonic_maf/mtmut"
input_format = r"*.maf"
grp_name = 'Caller'

mut_prof_df = pd.DataFrame(columns=[grp_name, 'mut_count'])

# tid, mut_count, grp

# print(os.listdir(root_dir))

for grp in os.listdir(root_dir):
    wd = os.path.join(root_dir, grp)
    input_lst = jm.get_input_path_list(wd, input_format, False)
    
    for i in range(len(input_lst)):
        input_maf_df = pd.read_csv(input_lst[i], sep='\t')
        # print(input_maf_df)
        row_count = input_maf_df.shape[0]
        t_id = input_maf_df['Tumor_Sample_Barcode'][0]
        mut_prof_df.loc[t_id] = [grp, row_count]

print(mut_prof_df)

sns.set(font_scale=1)

sns.barplot(data=mut_prof_df,
            x=mut_prof_df.index,
            y="mut_count",
            hue=grp_name)

plt.legend(loc=2, prop={'size': 15})
# plt.xticks(rotation=-90)
plt.show()

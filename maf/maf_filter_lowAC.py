

# maf파일 가져와서
# 뎁스, t ac, n ac정보에 대해 조건걸고
# drop 시키고 
# maf로 저장

import os
import pandas as pd 
from glob import glob


pd.set_option('display.max_seq_items', None)


input_dir = r'E:/UTUC_data/DH_ref/JM'
input_format = r'*.maf'

output_dir_name = r'ac_filtered'

min_total_depth = 30 # 30 이상
min_tumor_ac = 5 # 5 이상
max_normal_ac = 1 # 1 이하


output_dir = os.path.join(input_dir, output_dir_name)

if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)



input_maf_lst = glob(os.path.join(input_dir, input_format))
print(input_maf_lst)

for input_maf in input_maf_lst:

    maf_df = pd.read_csv(input_maf, sep='\t', low_memory=False)

    sample_name = maf_df['Tumor_Sample_Barcode'][0]
    print(sample_name)

    # print(maf_df.columns)
    print(maf_df.shape)

    filtering_idx = maf_df[(maf_df['t_depth'] < min_total_depth) | \
        (maf_df['t_alt_count'] < min_tumor_ac) | \
            (maf_df['n_alt_count'] > max_normal_ac)].index
    # print(filtering_idx)

    maf_df = maf_df.drop(filtering_idx)
    maf_df.reset_index(inplace=True, drop=True)

    print(maf_df.shape)

    output_name = sample_name + '_ac-filtered.maf'
    output_path = os.path.join(output_dir, output_name)

    maf_df.to_csv(output_path, index=False, sep='\t')

    print(maf_df)



    # break
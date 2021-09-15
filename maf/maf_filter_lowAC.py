

# maf파일 가져와서
# 뎁스, t ac, n ac정보에 대해 조건걸고
# drop 시키고 
# maf로 저장

import os
import pandas as pd 
from glob import glob


pd.set_option('display.max_seq_items', None)


# input_dir = r'D:/stemcell/hg38/passage_comp/hiPS29-A/filtered/p29/filter_som_germ_merge'
# input_dir = r'D:/stemcell/hg38/clone_comp/hiPS29/som_germ_merge'
# input_dir = r'E:/stemcell_ips/gdc/tech/29A/filtered/tech2/filter_som_germ_merge'
# input_dir = r'E:/stemcell_ips/gdc/passage/29A/29A_p29_muthap/filter_mut_hap_merge'
# input_dir = r'E:/stemcell_ips/gdc/clone/hips29/29E_muthap/filter_mut_hap_merge'
# input_dir = r'E:/UTUC_data/gdc_hg38/maf/1st_lynch'
input_dir = r'E:/stemcell_ips/gdc/clone/hips66/66C_muthap/filter_mut_hap_merge'



# input_format = r'*_filtered.maf'
input_format = r'*.maf'

target_sample_name = r'hiPS66-C' # Tumor_Sample_Barcode가 없는 maf파일인 경우

output_dir_name = r'DP_filtered_maf'
output_suffix = r'_filterComplete.maf'

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

    try:
        sample_name = maf_df['Tumor_Sample_Barcode'][0]
        print(sample_name)
    
    # Tumor_Sample_Barcode가 없는경우 커버
    except KeyError:
        sample_name = target_sample_name
        print(sample_name)

    # print(maf_df.columns)
    print(maf_df.shape)

    try:
        filtering_idx = maf_df[(maf_df['t_depth'] < min_total_depth) | \
            (maf_df['t_alt_count'] < min_tumor_ac) | \
                (maf_df['n_alt_count'] > max_normal_ac)].index
    
    # (tumor only sample일때) normal 정보가 없는경우임
    except KeyError:
        filtering_idx = maf_df[(maf_df['t_depth'] < min_total_depth) | \
        (maf_df['t_alt_count'] < min_tumor_ac)].index

            
    # print(filtering_idx)

    maf_df = maf_df.drop(filtering_idx)
    maf_df.reset_index(inplace=True, drop=True)

    print(maf_df.shape)

    output_name = sample_name + output_suffix
    output_path = os.path.join(output_dir, output_name)

    maf_df.to_csv(output_path, index=False, sep='\t')

    print(maf_df)



    # break
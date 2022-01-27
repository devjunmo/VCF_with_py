

# maf파일 가져와서
# 뎁스, t ac, n ac정보에 대해 조건걸고
# drop 시키고 
# maf로 저장

import os
import pandas as pd 
from glob import glob


pd.set_option('display.max_seq_items', None)

# input_dir = r'E:/stemcell_ips/gdc/clone/hips29/29E_muthap/filter_mut_hap_merge'
# input_dir = r'E:/UTUC_data/gdc_hg38/maf/3rd'
# input_dir = r'E:/stemcell_ips/gdc/tech/29A/filtered/tech1/filter_som_germ_merge/DP_filtered_maf'
# input_dir = r'D:/junmo/wd/utuc/maf/rmhd_maf/mutect2/sample2/'
# input_dir = r'E:/stemcell_ips/gdc/clone/hips66/66C_muthap/filter_mut_hap_merge'
# input_dir = r'E:/UTUC_data/gdc_hg38/germline/maf/LG_N'
input_dir = r'E:/stemcell/VQSR_MAF'



# input_format = r'*_filtered.maf'
input_format = r'*.maf'

target_sample_name = r'hiPS36-C' # Tumor_Sample_Barcode가 없는 maf파일인 경우

# output_dir_name = r'DP_filtered_maf'
output_dir_name = r'DP_AF_filtered_maf'

output_suffix = r'_DP_filtered.maf'

min_total_depth = 30 # 30 이상
min_tumor_ac = 5 # 5 이상
max_normal_ac = 1 # 1 이하

is_filter_lowAF = True
af_threshold = 0.05 # 이 수치 미만 컷

# do_not_filter_gene_lst = ['HRAS'] # 여기에 등록된 유전자 레코드는 필터되지 않는다
do_not_filter_gene_lst = []


output_dir = os.path.join(input_dir, output_dir_name)

if os.path.isdir(output_dir) is False:
    os.mkdir(output_dir)



input_maf_lst = glob(os.path.join(input_dir, input_format))
# print(input_maf_lst)


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
    # print(maf_df.shape)
    # print(maf_df['Hugo_Symbol'])
    # print(maf_df['Hugo_Symbol'] == 'HRAS')
    # exit(0)

    try:
        filtering_idx = maf_df[
            ((maf_df['t_depth'] < min_total_depth) | \
            (maf_df['t_alt_count'] < min_tumor_ac) | \
            (maf_df['n_alt_count'] > max_normal_ac))].index
    
    # (tumor only sample일때) normal 정보가 없는경우임
    except KeyError:
        filtering_idx = maf_df[(maf_df['t_depth'] < min_total_depth) | \
        (maf_df['t_alt_count'] < min_tumor_ac)].index

    # print(filtering_idx)

    unfilter_gene_idx = maf_df[maf_df['Hugo_Symbol'].isin(do_not_filter_gene_lst)].index
    # print(unfilter_gene_idx)
    # print(maf_df.loc[unfilter_gene_idx])
    # exit(0)

    diff_idx = filtering_idx.difference(unfilter_gene_idx)


    maf_df = maf_df.drop(diff_idx)
    maf_df.reset_index(inplace=True, drop=True)

    diff_idx = None


    if is_filter_lowAF:
        
        low_af_idx = maf_df[((maf_df['t_alt_count'] / maf_df['t_depth']) < af_threshold)].index

        # print(low_af_idx)

        # maf_df = maf_df.drop(low_af_idx)
        # maf_df.reset_index(inplace=True, drop=True)
        unfilter_gene_idx = maf_df[maf_df['Hugo_Symbol'].isin(do_not_filter_gene_lst)].index

        diff_idx = low_af_idx.difference(unfilter_gene_idx)


        maf_df = maf_df.drop(diff_idx)
        maf_df.reset_index(inplace=True, drop=True)




    print(maf_df.shape)
    # print(maf_df)

    

    output_name = sample_name + output_suffix
    output_path = os.path.join(output_dir, output_name)

    maf_df.to_csv(output_path, index=False, sep='\t')

    



    # break
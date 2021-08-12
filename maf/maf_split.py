
# 특정 컬럼 기준으로 merged_maf를 쪼개서 여러개의 maf파일로 나눠주는 코드 

import pandas as pd
import os


# colname 생략 없이 출력
pd.set_option('display.max_seq_items', None)
# row 생략 없이 출력
# pd.set_option('display.max_rows', None)
# col 생략 없이 출력
# pd.set_option('display.max_columns', None)

root_dir = r'E:/UTUC_data/DH_ref'

merged_maf_name = 'all.merged.maf.txt'

merged_maf_path = os.path.join(root_dir, merged_maf_name)

merged_maf_df = pd.read_csv(merged_maf_path, sep='\t')

# print(merged_maf_df.columns)


print(merged_maf_df['Tumor_Sample_Barcode'].unique())

sample_id_lst = merged_maf_df['Tumor_Sample_Barcode'].unique()


split_maf_lst = []

for sam_id in sample_id_lst:

    tmp_maf = merged_maf_df[merged_maf_df['Tumor_Sample_Barcode'] == sam_id]

    # print(tmp_maf)
    # print(tmp_maf['Tumor_Sample_Barcode'])

    output_name = sam_id + '.maf'

    output_path = os.path.join(root_dir, output_name)

    tmp_maf.to_csv(output_path, sep='\t', index=False)

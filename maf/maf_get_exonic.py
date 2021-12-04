
import pandas as pd
import os
from glob import glob

input_maf_dir = r'E:/UTUC_data/gdc_hg38/maf/1st_lynch/DP_AF_filtered_maf'
input_format = '*.maf'

output_suffix = '_exonic.maf'

coding_region_lst = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Del', \
                    'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Silent', 'Splice_Site', 'Translation_Start_Site']



input_maf_lst = glob(os.path.join(input_maf_dir, input_format))


for i in range(len(input_maf_lst)):
    print(input_maf_lst[i])
    file_path = input_maf_lst[i]

    f_path, f = os.path.split(file_path)
    print(f)
    f_name = f.split('.')[0]
    f_ext = f.split('.')[1]
    print(f_name)

    input_maf_df = pd.read_csv(file_path, sep='\t')
    print(input_maf_df.shape)

    exoic_df = input_maf_df[input_maf_df['Variant_Classification'].isin(coding_region_lst)]

    out_name = f_name + output_suffix

    print(out_name)

    out_path = os.path.join(input_maf_dir, out_name)

    exoic_df.to_csv(out_path, sep='\t', index=False)
    
    # break
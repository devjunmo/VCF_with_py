
# Filter에 AC관련 필터 tag를 추가해 주는 코드
# for opponent maf

import pandas as pd 
import os
from glob import glob


# input_dir = r'D:/stemcell/hg38/clone_comp/hiPS29/som_germ_merge/unfiltered'
# input_dir = r'E:/stemcell_ips/gdc/tech/29A/unfiltered/tech2/unfilter_som_germ_merge'
# input_dir = r'E:/stemcell_ips/gdc/tech/29B/29B_tech2_muthap/unfilter_mut_hap_merge'
# input_dir = r'E:/stemcell_ips/gdc/passage/29A/29A_p29_muthap/unfilter_mut_hap_merge'
# input_dir = r'E:/stemcell_ips/gdc/passage/29B/29B_p30_muthap/unfilter_mut_hap_merge'
# input_dir = r'E:/stemcell_ips/gdc/tech/29B/29B_tech2_muthap/unfilter_mut_hap_merge'
input_dir = r'E:/stemcell_ips/gdc/clone'


input_format = r'*.xlsx'

input_data_lst = glob(os.path.join(input_dir, input_format))

output_suffix = '_DPTag.xlsx'

# filter_lst = ['t_depth', 't_ref_count', 't_alt_count']
filter_dict = {'t_depth':[30, ';dp30'], 't_alt_count':[5, ';v5']}
print(filter_dict.keys())

for input_data in input_data_lst:

    # print(os.path.splitext(os.path.basename(input_data)))

    # break

    input_df = pd.read_excel(input_data, sheet_name='Gene data')
    input_df.reset_index(inplace=True, drop=True)

    input_df['tmp'] = ''

    for filter_name in filter_dict.keys():
        
        for filter_index, filter_value in input_df[filter_name].iteritems():
            # print(filter_index)
            # print(filter_value)

            if filter_value < filter_dict[filter_name][0]:
                # print(filter_dict[filter_name][0])

                input_df.loc[filter_index, 'tmp'] += filter_dict[filter_name][1]

    # print(input_df)

    input_df['FILTER'] = input_df[['FILTER', 'tmp']].apply(lambda row: ''.join(row.values.astype(str)), axis=1)

    input_df.pop('tmp')

    print(input_df)

    f_name = os.path.splitext(os.path.basename(input_data))[0] + output_suffix

    output_path = os.path.join(input_dir, f_name)

    input_df.to_excel(output_path, sheet_name='Gene data', index=False)





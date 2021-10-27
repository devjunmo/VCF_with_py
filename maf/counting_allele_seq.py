import pandas as pd
import os
import random
from matplotlib import pyplot as plt
import math

# long indel의 seq length 분포를 보기 위한 코드

input_maf_path = r'E:/UTUC_data/gdc_hg38/maf/3rd/DP_AF_filtered_maf/utuc_3rd_gdc_AF_filter_apply.xlsx'
# input_maf_path = r'E:/UTUC_data/gdc_hg38/maf/2nd/DP_AF_filtered_maf/utuc_2nd_gdc_acfilter.xlsx'

alt_1_col_name = r'Alt1'
alt_2_col_name = r'Alt2'


input_df = pd.read_excel(input_maf_path, sheet_name='Gene data')

# print(input_df.loc[:, [alt_1_col_name, alt_2_col_name]])

alt_df = input_df.loc[:, [alt_1_col_name, alt_2_col_name]]

alt_df['length'] = 0

# print(alt_df)

# for index_name, value in df['age'].iteritems():
#     print("인덱스 타입", type(index_name))
#     print("index명:", index_name)
#     print("value 타입:", type(value))  # <class 'pandas.core.series.Series'>
#     print("값: ", value)

len_lst = []

for idx, val in alt_df[alt_1_col_name].iteritems():
    # print(idx)
    # print(val)

    try:
        if val == '-':
            seq = alt_df.loc[idx, alt_2_col_name]
            seq_len = len(seq)
            len_lst.append(seq_len)
        
        else:
            seq = val
            seq_len = len(seq)
            len_lst.append(seq_len)
    except TypeError:
        print(f'{seq}: typeError exception --> Program ignores {seq}')
        len_lst.append(seq)


# print(len_lst)
len_lst_rm_nan = [x for x in len_lst if math.isnan(x) == False]
print(len_lst_rm_nan)
print(len(len_lst_rm_nan))

print(len_lst_rm_nan.index(32))

exit()
len_lst_rm_nan.sort()
print(len_lst_rm_nan)

plt.bar(range(len(len_lst_rm_nan)), len_lst_rm_nan)

plt.yticks(range(1, len_lst_rm_nan[-1], 5))

plt.grid(True, axis='y')

plt.show()

import venn
import os
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import numpy as np


# input_dir = r'E:/UTUC_data/WES/rmhd_maf/mutect2/sample2/'
input_dir = r'E:/UTUC_data/DH_ref/'
input_format = r'*.maf'

input_vcf_dir = r'E:/UTUC_data/VCF/rm_hd/mutect2/sample2/'
input_vcf_format = r'*.vcf'

vcf_header = ['Chromosome', 'Start_Position', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Normal', 'Tumor']

pair_info = r'E:/UTUC_data/utuc_NT_pair.csv'


pair_df = pd.read_csv(pair_info)
pair_df.set_index('Tumor', inplace=True)
pair_dict = pair_df.to_dict('index')


coding_region_lst = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', \
                'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', \
                    'Translation_Start_Site', 'Splice_Site']


print(pair_dict)


input_lst = glob(input_dir + input_format)
input_vcf_lst = glob(input_vcf_dir + input_vcf_format)

print(input_lst)
# print(input_vcf_lst)

# exit(0)

set_list = []

def get_normal_ac(sample_col):
    # print(sample_col)
    ac = int(sample_col.split(':')[1].split(',')[-1])
    # print(ac)
    # exit(0)
    if ac < 2:
        return_ac = ac
    else:
        return_ac = np.nan
    
    return return_ac

def get_tumor_ac(sample_col):
    # print(sample_col)
    ac = int(sample_col.split(':')[1].split(',')[-1])
    # print(ac)
    # exit(0)
    if ac >= 5:
        return_ac = ac
    else:
        return_ac = np.nan
    
    return return_ac

def get_ref_len(sample_col):
    length = len(sample_col)
    return length


for i in range(len(input_lst)):
    input_maf = input_lst[i]
    # print(input_maf)
    t_name = input_maf.split('\\')[-1].split(r'.')[0].split(r'_')[1] # 20S-14292-A1-7
    tumor_grade = pair_dict[t_name]['Tumor_Grade'] # low
    
    sample_tag = t_name + '-' + tumor_grade

    maf_df = pd.read_csv(input_maf, sep='\t', low_memory=False)

    input_vcf = input_vcf_lst[i]

    vcf_df = pd.read_csv(input_vcf, sep='\t', low_memory=False, names=vcf_header)
    vcf_df = vcf_df[['Chromosome', 'Start_Position', 'REF', 'ALT', 'Normal', 'Tumor']]

    normal_ac = vcf_df['Normal'].map(get_normal_ac)
    tumor_ac = vcf_df['Tumor'].map(get_tumor_ac)
    # ref_len = vcf_df['REF'].map(get_ref_len)

    vcf_df['Normal_ac'] = normal_ac
    vcf_df['Tumor_ac'] = tumor_ac
    # vcf_df['ref_len'] = ref_len

    print(vcf_df)

    # vcf_df['Start_Position'] = vcf_df['Start_Position'] + vcf_df['ref_len'] - 1
    vcf_df = vcf_df.dropna(axis=0)

    print(vcf_df)
    print(vcf_df.index)
    filter_idx = list(vcf_df.index)
    # exit()

    # print(maf_df)
    # print(len(vcf_df['Start_Position'].unique()))

    

    # print(t_name)
    # print(tumor_grade)

    # print(maf_df.columns)

    # break

    maf_df = maf_df[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', \
                     'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Variant_Classification']]

    print(maf_df)
    maf_df = maf_df.loc[filter_idx, :]

    # exit()

    # join_key = ['Chromosome', 'Start_Position']

    # joined_df = pd.merge(maf_df, vcf_df, how='left', on = join_key)

    # print(joined_df)
    
    # joined_df = joined_df.dropna(axis=0)
    # print(joined_df)

    # maf_df = joined_df

    maf_df.reset_index(inplace=True, drop=True)


    # exit(0)
    
    # coding region이 아닌 부분은 제외

    print(maf_df['Variant_Classification'].unique())
    print(pd.value_counts(maf_df['Variant_Classification']))

    coding_idx_lst = []

    for cd in coding_region_lst:
        idx = list(maf_df[maf_df['Variant_Classification'] == cd].index)
        coding_idx_lst.append(idx)

    coding_idx_lst = list(itertools.chain(*coding_idx_lst))
    coding_idx_lst.sort()
    # print(maf_df['Variant_Classification'].unique())
    
    maf_df = maf_df.iloc[coding_idx_lst, ]
    maf_df.reset_index(inplace=True, drop=True)

    # print(maf_df)



    # break
    
    # print(maf_df)

    sample_point_set = set()

    for row in maf_df.itertuples(index=False, name=None):
        # print(row)
        # print(type(row))
        sample_point_set.add(row)
        # break
    
    set_list.append([sample_tag, sample_point_set])
    # print(set_list)


# print(len(set_list)) # 6 샘플
# print(set_list[0][0]) # 샘플 태그
# print(set_list[0][1]) # 세트리스트 값



# print([f'{set_list[0][0]}', f'{set_list[1][0]}', f'{set_list[2][0]}', f'{set_list[3][0]}', f'{set_list[4][0]}', f'{set_list[5][0]}'])

# 6개일때

# labels = venn.get_labels([set_list[0][1], set_list[1][1], set_list[2][1], set_list[3][1], set_list[4][1], set_list[5][1]], \
#                             fill=['number']) # set으로 받아야함. list안됨
# # ['20S-14292-A1-7-low', '20S-31099-A4-14-low', '20S-31099-A4-1-high', '20S-31099-A4-2-high', '20S-31099-A4-3-low', '20S-31099-A4-5-intermediate']


# # venn.venn3(labels, names=[f'{set_list[0][0]}', f'{set_list[1][0]}', f'{set_list[2][0]}', \
# #                             f'{set_list[3][0]}', f'{set_list[4][0]}', f'{set_list[5][0]}'])

# venn.venn3(labels, names=['UTUC1-bx', 'UTUC1-LG2', 'UTUC1-HG1', 'UTUC1-HG2', 'UTUC1-LG3', 'UTUC1-Intermediate'])


# print(type(labels))



# exit(0)


# 5개일때

# ['20S-82978-A2-8-high', '20S-82978-A3-10-high', '20S-82978-A3-15-low', '20S-82978-A5-12-low', '20S-82978-A5-13-intermediate']

# print([f'{set_list[0][0]}', f'{set_list[1][0]}', f'{set_list[2][0]}', f'{set_list[3][0]}', f'{set_list[4][0]}'])

labels = venn.get_labels([set_list[0][1], set_list[1][1], set_list[2][1], set_list[3][1], set_list[4][1]], \
                            fill=['number']) # set으로 받아야함. list안됨

# venn.venn3(labels, names=[f'{set_list[0][0]}', f'{set_list[1][0]}', f'{set_list[2][0]}', \
#                             f'{set_list[3][0]}', f'{set_list[4][0]}'])

venn.venn3(labels, names=['UTUC2-HG1', 'UTUC2-HG2', 'UTUC2-LG2', 'UTUC2-LG1', 'UTUC2-Intermediate'])




plt.show()


exit(0)


from itertools import chain, combinations
 

def powerset(iterable):

    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))
 


# print(len(list(powerset(set_list)))) # 32

# print(len(list(powerset(set_list))[0])) # 0
# print(len(list(powerset(set_list))[1])) # 1
# print(len(list(powerset(set_list))[10])) # 2
# print(len(list(powerset(set_list))[31])) # 5



sub_lst = list(powerset(set_list)) # [(0), ..,   ( [ tag, {(),()..()} ], [ ] ... [ ] ),    () ..., ()]

# print(len(sub_lst[10]))

print('----------------------')
# for i in range(len(sub_lst)):
#     print(len(sub_lst[i]))


for i in range(len(sub_lst)):
    try:
        # print(len(sub_lst[0]))
        if len(sub_lst[0]) < 2:
            del sub_lst[0]
    except IndexError as e:
        break


print('----------------------')



# print(type(sub_lst[0])) # 튜플 <class 'tuple'> / 이 튜플 안에 [태그, 튜플셋]형태의 N개의 리스트가 부분집합으로써 들어있다. 
#                         # 현재 0번째에는 최소 매칭 부분집합인 2개 샘플이 들어가있음. 
#                         # 0번째를 포함하여 2개로 짝지을수 있는 모든 경우의수는 0, 1, 2...번째에 쭉 들어가있음.

# print(type(sub_lst[0][0])) # list <class 'list'> 부분집합 구성품. tag와 튜플세트로 구성.
#                             # sub_lst[i][j] = 샘플 꾸러미 i번째의 j번째 구성품 

# print(sub_lst[0][0][0]) # 20S-14292-A1-7-low
# print(type(sub_lst[0][0][1])) # 튜플 세트 <class 'set'>  # 연산 해야하는 




# 일단 자료구조를 tag: set의 딕셔너리 폼으로 바꾸자
# 그리고 태그:set의 딕셔너리 N개를 지정해두고
# 만약 ABCDE샘플의 CD specific region을 구하고 싶다면, CD를 &연산 한거에 CD에 없는 ABD의 합집합을 빼주면 됨.


for i in range(len(sub_lst)):

    isec_data = eval("&".join([f"sub_lst[{i}][{j}][1]" for j in range(len(sub_lst[i]))]))

    # concat_tag = eval("_".join([f"sub_lst[{i}][{j}][0]" for j in range(len(sub_lst[i]))]))

    print(len(isec_data))
    # print(isec_data)

    # print(sub_lst[i][0][0]) # 20S-14292-A1-7-low
    # print(sub_lst[i][1][0]) # 20S-31099-A4-14-low

    # print(len(sub_lst[i][0][1])) # 818
    # print(len(sub_lst[i][1][1])) # 1033
    # print(concat_tag)
    # [print(j) for j in range(len(sub_lst[i]))]
    # print(sub_lst[i])
    # break

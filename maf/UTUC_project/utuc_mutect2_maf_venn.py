import venn
import os
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import numpy as np


# input_dir = r'E:/UTUC_data/WES/maf/mutect2/sample2/'
# input_dir = r'E:/UTUC_data/gdc_hg38/maf/1st_lynch/DP_filtered_maf/'
# input_dir = r'E:/UTUC_data/gatk_hg38/DH_ref/JM/ac_filtered/sample2/'
# input_dir = r'E:/UTUC_data/gdc_hg38/maf/2nd/DP_AF_filtered_maf/'
input_dir = r'D:/junmo/wd/utuc/maf/rmhd_maf/mutect2/sample2/'
# input_dir = r'E:/UTUC_data/WES/rmhd_maf/mutect/mutect2/filtered_maf/sample2/'
input_format = r'*.maf'
save_gene_df_path = input_dir + r'utuc_2nd_gdc_acfilter.xlsx'
# save_gene_df_path = input_dir + r'utuc_3rd_gdc_AF_filter_apply.xlsx'
# save_gene_df_path = input_dir + r'utuc_sample2.xlsx'

pair_info = r'E:/UTUC_data/utuc_NT_pair_ver_210910.csv'


pair_df = pd.read_csv(pair_info)
pair_df.set_index('Tumor', inplace=True)
pair_dict = pair_df.to_dict('index')




# coding_region_lst = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', \
#                 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', \
#                     'Translation_Start_Site', 'Splice_Site']

coding_region_lst = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Del', \
                    'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Silent', 'Splice_Site']


# print(pair_dict)


input_lst = glob(input_dir + input_format)

print(input_lst)

set_list = []

for i in range(len(input_lst)):
    input_maf = input_lst[i]
    # print(input_maf)
    # t_name = input_maf.split('\\')[-1].split(r'.')[0].split(r'_')[1] # 20S-14292-A1-7
    t_name = input_maf.split('\\')[-1].split(r'.')[0].split(r'_')[0] # filtered_mutect2_20S-14292-A1-7 -> 20S-14292-A1-7
    tumor_grade = pair_dict[t_name]['Tumor_Grade'] # low
    
    sample_tag = t_name + '-' + tumor_grade

    maf_df = pd.read_csv(input_maf, sep='\t', low_memory=False)

    print(maf_df.shape)

    # exit(0)

    # print(t_name)
    # print(tumor_grade)

    # print(maf_df.columns)
    # print(maf_df['Variant_Classification'].unique())

    # break

    maf_df = maf_df[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', \
                     'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Variant_Classification']]
    
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
    
    set_list.append([sample_tag, sample_point_set]) # 내부를 리스트로
    # set_list.append((sample_tag, sample_point_set)) # 내부를 튜플로
    # print(set_list)


# print(len(set_list)) # 6 샘플
# print(set_list[0][0]) # 샘플 태그
# print(set_list[0][1]) # 세트리스트 값


# # 7개 일 때 -> X

# labels = venn.get_labels([set_list[0][1], set_list[1][1], set_list[2][1], set_list[3][1], set_list[4][1], set_list[5][1], set_list[6][1]], \
#                             fill=['number']) # set으로 받아야함. list안됨

# venn.venn3(labels, names=[f'{set_list[0][0]}', f'{set_list[1][0]}', f'{set_list[2][0]}', \
#                             f'{set_list[3][0]}', f'{set_list[4][0]}', f'{set_list[5][0]}', f'{set_list[6][0]}'])

# # 6개일때

# labels = venn.get_labels([set_list[0][1], set_list[1][1], set_list[2][1], set_list[3][1], set_list[4][1], set_list[5][1]], \
#                             fill=['number']) # set으로 받아야함. list안됨

# venn.venn3(labels, names=[f'{set_list[0][0]}', f'{set_list[1][0]}', f'{set_list[2][0]}', \
#                             f'{set_list[3][0]}', f'{set_list[4][0]}', f'{set_list[5][0]}'])


# print(type(labels))



# exit(0)


# 5개일때


labels = venn.get_labels([set_list[0][1], set_list[1][1], set_list[2][1], set_list[3][1], set_list[4][1]], \
                            fill=['number']) # set으로 받아야함. list안됨

venn.venn3(labels, names=[f'{set_list[0][0]}', f'{set_list[1][0]}', f'{set_list[2][0]}', \
                            f'{set_list[3][0]}', f'{set_list[4][0]}'])




plt.show()


# exit(0)

############################# 여기서 부터 클래스화 시켜야 함 ################################################


from itertools import chain, combinations
 

def powerset(iterable):

    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))
 


# print(len(list(powerset(set_list)))) # 32

# print(len(list(powerset(set_list))[0])) # 0
# print(len(list(powerset(set_list))[1])) # 1
# print(len(list(powerset(set_list))[10])) # 2
# print(len(list(powerset(set_list))[31])) # 5


# 모든 경우 부분집합 구하기
sub_lst = list(powerset(set_list)) # [(0), ..,   ( [ tag, {(),()..()} ], [ ] ... [ ] ),    () ..., ()]

# print(len(sub_lst[10]))

print('----------------------')
# for i in range(len(sub_lst)):
#     print(len(sub_lst[i]))


for i in range(len(sub_lst)):
    try:
        # print(len(sub_lst[0]))
        if len(sub_lst[0]) < 1:
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

# i는 서브 리스트의 요소 수
# j는 한 요소 내 세트 수

whole_case_lst = list(sub_lst[-1])
# print(len(whole_case_lst))
# print(type(whole_case_lst))


final_df = pd.DataFrame(columns=['CaseNumber', 'Gene', 'Chr', 'Start', 'End', \
    'Ref', 'Alt1', 'Alt2', 'Type'])

# info_df = pd.DataFrame(columns=['Case', 'Symbol_number'])
info_df = pd.DataFrame(columns=['Case_number', 'Case'])


symbol_num = 1


for i in range(len(sub_lst)):

    etc_case = []

    # print(whole_case)

    for one_case in whole_case_lst:
        if one_case not in sub_lst[i]:
            etc_case.append(one_case)   


    # print(etc_case)
    # print(len(etc_case))

    # print(whole_case_lst)

    # sub_lst_whole_set = set(sub_lst[-1])
    # sub_lst_current_set = set(sub_lst[i])

    # etc_set = sub_lst_whole_set - sub_lst_current_set

    # print(etc_set)

    name_lst = []

    for n in range(len(sub_lst[i])):
        # print(sub_lst[i][n][0])
        name_lst.append(sub_lst[i][n][0])

    concat_name = '___'.join(name_lst)
    # print(concat_name)

    # info_df = info_df.append(pd.Series([symbol_num, concat_name], index=info_df.columns), ignore_index=True)

    # symbol_num = symbol_num + 1



    isec_data = eval("&".join([f"sub_lst[{i}][{j}][1]" for j in range(len(sub_lst[i]))]))
    # print(len(isec_data))

    try:
        etc_isec_data = eval("|".join([f"etc_case[{j}][1]" for j in range(len(etc_case))]))

        # print(len(etc_isec_data))

        tmp_isec = isec_data & etc_isec_data

        # print(len(tmp_isec))

        specific_gene_set = isec_data - tmp_isec

    except SyntaxError as e:
        specific_gene_set = isec_data
    


    # print(len(specific_gene_set))
    # print(specific_gene_set)

    # print(final_df)

    for row in specific_gene_set:
        # print(list(row))

        input_row = list(row)
        # print(input_row)
        # input_row.insert(0, concat_name)
        input_row.insert(0, str(symbol_num))
        # print(input_row)
        # print(type(input_row))
        final_df = final_df.append(pd.Series(input_row, index=final_df.columns), ignore_index=True)
        info_df = info_df.append(pd.Series([symbol_num, concat_name], index=info_df.columns), ignore_index=True)
    
    if len(specific_gene_set) == 0: # 교집합이 없을때 처리
            input_row = [np.nan]*(len(final_df.columns)-1)
            input_row.insert(0, str(symbol_num))
            final_df = final_df.append(pd.Series(input_row, index=final_df.columns), ignore_index=True)
    

    symbol_num = symbol_num + 1


    # concat_tag = eval("_".join([f"sub_lst[{i}][{j}][0]" for j in range(len(sub_lst[i]))]))
    # print(isec_data)

    # print(len(isec_data))
    # print(isec_data)

    # print(sub_lst[i][0][0]) # 20S-14292-A1-7-low
    # print(sub_lst[i][1][0]) # 20S-31099-A4-14-low

    # print(len(sub_lst[i][0][1])) # 818
    # print(len(sub_lst[i][1][1])) # 1033
    # print(concat_tag)
    # [print(j) for j in range(len(sub_lst[i]))]
    # print(sub_lst[i])
    # break



print(final_df)

info_df.drop_duplicates(['Case_number', 'Case'], inplace=True)

writer = pd.ExcelWriter(save_gene_df_path, engine='xlsxwriter')


final_df.to_excel(writer, sheet_name='Gene data', index=False, na_rep='NaN')

info_df.to_excel(writer, sheet_name='info', index=False)

writer.save()
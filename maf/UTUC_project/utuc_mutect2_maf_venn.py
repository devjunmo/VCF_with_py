import venn
import os
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
import itertools


input_dir = r'E:/UTUC_data/WES/maf/mutect2/sample1/'
input_format = r'*.maf'

pair_info = r'E:/UTUC_data/utuc_NT_pair.csv'


pair_df = pd.read_csv(pair_info)
pair_df.set_index('Tumor', inplace=True)
pair_dict = pair_df.to_dict('index')


coding_region_lst = ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', \
                'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', \
                    'Translation_Start_Site', 'Splice_Site']


print(pair_dict)


input_lst = glob(input_dir + input_format)

print(input_lst)

set_list = []

for i in range(len(input_lst)):
    input_maf = input_lst[i]
    # print(input_maf)
    t_name = input_maf.split('\\')[-1].split(r'.')[0].split(r'_')[1] # 20S-14292-A1-7
    tumor_grade = pair_dict[t_name]['Tumor_Grade'] # low
    
    sample_tag = t_name + '-' + tumor_grade

    maf_df = pd.read_csv(input_maf, sep='\t', low_memory=False)

    # print(t_name)
    # print(tumor_grade)

    # print(maf_df.columns)

    # break

    maf_df = maf_df[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', \
                     'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Variant_Classification']]
    
    # coding region이 아닌 부분은 제외

    # print(maf_df['Variant_Classification'].unique())

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



# # 6개일때

# labels = venn.get_labels([set_list[0][1], set_list[1][1], set_list[2][1], set_list[3][1], set_list[4][1], set_list[5][1]], \
#                             fill=['number']) # set으로 받아야함. list안됨

# venn.venn3(labels, names=[f'{set_list[0][0]}', f'{set_list[1][0]}', f'{set_list[2][0]}', \
#                             f'{set_list[3][0]}', f'{set_list[4][0]}', f'{set_list[5][0]}'])



# # 5개일때

# # labels = venn.get_labels([set_list[0][1], set_list[1][1], set_list[2][1], set_list[3][1], set_list[4][1]], \
# #                             fill=['number']) # set으로 받아야함. list안됨

# # venn.venn3(labels, names=[f'{set_list[0][0]}', f'{set_list[1][0]}', f'{set_list[2][0]}', \
# #                             f'{set_list[3][0]}', f'{set_list[4][0]}'])




# plt.show()


# exit(0)


from itertools import chain, combinations
 

def powerset(iterable):

    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))
 


print(len(list(powerset(set_list)))) # 32

print(len(list(powerset(set_list))[0])) # 0
print(len(list(powerset(set_list))[1])) # 1
print(len(list(powerset(set_list))[10])) # 2
print(len(list(powerset(set_list))[31])) # 5



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



print(type(sub_lst[0])) # 튜플 <class 'tuple'> / 이 튜플 안에 [태그, 튜플셋]형태의 N개의 리스트가 부분집합으로써 들어있다. 
                        # 현재 0번째에는 최소 매칭 부분집합인 2개 샘플이 들어가있음. 
                        # 0번째를 포함하여 2개로 짝지을수 있는 모든 경우의수는 0, 1, 2...번째에 쭉 들어가있음.

print(type(sub_lst[0][0])) # list <class 'list'> 부분집합 구성품. tag와 튜플세트로 구성.
                            # sub_lst[i][j] = 샘플 꾸러미 i번째의 j번째 구성품 

print(sub_lst[0][0][0]) # 20S-14292-A1-7-low
print(type(sub_lst[0][0][1])) # 튜플 세트 <class 'set'>  # 연산 해야하는 




for i in range(len(sub_lst)):
    print(i)


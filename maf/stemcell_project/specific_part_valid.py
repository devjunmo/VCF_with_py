
# 벤다이어그램에서 스페시픽한 부분을 상대방 원본 maf파일과 비교하는 코드

# input: 생성된 genelist 엑셀파일, case넘버와 매칭되는 maf 경로가 적힌 csv파일(maf가 여러개일때는 ';'로 구분) 
# maf file 매칭되는게 없을때는 maf 경로자리에 "N" 입력
# format: {col1: CaseNum, col2: oppMafPath}



from numpy.core.fromnumeric import shape
import pandas as pd
import os
import numpy as np


root_dir = r'E:/stemcell_ips/gdc/tech/29A/filtered/tech_comp/exclude_filterTag_tech_comp'

# variant_info = r'D:/stemcell/hg38/passage_comp/hiPS29-A/filtered/passage_comp/exclude_filterTag_passage_comp/hiPS29-A_genelst_filtered_passageComp.xlsx'
variant_info = os.path.join(root_dir, 'hiPS29-A-p49_varinat_filtered.xlsx')

data_sheet_name = 'Gene data'
info_sheet_name = 'info'

# pair_info = r'E:/stemcell_ips/somatic_call/vardict/hiPS29/tech_compare/B_p49/tech_comp_pair_info.csv'
# pair_info = r'D:/stemcell/hg38/tech_comp/ips29-B-p49/filtered/final_tech_comp/exclude_filterTag_variant/hiPS29-B-tech_comp_pair_info.csv'
# pair_info = r'D:/stemcell/hg38/passage_comp/p49-2/filtered/final_tech_comp/exclude_filterTag_variant/hiPS29-2_p49_comp_pair-info.csv'
pair_info = os.path.join(root_dir, 'hiPS29-A-tech_comp_pair_info.csv')

# output_path = r'E:/stemcell_ips/HAP_VAD_compare/hIPS29/technical/B-2/hIPS29-B-2-specific_caller_compare.xlsx'
# output_path = r'D:/stemcell/hg38/tech_comp/ips29-B-p49/filtered/final_tech_comp/exclude_filterTag_variant/hiPS29-B-p49_techComp_filterComplete_specific.xlsx'
output_path = os.path.join(root_dir, 'hiPS29-A-tech_comp_specific.xlsx')

# case_num_name = ['']


var_df = pd.read_excel(variant_info, sheet_name = data_sheet_name)
var_df_info = pd.read_excel(variant_info, sheet_name = info_sheet_name)

# print(var_df.head())
print(var_df.shape)


pair_df = pd.read_csv(pair_info)

# print(pair_df)

pair_df.set_index('CaseNum', inplace=True)
pair_dict = pair_df.to_dict('index')

print(pair_dict)

# print(pair_dict[1]['oppMafPath']) # 정상


# pair_dict의 key만 뽑아서, 그 키에 해당하는 var_df의 case num에 대해 df를 추출하고, 이너조인하고, 리스트에 저장 후 마지막에 리스트에 있는거 한번에 concat

split_df_lst = []

for case_key in pair_dict.keys():

    if pair_dict[case_key]['oppMafPath'] == 'N':
        target_var_df = var_df[var_df['Case_number'] == case_key]
        continue

    # print(case_key)

    opp_maf_paths = pair_dict[case_key]['oppMafPath']
    print(opp_maf_paths)
    opp_maf_path_lst = opp_maf_paths.split(';')

    res_df = None

    for i in range(len(opp_maf_path_lst)):
        
        opp_maf_path = opp_maf_path_lst[i]
    
        opp_maf_df = pd.read_csv(opp_maf_path, sep='\t', low_memory=False)

        # print(opp_maf_df)
        

        opp_maf_df = opp_maf_df[['Chromosome', 'Start_Position', 'End_Position', \
                                'Reference_Allele', 'Tumor_Seq_Allele2', \
                                'FILTER', 't_depth', 't_ref_count', 't_alt_count']]

        opp_maf_df.columns = ['Chr', 'Start', 'End', 'Ref', 'Alt2', f'opp_{i}_FILTER', \
                            f'opp_{i}_t_depth', f'opp_{i}_t_ref_count', f'opp_{i}_t_alt_count']

        # print(opp_maf_df.head())
        # print(opp_maf_df)

        target_var_df = var_df[var_df['Case_number'] == case_key]

        # print(target_var_df)

        target_var_df = target_var_df.sort_values(by=['Chr', 'Start'])
        target_var_df.reset_index(inplace=True, drop=True)

        print(target_var_df.head())
        print(target_var_df.tail())
        print(target_var_df.shape)

        if res_df is None:
            res_df = pd.merge(target_var_df, opp_maf_df, on=['Chr', 'Start', 'End', 'Ref', 'Alt2'], how='left')
            # print(res_df)
            # exit(0)
        else:
            tmp_df = pd.merge(target_var_df, opp_maf_df, on=['Chr', 'Start', 'End', 'Ref', 'Alt2'], how='left')
            res_df = pd.merge(res_df, tmp_df, on=['Case_number', 'Gene', 'Chr', 'Start', 'End', 'Ref', 'Alt2', \
                'Type', 'FILTER', 't_depth', 't_ref_count', 't_alt_count'], how='inner')

        # print(res_df.head())
        # print(res_df.tail())
        # print(res_df.shape)

        # exit()

        # res_df_NAN = res_df[res_df['opp_1_FILTER'].isnull()]


        # print(res_df_NAN.head())
        # print(res_df_NAN.tail())
        # print(res_df_NAN.shape)

        # res_df_NAN_dp30 = res_df_NAN[res_df_NAN['t_depth'] >= 30]

        # print(res_df_NAN_dp30.head())
        # print(res_df_NAN_dp30.tail())
        # print(res_df_NAN_dp30.shape)

        if i == (len(opp_maf_path_lst) - 1):
            split_df_lst.append(res_df)


        # break

# print(split_df_lst[2])
# print(split_df_lst[3])
# exit(0)


fin_df = pd.concat(split_df_lst)

print(fin_df.shape)


    

writer = pd.ExcelWriter(output_path, engine='xlsxwriter')

fin_df.to_excel(writer, sheet_name='Specific data', index=False, na_rep='NaN')

var_df_info.to_excel(writer, sheet_name='info', index=False, na_rep='NaN')

writer.save()







# var_df.to_excel(na_rep='NaN')


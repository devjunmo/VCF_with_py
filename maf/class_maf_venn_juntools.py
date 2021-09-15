import subprocess as sp
import os
from glob import glob
import pandas as pd
import natsort
import venn
import matplotlib.pyplot as plt
import itertools
import numpy as np
from itertools import chain, combinations


# 사용자에게 노출: venn 그리기, specific 


class MafVennJuntools():
    

    def __init__(self, _filter_maf_lst=[], _unfilter_maf_lst=[], \
        _coding_regions=[]):
        self.filter_maf_lst = _filter_maf_lst
        self.unfilter_maf_lst = _unfilter_maf_lst
        self.merged_maf = None
        self.coding_regions = _coding_regions


    def mk_setlist(self, _maf_lst, _exclude_filter_tag=True, \
        _is_just_exonic=True, _is_inc_germline=True):

        set_list = []

        # maf_count_dict = {sample_tag:{(mutid_key_tup):(count_value_tup)}}
        maf_count_dict = dict()

        for i in range(len(_maf_lst)):
            input_maf = _maf_lst[i]

            sample_tag = input_maf.split('\\')[-1].split(r'.')[0].split(r'_')[0] # hiPS66-C-P10_rmHd.maf -> hiPS66-C-P10
            # t_name = input_maf.split('\\')[-1].split(r'.')[0]

            print(sample_tag)

            maf_df_raw = pd.read_csv(input_maf, sep='\t', low_memory=False)

            print(maf_df_raw.shape)

            maf_df = maf_df_raw[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', \
                                    'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', \
                                    'FILTER', 't_depth', 't_ref_count', 't_alt_count']]

            del maf_df_raw
            
            # maf_count_dict[sample_tag] = {}
            # print(type(maf_df.itertuples(index=False, name=None)))

            
            ##############################################################

            # print(maf_df['Variant_Classification'].unique())
            # print(pd.value_counts(maf_df['Variant_Classification']))

            ##############################################################

            coding_idx_lst = []

            if _is_just_exonic:
                for cd in self.coding_regions:
                    idx = list(maf_df[maf_df['Variant_Classification'] == cd].index)
                    coding_idx_lst.append(idx)

                coding_idx_lst = list(itertools.chain(*coding_idx_lst))
                coding_idx_lst.sort()
                # print(maf_df['Variant_Classification'].unique())
                
                maf_df = maf_df.iloc[coding_idx_lst, ]
                maf_df.reset_index(inplace=True, drop=True)



            if _exclude_filter_tag:
                if _is_inc_germline:
                    non_pass_idx = maf_df[(maf_df['FILTER'] != 'PASS') & \
                        (maf_df['FILTER'] != 'common_variant') & (maf_df['FILTER'] != 'germline')].index
                else:
                    non_pass_idx = maf_df[(maf_df['FILTER'] != 'PASS') & (maf_df['FILTER'] != 'common_variant')].index
                maf_df = maf_df.drop(non_pass_idx)
                maf_df.reset_index(inplace=True, drop=True)

            # print(pd.value_counts(maf_df['Variant_Classification']))

            maf_df_base = maf_df[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', \
                                'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification']]

            maf_df_for_dict = maf_df[['Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', \
                                'Reference_Allele', 'Tumor_Seq_Allele2', 'Variant_Classification', \
                                    'FILTER', 't_depth', 't_ref_count', 't_alt_count']]

            del maf_df

            # print(maf_df)


            maf_count_dict[sample_tag] = dict()

            for raw_row in maf_df_for_dict.itertuples(index=False, name=None):
                
                # print(raw_row)
                mut_id = raw_row[:7]
                count_val = raw_row[-4:]

                # print(mut_id) # ('DDX11L1', '1', 14653, 14653, 'C', 'T', "3'Flank")
                # print(count_val) # ('PASS', 12, 9, 3)
                # exit(0)

                maf_count_dict[sample_tag][mut_id] = count_val


            sample_point_set = set()

            for row in maf_df_base.itertuples(index=False, name=None):

                sample_point_set.add(row)

            
            set_list.append([sample_tag, sample_point_set]) # 내부를 리스트로
            # set_list.append((sample_tag, sample_point_set)) # 내부를 튜플로

        return set_list



    def powerset(self, iterable):
    
        s = list(iterable)
        return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


    def mk_gene_list(self, _set_list):
    
        
        # 모든 경우 부분집합 구하기
        sub_lst = list(self.powerset(_set_list)) # [(0), ..,   ( [ tag, {(),()..()} ], [ ] ... [ ] ),    () ..., ()]

        # print(len(sub_lst[10]))

        print('----------------------')
        print('|   mk gene list ... |')
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


        

        whole_case_lst = list(sub_lst[-1])
        # print(len(whole_case_lst))
        # print(type(whole_case_lst))


        final_df = pd.DataFrame(columns=['Case_number', 'Gene', 'Chr', 'Start', 'End', \
                    'Ref', 'Alt2', 'Type', 'FILTER', 't_depth', 't_ref_count', 't_alt_count'])

        # info_df = pd.DataFrame(columns=['Case', 'Symbol_number'])
        info_df = pd.DataFrame(columns=['Case_number', 'Case'])


        symbol_num = 1


        for i in range(len(sub_lst)):

            etc_case = []

            # print(whole_case)

            for one_case in whole_case_lst:
                if one_case not in sub_lst[i]:
                    etc_case.append(one_case)   


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
                etc_isec_data = eval("|".join([f"etc_case[{j}][1]" for j in range(len(etc_case))])) # etc의 합집합

                # print(len(etc_isec_data))

                tmp_isec = isec_data & etc_isec_data

                # print(len(tmp_isec))

                specific_gene_set = isec_data - tmp_isec

            except SyntaxError as e:
                specific_gene_set = isec_data
            


            # print(len(specific_gene_set))
            # print(specific_gene_set)

            # print(final_df)

            # print(type(specific_gene_set)) # <class 'set'>
            # print(len(specific_gene_set)) # 길이 변화 x
            # exit(0)

            for row in specific_gene_set:
                # print(list(row))

                # print(row) 
                # print(type(row)) # tuple

                # dic_key = row[]

                input_row = list(row)
                # print(input_row)

                samp_name = concat_name.split('___')[0]

                # print(maf_count_dict[samp_name][mut_id]) # (2,0,2)
                # print(samp_name)


                t_count_val = list(maf_count_dict[samp_name][row])

                # print(t_count_val)
                # exit(0)

                input_row = input_row + t_count_val

                # print(input_row)

                # exit(0)

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



        final_df = final_df.sort_values(by=['Case_number', 'Chr', 'Start'])

        print(final_df)

        # exit(0)

        info_df.drop_duplicates(['Case_number', 'Case'], inplace=True)

        writer = pd.ExcelWriter(save_gene_df_path, engine='xlsxwriter')


        final_df.to_excel(writer, sheet_name='Gene data', index=False)

        info_df.to_excel(writer, sheet_name='info', index=False)

        writer.save()


    # def _mk_maf_subset(_maf_df, exclude_filter_tag=True, \
    #     is_just_exonic=True, is_inc_germline=True, \
    #     is_just_coding_regions=False):

    #     pass


    @staticmethod
    def draw_venn_diagram(_setlist_lst):
        pass

    
    

    def gene_lst_to_maf(self, _genelist_lst):
        
        return merged_maf_lst
        pass


    def test_inst_super_mtd(self):
        pass


    def __mk_filtered_maf(self):
        pass


    def mk_unfiltered_maf(self):
        pass


    def run_general():
        pass
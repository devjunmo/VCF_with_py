import subprocess as sp
import os
import sys
from glob import glob
import pandas as pd
import natsort
import venn
import matplotlib.pyplot as plt
import itertools
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))

from class_maf_venn_juntools import MafVennJuntools



# 사용자에게 노출: venn 그리기, specific 


class MafVennStemcell(MafVennJuntools):
    

    def __init__(self, _muthap_lst, _filter_maf_lst=[], _unfilter_maf_lst=[], \
        _coding_regions=[]):
        
        
        super(MafVennStemcell, self).__init__(_filter_maf_lst, _unfilter_maf_lst, \
            _coding_regions) # 부모 클래스 생성자 call

        self.muthap_lst = natsort.natsorted(_muthap_lst)

        self.muthap_setlist_lst_filter = []
        self.muthap_setlist_lst_unfilter = []

        self.muthap_gene_list_filter = None
        self.muthap_gene_list_unfilter = None

        self.muthap_merged_maf_filter = None
        self.muthap_merged_maf_unfilter = None
        


    def run_stemcell_pipe(self):

        self._mk_muthap_gene_list()

        self.muthap_merged_maf_filter = super().gene_lst_to_maf(self.muthap_gene_list_filter)
        del self.muthap_gene_list_filter

        dp_tag_unfilter_genelst = self.add_DPtag_unfilter(self.muthap_gene_list_unfilter)
        self.muthap_merged_maf_unfilter = super().gene_lst_to_maf(dp_tag_unfilter_genelst)
        del self.muthap_gene_list_unfilter



    def _mk_muthap_gene_list(self): # 디렉토리에 두개 파일로 시작 (mut, hap)

        maf1 = self.muthap_lst[0] # mut
        maf2 = self.muthap_lst[1] # hap
        maf_lst = [maf1, maf2]

        muthap_setlist_filter = super().mk_setlist(maf_lst, ) # 옵션 지정
        muthap_setlist_unfilter = super().mk_setlist(maf_lst, )
    
        self.muthap_gene_list_filter = super().mk_gene_list(muthap_setlist_filter)
        self.muthap_gene_list_unfilter = super().mk_gene_list(muthap_setlist_unfilter)


    
    # def _mk_muthap_gene_list(self):

    #     for i in range(len(self.muthap_lst)):
    #         if i%2 == 0:
    #             maf1 = self.muthap_lst[i] # mut
    #             maf2 = self.muthap_lst[i+1] # hap
    #             maf_lst = [maf1, maf2]

    #             muthap_setlist_lst_filter = super().mk_setlist(maf_lst)
    #             muthap_setlist_lst_unfilter = super().mk_setlist(maf_lst)
                
    #             self.muthap_setlist_lst_filter.append(muthap_setlist_lst_filter)
    #             self.muthap_setlist_lst_unfilter.append(muthap_setlist_lst_unfilter)
        
    #     for i in range(len(muthap_setlist_lst_filter)):
    #         assert len(muthap_setlist_lst_filter) == len(muthap_setlist_lst_unfilter)

    #         muthap_gene_list_filter = super().mk_gene_list(muthap_setlist_lst_filter[i])
    #         muthap_gene_list_unfilter = super().mk_gene_list(muthap_setlist_lst_unfilter[i])

    #         self.muthap_gene_list_lst_filter.append(muthap_gene_list_filter)
    #         self.muthap_gene_list_lst_unfilter.append(muthap_gene_list_unfilter)
            

    def add_DPtag_unfilter():
        
        # filter_lst = ['t_depth', 't_ref_count', 't_alt_count']
        filter_dict = {'t_depth':[30, ';dp30'], 't_alt_count':[5, ';v5']}
        # print(filter_dict.keys())

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



    def _mk_comp_gene_list(self, _draw_venn=False):
        muthap_setlist_lst = super().mk_setlist()

        if _draw_venn:
            super().draw_venn_diagram(muthap_setlist_lst)
        
        self.muthap_gene_list_lst = super().mk_gene_list(muthap_setlist_lst)





 
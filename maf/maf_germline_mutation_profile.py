from ast import Try
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from jun_tools import jun_mtd as jm  # pip install YjmTool
import os
from typing import List, Tuple, Dict, Union
from typing_extensions import Final


"""
germline mutation profile plotting

sample pair.csv를 받음: 

# 첫 열 = Normal(origin) / 다음 열: diff_sample - 여러개일시 구분자로 ':' 사용 / 한번 더 분화했다면 열 추가
# pair 정보의 동일 행의 데이터들은 같은 그룹으로 묶여서 플로팅


"""


class mutation_profiling():

    def __init__(self, input_path_list: List[str], pair_info_csv_path: str) -> None:
        self.INPUT_PATH_LIST: Final = input_path_list
        self.INPUT_ID_LIST: Final = self.make_sample_id_list(
            self.INPUT_PATH_LIST)
        print(self.INPUT_PATH_LIST[0])
        print(self.INPUT_ID_LIST[0])

        self.mut_count_df = pd.DataFrame(
            columns=['mutation_count', 'sample_group'])

        self.pair_path_dict = self.read_pair_info(pair_info_csv_path)

        print(self.mut_count_df)
        self.plotting(self.mut_count_df)

    def make_sample_id_list(self, _maf_list: List[str]) -> List[str]:
        id_list = []
        for i in range(len(_maf_list)):
            maf_df = pd.read_csv(_maf_list[i], sep='\t')
            id_list.append(maf_df.loc[1, 'Tumor_Sample_Barcode'])
        return id_list

    def check_isin_input_list(self, _target_sample: str) -> List[Union[bool, int]]:
        try:
            _target_idx = self.INPUT_ID_LIST.index(_target_sample)
        except ValueError:
            _target_idx = None

        return [_target_sample in self.INPUT_ID_LIST, _target_idx]

    def make_mutation_count_df(self, _sample_name, _target_info, _grp) -> None:
        if _target_info[0]:
            target_maf_path = self.INPUT_PATH_LIST[_target_info[1]]
        else:
            target_maf_path = None

        print(target_maf_path)
        try:
            target_maf_df = pd.read_csv(target_maf_path, sep='\t')
            mut_count = target_maf_df.shape[0]

        except ValueError:
            mut_count = 0

        self.mut_count_df.loc[_sample_name] = [mut_count, _grp]

    def read_pair_info(self, _pair_info_path: str) -> Dict[str, int]:
        pair_df = pd.read_csv(_pair_info_path)
        for index, rows in pair_df.iterrows():
            print(index)
            print(rows)
            print(len(rows))
            for i in range(len(rows)):
                if i == 0:
                    origin_name = rows[i]
                samples = rows[i].split(':')
                for sample in samples:
                    print(sample)
                    target_info: List[Union[bool, int]
                                      ] = self.check_isin_input_list(sample)
                    self.make_mutation_count_df(
                        sample, target_info, origin_name)

    def plotting(self, _df) -> None:
        sns.set(font_scale=1)

        sns.barplot(data=_df,
                    x=_df.index,
                    y="mutation_count",
                    hue="sample_group",
                    errwidth=20)

        plt.legend(loc=2, prop={'size': 10})
        plt.xticks(rotation=-90)
        plt.show()


input_dir = r'E:/stemcell/VQSR_MAF/DP_AF_filtered_maf/exonic_maf/'
input_format = r"*.maf"

input_pair_info = r'E:/stemcell/stemcell_pair_wide_220126.csv'

input_lst = jm.get_input_path_list(input_dir, input_format, False)
output_dir_name = 'mut_profile'
output_dir = jm.set_output_dir(input_dir, output_dir_name)

mut_prof = mutation_profiling(input_lst, input_pair_info)


mut_prof.mut_count_df.to_csv(os.path.join(output_dir, 'WES65samples_germline.csv'),\
    header=True, index=True)

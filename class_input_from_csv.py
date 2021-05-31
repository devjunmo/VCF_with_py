import pandas as pd
import os
import glob

class MakePairInputList():

    def __init__(self, _pair_info_path, _path_vcf1, _path_vcf2, _enum):
        self.enum_data = _enum
        self.pair_info_df = pd.read_csv(_pair_info_path)
        self.path_vcf1 = pd.DataFrame(glob.glob(_path_vcf1))
        self.path_vcf2 = pd.DataFrame(glob.glob(_path_vcf2))
        self.enum_comp = []
        for e in self.enum_data:
            self.enum_comp.append(e.name)

    def trim_pair_df(self):
        # print(self.enum_comp)
        for col_name in self.pair_info_df.columns:
            # print(col_name)
            if col_name not in self.enum_comp:
                self.pair_info_df.drop(columns=[f"{col_name}"], inplace=True)
        
        # print(self.pair_info_df)

        MakePairInputList.__rm_list_not_exist_sample(self.pair_info_df, self.enum_data, self.path_vcf1)
        MakePairInputList.__rm_list_not_exist_sample(self.pair_info_df, self.enum_data, self.path_vcf2)


    @staticmethod
    def __rm_list_not_exist_sample(_pair_info_df, _enum_data, _path_df):
        for i in range(len(_enum_data)):
            for rows in _pair_info_df.itertuples():
                _count = _path_df[0].str.contains(rows[i+1])
                # print(rows)
                # print(_path_df)
                # print(rows[i+1])
                # print("카운트 =", _count)
                MakePairInputList.__count_ctrl(_pair_info_df, _pair_info_df.columns[i], rows[i+1], _count.sum())


    @staticmethod
    def __count_ctrl(_df, _col_name, rm_content, _data_count):
        if _data_count == 0:
            print('없는 페어 데이터 목록 삭제')
            rm_idx = _df[_df[_col_name] == rm_content].index
            _df.drop(rm_idx[0], inplace=True)

        elif _data_count > 1:
            print('데이터 무결성 검정')
            exit(1)



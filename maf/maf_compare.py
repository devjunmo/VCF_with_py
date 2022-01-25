from msilib.schema import Class
import seaborn as sns
import os
from random import sample
from turtle import shape
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
import venn
import itertools
from itertools import chain, combinations
import matplotlib.pyplot as plt

# sys.path.append(os.path.join(os.getcwd(), 'jun_tools'))
from jun_tools import jun_mtd as jm  # pip install YjmTools

print(sys.path)


input_dir = r"E:/stemcell/somatic_analysis/maf/mutect2"
input_format = r"*.maf"
output_dir_name = r"comp_result"
output_dir = jm.set_output_dir(input_dir, output_dir_name)
pair_info_path = ""


input_lst = jm.get_input_path_list(input_dir, input_format, False)

dic_var: Dict[str, int] = {"hello": 47}


class MafCompare:
    def __init__(
        self, _input_mafs_paths: List, _input_maf: str, _pair_info_path: str
    ) -> None:

        self.maf_list = _input_maf_list
        self.pair_info: Dict[str, str] = self.__mk_pair_dict(_pair_info_path)

    def __mk_pair_dict(self) -> dict:
        pass  # sample name에 pre/suf-fix 붙여서 path dict 형성

    def draw_venn_diagram(self, _save_path: str, _sample_num: int, _is_show=False):
        pass
    
    def powerset(_iterable):
        s = list(_iterable)
        return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


from msilib.schema import Class
import seaborn as sns
import os
from random import sample
from turtle import shape
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
# sys.path.append(os.path.join(os.getcwd(), 'jun_tools'))
from jun_tools import jun_mtd as jm # pip install YjmTools

print(sys.path)


input_dir = r'E:/stemcell/somatic_analysis/maf/mutect2'
input_format = r'*.maf'
output_dir_name = r'comp_result'
output_dir = jm.set_output_dir(input_dir, output_dir_name)

input_lst = jm.get_input_path_list(input_dir, input_format, False)


class MafCompare():
    
    def __init__(self, _input_maf_list) -> None:
        self.maf_list = _input_maf_list
    
    
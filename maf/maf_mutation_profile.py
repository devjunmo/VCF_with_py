import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from jun_tools import jun_mtd as jm  # pip install YjmTool


input_dir = r""
input_format = r"*.maf"
input_lst = jm.get_input_path_list(input_dir, input_format, False)

output_dir_name = r"comp_result"
output_dir = jm.set_output_dir(input_dir, output_dir_name)







# 디렉토리의 마프파일 가져오기
# {Tumor_id: [N_id, row_count]} 딕셔너리 만들기


input_dir = r''

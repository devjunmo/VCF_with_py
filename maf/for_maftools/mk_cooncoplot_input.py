"""
R에서 co-oncoplot을 사용할때 필요한 maf 생성

1. maf rbind
2. 필요 컬럼 추출
3. cohort 정보 추가

"""

import pandas as pd
from jun_tools import jun_mtd as jm  # pip install YjmTool
import os
from typing import List, Tuple, Dict, Union
from typing_extensions import Final


# input_dir = r'E:/stemcell/somatic_analysis/maf/mutect2/DP_AF_filtered_maf/exonic_maf/diff_comp'
input_dir =r'E:/stemcell/somatic_analysis/maf/mutect2/DP_AF_filtered_maf/Diff_whole_regions'
input_format = r"*.maf"

cohort_info_path = r'E:/stemcell/somatic_analysis/maf/mutect2/DP_AF_filtered_maf/exonic_maf/diff_comp/diff_cohort_info.csv'

input_maf_lst = jm.get_input_path_list(input_dir, input_format, False)

output_dir_name = 'cohort_comp'
output_dir = jm.set_output_dir(input_dir, output_dir_name)
output_maf_name = 'diff_cohort_for_oncoplot.maf'

# print(input_lst)


def join_with_info(_maf_df, _cohort_info_df):
    joined_df = pd.merge(_maf_df, _cohort_info_df, on='Tumor_Sample_Barcode')
    return joined_df


def rbind_mafs(_target_maf, _curr_maf):
    if _target_maf is None:
        _target_maf = _curr_maf
        print(f'stacked rows: {_target_maf.shape}')
    else:
        _target_maf = pd.concat([_target_maf, _curr_maf])
        print(f'stacked rows: {_target_maf.shape}')

    return _target_maf


def select_maf_col(_raw_maf_df):
    select_column = ['Hugo_Symbol', 'Center', 'NCBI_Build',
                     'Chromosome', 'Start_Position', 'End_Position', 'Variant_Classification',
                     'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'HGVSp_Short',
                     'Tumor_Sample_Barcode']

    rename_column = ['Hugo_Symbol', 'Center', 'NCBI_Build',
                     'Chromosome', 'Start_Position', 'End_Position', 'Variant_Classification',
                     'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2', 'amino_acid_change',
                     'Tumor_Sample_Barcode']

    _raw_maf_df = _raw_maf_df[select_column]
    _raw_maf_df.columns = rename_column

    return _raw_maf_df


def make_input_maf(_maf_list, _cohort_df):
    cohort_lst = _cohort_df.iloc[:, 0].tolist()
    rbind_maf_df = None

    for i in range(len(_maf_list)):
        maf_df = pd.read_csv(_maf_list[i], sep='\t')
        try:
            sample_name = maf_df.loc[0, 'Tumor_Sample_Barcode']
            if sample_name in cohort_lst:
                print(sample_name)
            else:
                continue
        except KeyError:
            print(
                '\033[32m' + f'-- There is no mutation in this sample. --' + '\033[0m')
            print('\033[33m' + f'{_maf_list[i]}' + '\033[0m')
            print(
                '\033[32m' + f'------------------------------------------' + '\033[0m')
            continue

        maf_df = select_maf_col(maf_df)
        rbind_maf_df = rbind_mafs(rbind_maf_df, maf_df)

    rbind_maf_df.reset_index(inplace=True, drop=True)
    processed_maf_df = join_with_info(rbind_maf_df, _cohort_df)
    print(processed_maf_df)

    return processed_maf_df


COHORT_DF: Final = pd.read_csv(cohort_info_path)
print(COHORT_DF)

input_maf = make_input_maf(input_maf_lst, COHORT_DF)

output_path = os.path.join(output_dir, output_maf_name)

input_maf.to_csv(output_path, index=False, header=True, sep='\t')


# T only파일에 bad GT를 제거하지 않았어서 작성한 script 파일 # 

# 확인해 보니 Tonly 만들때 GT 정리한 데이터로 만들었었음!! left join시 NaN 존재 X

import pandas as pd
from glob import glob
# import natsort


# Tonly
tonly_bed_dir = r'/data_244/WES/T-only_bed/'

# Tsp, DP30 적용
dp_applied_tsp_processed_bed_dir= r'/data_244/WES/T_DP_O_NDP_samples/pass_only_Teratoma_specifics/teratoma_specific_processed_bed_210606/'


t_only_bed_lst = glob(tonly_bed_dir + 'SNP*') 
t_sp_bed_lst = glob(dp_applied_tsp_processed_bed_dir + 'SNP*')
# print(t_only_bed_lst)


t_only_header = ['CHROM', 'START', 'END']
t_sp_header = ['CHROM', 'START', 'END', 'REF', 'ALT']


for i in range(len(t_sp_bed_lst)):

    t_only_bed_df = pd.read_csv(t_only_bed_lst[i], names=t_only_header, sep='\t')
    # t_only_bed_df['CHROM'] = t_only_bed_df['CHROM'].apply(lambda _: str(_))
    # t_only_bed_df['START'] = t_only_bed_df['START'].apply(lambda _: str(_))
    # t_only_bed_df['END'] = t_only_bed_df['END'].apply(lambda _: str(_))
    # dp_aply_pced_bed_df['pk'] = dp_aply_pced_bed_df['CHROM'] + '_' + dp_aply_pced_bed_df['START'] + '_' + dp_aply_pced_bed_df['END']
    print(t_only_bed_df.head())
    print(t_only_bed_df.shape)


    t_sp_bed_df = pd.read_csv(t_sp_bed_lst[i], names=t_sp_header, sep='\t')
    # t_sp_bed_df['CHROM'] = t_sp_bed_df['CHROM'].apply(lambda _: str(_))
    # t_sp_bed_df['START'] = t_sp_bed_df['START'].apply(lambda _: str(_))
    # t_sp_bed_df['END'] = t_sp_bed_df['END'].apply(lambda _: str(_))
    # t_sp_bed_df['pk'] = t_sp_bed_df['CHROM'] + '_' + t_sp_bed_df['START'] + '_' + t_sp_bed_df['END']
    print(t_sp_bed_df.head())
    print(t_sp_bed_df.shape)


    # left join할껀데, 내 예상 = Tonly랑 컬럼수는 똑같고(포함관계라서), NaN이 섞인곳이 존재(Tonly에 badGT부분)
    res = pd.merge(t_only_bed_df, t_sp_bed_df, how='left', \
                        on = ['CHROM', 'START', 'END'])

    print(res.shape)

    print(res[res['END'].isnull()])
    
    # nan_res = res.isnull()
    # print(nan_res.head())
    # print(nan_res.shape)


    
    
    # print(t_sp_bed_df.shape)
    # print(t_only_bed_df.shape)

    # left_outer_join = pd.merge(t_sp_bed_df, t_only_bed_df, on=['CHROM', 'START', 'END'], how='right')
    # print(left_outer_join.head)
    # print(left_outer_join.shape)
    # left_outer_join.dropna(axis=0, inplace=True)
    # print(left_outer_join.shape)

    break
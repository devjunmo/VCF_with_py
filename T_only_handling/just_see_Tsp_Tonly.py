
import pandas as pd
from glob import glob

from pandas.core.tools.numeric import to_numeric
# import natsort

t_sp_bedfile = r'/home/jun9485/data/WES/T_DP_O_NDP_samples/pass_only_Teratoma_specifics/teratoma_specific_processed_bed_210606/'
tonly_bed_dir = r'/home/jun9485/data/WES/T-only_bed/'

t_sp_bed_lst = glob(t_sp_bedfile + 'SNP*') 
tonly_bed_lst = glob(tonly_bed_dir + '*')


bed_header = ['CHROM', 'START', 'END']

show_df = pd.DataFrame(columns=['Tsp', 'Tonly'])

for i in range(len(t_sp_bed_lst)):

    print(t_sp_bed_lst[i])
    print(tonly_bed_lst[i])
    s_name = t_sp_bed_lst[i].split(r'/')[-1].split('.')[0].split('_')[1]

    t_sp_bed_df = pd.read_csv(t_sp_bed_lst[i], names=bed_header, sep='\t')
    print(t_sp_bed_df.head())

    tonly_bed_df = pd.read_csv(tonly_bed_lst[i], names=bed_header, sep='\t')
    print(tonly_bed_df.head())
    
    print(t_sp_bed_df.shape)
    print(tonly_bed_df.shape)
    a = str(t_sp_bed_df.shape[0])
    b = str(tonly_bed_df.shape[0])
    # print(type(a))

    show_df.loc[s_name] = [a, b]

    
    
    

    # break

show_df = show_df.apply(pd.to_numeric, axis=1)
show_df["%"] = round((show_df["Tonly"] / show_df["Tsp"]) * 100, 2)
print(show_df)
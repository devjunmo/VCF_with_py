import pandas as pd
from pandas.core.algorithms import unique


pd.set_option('display.max_rows', None)


maf_path = r'E:/stemcell_ips/somatic_call/vardict/hiPS65/hiPS65-A-p10_rmHd.maf'

maf_df = pd.read_csv(maf_path, sep='\t')

# print(maf_df.columns)

print(maf_df.shape) # (137918, 118)

print(pd.value_counts(maf_df['Variant_Classification'])) # 

# print(unique(maf_df['FILTER']))
# print(pd.value_counts(maf_df['FILTER']))


# pass만 남기기

non_pass_idx = maf_df[maf_df['FILTER'] != 'PASS'].index
maf_df = maf_df.drop(non_pass_idx)
maf_df.reset_index(inplace=True, drop=True)

print(maf_df.shape) # (17572, 118)

print(pd.value_counts(maf_df['Variant_Classification']))
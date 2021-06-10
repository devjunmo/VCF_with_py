

import pandas as pd


df = pd.DataFrame()

df['test'] = [1, 2, 3, 4]

print(df)


out = '/home/jun9485/data/WES/test.csv'

df.to_csv(out)
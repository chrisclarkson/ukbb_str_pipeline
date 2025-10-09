import pandas as pd
import numpy as np
import os
os.system('rm tmp.csv')
data_full=pd.read_csv('ukbb_informal_dated_main_icd10_codes.tsv',sep='\t')
icd10_codes=pd.read_csv('icds_gt36.tsv',sep='\t',header=None)

for data in np.array_split(data_full,1000):
    print(data.shape)
    out=pd.DataFrame(index=data['eid'],columns=icd10_codes[0])
    for idx, row in data.iterrows():          # idx is the index label, row is a Series
        if pd.isna(row['codes']):
            continue
        codes = [c.strip() for c in str(row['codes']).split(',')]
        for j, code in enumerate(codes):
            if code in out.columns:
                # option A: positional in the row (works if the j+2 logic is correct)
                val = row.iloc[j + 2]

                # option B (safer): use the actual column name for j+2
                # colname = data.columns[j + 2]
                # val = row[colname]

                out.at[row['eid'], code] = val   # .at is a bit faster for scalar writes

    if os.path.exists('tmp.csv'):
        out.to_csv('tmp.csv',mode='a',header=False)
    else:
        out.to_csv('tmp.csv')


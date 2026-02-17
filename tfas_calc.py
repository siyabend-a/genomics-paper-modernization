# TEST
import pandas as pd
import numpy as np

# 1: TSS start, 2: NCBI Reference Sequence, 13: Peak summit, 16: Intensity(g_k)
df = pd.read_csv("output/pre_TFAS.txt", sep="\t", header=None, usecols=[1, 3, 13, 16]) 
df.columns = ['tss_start', 'ref_seq', 'peak_summit', 'intensity']
d = np.abs(df['tss_start'] - df['peak_summit'])

# g_k*e^(-d/d_0) exponential decay as provided in the study.
# d_0 = 500bps for E2f1, 5000bps for others. (Ouyang et al., 2009)
df['tfas'] = df['intensity']*np.exp(-d/500)
tfas = df.groupby(['ref_seq'])['tfas'].sum().reset_index()

# Normalize the scores
col_means = tfas['tfas'].mean()
col_stds = tfas['tfas'].std()
tfas['normalized_tfas'] = (tfas['tfas']-col_means)/col_stds

print(tfas['normalized_tfas'])
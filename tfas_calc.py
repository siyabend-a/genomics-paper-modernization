# Script for TFAS and binary scores calculation
import pandas as pd
import numpy as np
from functools import reduce
import qnorm
from sklearn.preprocessing import StandardScaler

dfs_tfas=[]
dfs_bin=[]
tfs=("Oct4", "Sox2", "Nanog", "Esrrb", "Kfl4", "Tcfcp2l1", "Stat3", "Smad1", "Myc", "Mycn", "Zfx", "E2f1")
for tf in tfs:
	# 1: TSS start, 2: NCBI Reference Sequence, 13: Peak summit, 16: Intensity(g_k)
	df = pd.read_csv(f"output/{tf}_pretfas.txt", sep="\t", header=None, usecols=[1, 2, 3, 5, 13, 18, 21]) 
	df.columns = ['tss_start_1', 'tss_start_2', 'gene', 'strand', 'peak_summit', 'intensity', 'rel_summit']
	df['tss_start_1'] = np.where(df['strand'] == '+', df['tss_start_1'], df['tss_start_2'])
	df['peak_summit'] = df['peak_summit'] + df['rel_summit']
	d = np.abs(df['tss_start_1'] - df['peak_summit'])
	# TFAS calculation
	# g_k*e^(-d/d_0) exponential decay as provided in the study.
	# d_0 = 500bps for E2f1, 5000bps for others. (Ouyang et al., 2009)
	df['tfas'] = df['intensity']*np.exp(-d/500) if tf=="E2f1" else df['intensity']*np.exp(-d/5000)
	tfas = df.groupby(['gene'])['tfas'].sum().reset_index()
	
	# Normalize the scores
	#col_means = tfas['tfas'].mean()
	#col_stds = tfas['tfas'].std()
	#tfas[f'{tf}'] = ((tfas['tfas']-col_means)/col_stds).fillna(0)
	tfas[f'{tf}'] = tfas['tfas']# np.log(tfas['tfas'] + 0.1)
	dfs_tfas.append(tfas[['gene', f'{tf}']])

	# Binary calculation
	#df[f'b{tf}'] = (d<=2000).astype(int) # 2000bps distance to TSS
	#bin_score = df.groupby('ref_seq')[f'b{tf}'].max().reset_index() # Unique transcripts
	#dfs_bin.append(bin_score[f'b{tf}'])

#final_tfas = pd.concat(dfs_tfas, axis=1)

#dfs_tfas = dfs_tfas.dropna()
scaler = StandardScaler()
final_df = reduce(lambda left, right: pd.merge(left, right, on='gene', how='outer'), dfs_tfas).fillna(0)
final_df[list(tfs)] = np.log(final_df[list(tfs)] + 0.1)
#final_df[list(tfs)] = scaler.fit_transform(final_df[list(tfs)])
#final_df[list(tfs)] = qnorm.quantile_normalize(final_df[list(tfs)])
#print(dfs_tfas)

#rpkm = pd.read_csv(f"sd1.txt", sep="\t").drop(['Entrez_ID', 'EB_RPKM'], axis=1)
#rpkm['ESC_RPKM'] = np.log(rpkm[['ESC_RPKM']] + 0.1)
#rpkm['ESC_RPKM'] = qnorm.quantile_normalize(rpkm[['ESC_RPKM']])
#col_means = rpkm['rpkm'].mean()
#col_stds = rpkm['rpkm'].std()
#rpkm['rpkm'] = ((rpkm['rpkm']-col_means)/col_stds).fillna(0)

rpkm = pd.read_csv(f"final_rpkm.csv", sep="\t")
#rpkm['rpkm'] = scaler.fit_transform(rpkm[['rpkm']])
#sd1=pd.read_csv('sd1.txt', sep='\t')

final_df=final_df.merge(rpkm, on='gene', how='inner')
#final_df=final_df.merge(sd1, on='gene', how='inner').drop(['ESC_RPKM', 'EB_RPKM', 'Entrez_ID'], axis=1)


#final_df.to_csv("total.txt", index=False, sep='\t')

final_df.to_csv(f"total2.txt", index=False, sep='\t')
#pd.concat(final_df, axis=1).to_csv(f"TFAS.txt", index=False, sep='\t')
#pd.concat(dfs_bin, axis=1).to_csv(f"BinScore.txt", index=False)

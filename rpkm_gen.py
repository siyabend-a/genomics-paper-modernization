import pandas as pd

df = pd.read_csv("counts.txt", sep="\t", skiprows=1)

# RPKM calculation
nm = df['merged_rnaseq.bam'].sum()/1000000
df['rpkm']=(df['merged_rnaseq.bam']/nm)/(df['Length']/1000)

df[['gene', 'rpkm']].to_csv('final_rpkm.csv', sep='\t', index=False)

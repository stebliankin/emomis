import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
import os

np.random.seed(7272)

if not os.path.exists('lists'):
    os.mkdir('lists')

out_list = 'lists/all_ppi_nonredundant.txt'
summary_updated = 'metadata/summary_updated.tsv'

df = pd.read_csv('metadata/20210511_0686847_summary.tsv', sep='\t')
#df = df[["pdb", "Hchain", "Lchain", "antigen_chain"]]


df["Hchain"] = df["Hchain"].fillna("")
df["Lchain"] = df["Lchain"].fillna("")

df = df.dropna(subset=["pdb", "Hchain", "Lchain", "antigen_chain"]).reset_index(drop=True)
df["antigen_chain"] = df["antigen_chain"].apply(lambda x: x.replace(" | ", ""))

df['Hchain'] = df.apply(lambda row: row.Hchain if row.Hchain!=row.Lchain else '', axis=1)

df["PPI"] = df.apply(lambda x: "{}_{}_{}".format(x.pdb.upper(), x.Hchain+x.Lchain, x.antigen_chain), axis=1)

df = df.drop_duplicates(subset=["pdb", "Hchain", "Lchain", "antigen_chain"]).reset_index(drop=True)


unique_ppi = df["PPI"].unique()

df_ppi = pd.DataFrame({'PPI':unique_ppi})
df_ppi = df_ppi.drop_duplicates()
print(df_ppi)

df_ppi.to_csv(out_list, header=False, index=False)
df.to_csv(summary_updated, sep='\t', index=False)


# Allocate SARS2-complexes as test

# Split training and testing
all_pid = df["pdb"].unique()
train, test = train_test_split(all_pid, test_size=0.2)
print("Number of unique complexes: {}".format(len(all_pid)))

# train_df = df[df["pdb"].isin(train)]
# train_df["PPI"].to_csv("lists/train_ppi.txt",header=False, index=False)
# print(len(train_df))
#
# test_df = df[df["pdb"].isin(test)]
# test_df["PPI"].to_csv("lists/test_ppi.txt",header=False, index=False)
# print(len(test_df))

############
# STDOUT:
# Number of unique complexes: 536
############

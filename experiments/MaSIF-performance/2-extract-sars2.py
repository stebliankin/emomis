import pandas as pd

out_sars = 'lists/sars2_nonredundant.txt'

all_ppi = [x.strip('\n') for x in open('lists/all_ppi_nonredundant.txt').readlines()]
metadata_file = 'metadata/summary_updated.tsv'
keyword = 'severe acute respiratory syndrome coronavirus2'

meta_df = pd.read_csv(metadata_file, sep='\t')

meta_df["Hchain"] = meta_df["Hchain"].fillna("")
meta_df["Lchain"] = meta_df["Lchain"].fillna("")

meta_df["antigen_chain"] = meta_df["antigen_chain"].apply(lambda x: x.replace(" | ", ""))


meta_df = meta_df[['pdb', 'Hchain', 'Lchain', 'antigen_chain', "antigen_species"]]
meta_df = meta_df.dropna()
meta_df = meta_df[meta_df['antigen_species'].str.contains(keyword)]

meta_df["PPI"] = meta_df.apply(lambda x: "{}_{}_{}".format(x.pdb.upper(), x.Hchain+x.Lchain, x.antigen_chain), axis=1)

# meta_df = meta_df[meta_df['PPI'].isin(all_ppi)]
#
# print(meta_df)

print('Number of unique SARS2: {}'.format(len(meta_df['pdb'].unique())))
meta_df = meta_df['PPI']
meta_df=meta_df.drop_duplicates()

meta_df.to_csv(out_sars, index=False, header=False)

# STDOUT:
# Number of unique SARS2: 19



import pandas as pd

df = pd.read_csv('output/variants/4-MaSIF_scores.tsv', sep='\t')
df_unique_indx = df.sort_values('score_target_subject', ascending=True).drop_duplicates(['query_target','match','subject',
                                                                                         'PDB_target','PDB_subject'])

metadata_df = pd.read_csv('../demo-12-2021/metadata/sabdab_summary_all.tsv', sep='\t')
metadata_df['PDB_subject'] = metadata_df['pdb'].apply(lambda x: x.upper())
metadata_df = metadata_df[['PDB_subject', 'antigen_name', 'organism', 'authors']]
df_unique_indx = df_unique_indx.merge(metadata_df, on='PDB_subject', how='left')

df_unique_indx = df_unique_indx.drop(['eval','filter1_flag','filter2_flag','filter3_flag'], axis=1)
df_unique_indx = df_unique_indx.drop_duplicates()
df_unique_indx.to_csv('output/variants/5-MaSIF-unique.tsv', index=False, sep='\t')

print(df)


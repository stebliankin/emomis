import pandas as pd

df = pd.read_csv('output/4-MaSIF_scores.tsv', sep='\t')
df_unique_indx = df.sort_values('score_target_subject', ascending=True).drop_duplicates(['query_target','match','subject'])

metadata_df = pd.read_csv('metadata/sabdab_summary_all.tsv', sep='\t')
metadata_df['PDB_target'] = metadata_df['pdb'].apply(lambda x: x.upper())
metadata_df = metadata_df[['PDB_target', 'antigen_name', 'organism', 'authors']]
df_unique_indx = df_unique_indx.merge(metadata_df, on='PDB_target', how='left')

df_unique_indx = df_unique_indx.drop(['eval','filter1_flag','filter2_flag','filter3_flag'], axis=1)
df_unique_indx = df_unique_indx.drop_duplicates(['query_target','match','subject'])
df_unique_indx.to_csv('output/5-MaSIF-unique.tsv', index=False, sep='\t')

print(df)

spike_list = list(set([x.strip('\n').split('_')[0] for x in open('lists/target.txt', 'r').readlines()]))
db_list = list(set([x.strip('\n').split('_')[0] for x in open('lists/db.txt', 'r').readlines()]))

print(len(spike_list))
print(', '.join(spike_list))
print()

print(len(db_list))
print(', '.join(db_list))
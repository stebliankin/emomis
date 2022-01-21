#from sklearn.model_selection import train_test_split
import random

random.seed(7272)
all_pids = set([x.strip('\n') for x in open('lists/processed_pids.txt', 'r').readlines()])
all_ppis = [x.strip('\n') for x in open('lists/processed_ppis.txt', 'r').readlines()]

# PIDs to exclude from testing:
intersection_pids = set([x.strip('\n') for x in open('lists/intersection_masif_sabdab.txt', 'r').readlines()])
homolog_pids = set([x.strip('\n') for x in open('lists/homologs_masif_sabdab.txt', 'r').readlines()])

# PIDs to exclude from training and testing
sars2_pids = set([x.strip('\n').split('_')[0] for x in open('lists/sars2_nonredundant.txt', 'r').readlines()])

print("Total Number of PDBs: {}".format(len(all_pids)))
print("Number of complexes that intersects with MaSIF train: {}".format(len(intersection_pids)))
print("Number of homolog complexes with MaSIF train: {}".format(len(homolog_pids)))
print("Number of SARS2 complexes: {}".format(len(sars2_pids)))

pids_to_select_test = all_pids.difference(sars2_pids)
pids_to_select_test = pids_to_select_test.difference(intersection_pids)
pids_to_select_test = pids_to_select_test.difference(homolog_pids)

print("Number of non-redundant complexes for test: {}".format(len(pids_to_select_test)))
print("Number of SARS2 PIDs: {}".format(len(sars2_pids)))
#print(homolog_pids.intersection(sars2_pids))

with open('lists/train.txt', 'w') as out:
    out.write('')

with open('lists/test_pids.txt', 'w') as out:
    for pid in pids_to_select_test:
        out.write(pid+'\n')

with open('lists/sars2_pids.txt', 'w') as out:
    for pid in sars2_pids:
        out.write(pid+'\n')

with open('lists/test.txt', 'w') as out:
    for ppi in all_ppis:
        pid = ppi.split('_')[0]
        if pid in pids_to_select_test:
            out.write(ppi+'\n')


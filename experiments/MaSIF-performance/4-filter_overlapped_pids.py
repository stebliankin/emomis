# Here we will select PDB IDs of the SAbDab that is not overlapped with MaSIF train


sabdab_pid_list = set([x.strip('\n') for x in open('lists/processed_pids.txt', 'r').readlines()])
masif_pid_list = set([x.strip('\n').split('_')[0] for x in open('lists/masif_original_train.txt', 'r').readlines()])

pid_overlapped_file = 'lists/intersection_masif_sabdab.txt'
pid_unique_file = 'lists/unique_sabdab.txt'

pid_intersection = sabdab_pid_list.intersection(masif_pid_list)
pid_unique = sabdab_pid_list.difference(masif_pid_list)

print('Length of intersection of SAbDab and MaSIF: {}'.format(len(pid_intersection)))
print('Length of unique SAbDab complexes: {}'.format(len(pid_unique)))

with open(pid_overlapped_file, 'w') as out:
    for pid in pid_intersection:
        out.write(pid+'\n')

with open(pid_unique_file, 'w') as out:
    for pid in pid_unique:
        out.write(pid+'\n')
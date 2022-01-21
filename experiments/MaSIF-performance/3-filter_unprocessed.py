# Select complexes that we were able to process

cached_list = '../04-2021-SAbDab/lists/cached_standard-sc05-min3-nomax.txt' # the list of PPIs from which we extracted patches
all_ppis = 'lists/all_ppi_nonredundant.txt'
out_ppis = 'lists/processed_ppis.txt'
out_pids = 'lists/processed_pids.txt'

cached_ppis = [x.strip('\n') for x in open(cached_list, 'r').readlines()]
current_ppis = [x.strip('\n') for x in open(all_ppis, 'r').readlines()]
processed_pids = []

with open(out_ppis, 'w') as out:
    for ppi in current_ppis:
        if ppi in cached_ppis:
            out.write(ppi+'\n')
            pid = ppi.split('_')[0]
            if pid not in processed_pids:
                processed_pids.append(pid)

with open(out_pids, 'w') as out:
    for pid in processed_pids:
        out.write(pid+'\n')



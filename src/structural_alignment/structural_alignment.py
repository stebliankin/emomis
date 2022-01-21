import pandas as pd
import os
import pdb
import subprocess
from tqdm import tqdm
from scipy.stats import zscore
from scipy.stats import norm
def convert_threeAA_oneAA(one_code):
    # convert one amino acid code to three
    aa_dict = {'CYS':'C', 'ASP':'D', 'SER': 'S', 'GLN':'Q', 'LYS':'K',
         'ILE':'I', 'PRO':'P', 'THR':'T', 'PHE':'F', 'ASN':'N',
         'GLY':'G', 'HIS':'H', 'LEU':'L', 'ARG':'R', 'TRP':'W',
         'ALA':'A', 'VAL':'V', 'GLU':'E', 'TYR':'Y', 'MET':'M'}
    return aa_dict[one_code]


def extract_motif_pdb(pid, ch, start_res, end_res, motif, config):
    raw_pdb_dir = config['dirs']['raw_pdb']
    pdb_file = raw_pdb_dir + 'pdb{}.ent'.format(pid.lower())

    out_motif_file = config['dirs']['motifs_pdb'] + '{}_{}.pdb'.format(pid, motif)
    sequence_str = ''

    prev_resid_i = 0
    with open(out_motif_file, 'w') as out:
        with open(pdb_file, 'r') as f:
            for line in f.readlines():
                if line[:4]=='ATOM':
                    resid_i = int(line[22:26])
                    ch_i = line[21]
                    resid_i = int(line[22:26])
                    resname_i = line[17:20].strip(' ')
                    if ch_i==ch and resid_i>=start_res and resid_i<end_res:
                        out.write(line)
                        if resid_i!=prev_resid_i:
                            res_letter = convert_threeAA_oneAA(resname_i)
                            sequence_str+=res_letter
                        prev_resid_i=resid_i

    return out_motif_file, sequence_str

def parse_TM_align(stdout):
    rmsd, tm_score = None, None
    for line in str(stdout).split('\n'):
        if "RMSD=" in line:
            rmsd = float(line.split("RMSD=")[1].split(',')[0].strip(' ').strip('\t'))
        if 'TM-score=' in line:
            tm_score = float(line.split('TM-score=')[1].split('(if')[0].strip(' ').strip('\t'))
    return rmsd, tm_score

def save_fasta_file(seq_target, seq_subject, fasta_file):
    with open(fasta_file, 'w') as out:
        out.write('>target\n')
        out.write(seq_target+'\n')
        out.write('>subject\n')
        out.write(seq_subject+'\n')

def compute_exact_match_res(match, hit_i):
    # match example: "E  +R  NITN"
    # hit_i example: "9,10"
    # Output: residues of the exact match - 7-10

    #res_contact_i = int(contact_i.split(',')[0].split(':')[1]) # 277
    res_hit_i = int(hit_i.split(',')[0].split('.')[0])

    # go to the left of the match
    left_i = res_hit_i
    while left_i>=0:
        if match[left_i]==' ' or match[left_i]=='+':
            break
        left_i-=1

    # go to the right:
    right_i =res_hit_i
    while right_i<=len(match)-1:
        if match[right_i]==' ' or match[right_i]=='+':
            break
        right_i+=1

    # account for the last iteration in the loop (+1).
    left_i+=1
    right_i-=1
    # if match=="SNNS  +PT":
    #     pdb.set_trace()
    return '{},{}'.format(left_i, right_i)


def get_all_rmsd(config):
    # output  RMSDs for each motif length
    curr_path_split = os.path.realpath(__file__).split('/')
    curr_path = '/'.join(curr_path_split[0:len(curr_path_split)-1])
    rmsd_df = pd.read_csv(curr_path+'/rmsd_distribution.csv')
    all_rmsd = {}
    for i in range(3, 33):
        tmp_df = rmsd_df[rmsd_df['Len'] == i]
        rmsd_values = list(tmp_df['RMSD'])
        all_rmsd[i] = list(rmsd_values)
    return all_rmsd

def closest(lst, K):
    # find the closest number in a list
    return min(range(len(lst)), key=lambda i: abs(lst[i] - K))

def filter_rmsd(config):
    all_rmsd_dict = get_all_rmsd(config)

    def compute_pval(rmsd, len_i, all_rmsd_dict):
        # Obtain Z-score of the 'rmsd' for a given 'len_i' and convert it to p-value
        all_rmsd_tmp = all_rmsd_dict[len_i].copy()
        all_rmsd_tmp.append(rmsd)
        #all_rmsd_dict_tmp[len_i] = sorted(all_rmsd_dict_tmp[len_i], reverse=True)
        #rmsd_indx = all_rmsd_dict_tmp[len_i].index(rmsd)

        #all_rmsd_tmp = [-x for x in all_rmsd_tmp]

        # compute z-scores
        z_scores = zscore(all_rmsd_tmp)
        pvals = norm.sf(z_scores)

        return 1-pvals[-1]

    def compute_zscore(rmsd, len_i, all_rmsd_dict):
        # Obtain Z-score of the 'rmsd' for a given 'len_i' and convert it to p-value
        all_rmsd_tmp = all_rmsd_dict[len_i].copy()
        all_rmsd_tmp.append(rmsd)
        z_scores = zscore(all_rmsd_tmp)
        return z_scores[-1]


    df = pd.read_csv(config['out_files']['TM_align_all'], sep='\t')
    df['RMSD_pval'] = df.apply(lambda row: compute_pval(row['RMSD_motifs'], row['len_exact_match'],all_rmsd_dict), axis=1)
    df['RMSD_ZSCORE'] = df.apply(lambda row: compute_zscore(row['RMSD_motifs'], row['len_exact_match'],all_rmsd_dict), axis=1)


    df.to_csv(config['out_files']['TM_align_all'], index=False, sep='\t')

    df['PPI_target'] = df.apply(lambda row: get_ppi(row['PDB_target'], row['target_resid_ag'].split(':')[0], config), axis=1)
    df['PPI_subject'] = df.apply(lambda row: get_ppi(row['PDB_subject'], row['contact_subject_ag'].split(':')[0], config), axis=1)
    # threshold by pval=0.1
    df = df[df['RMSD_pval']<=config['p_val_threshold']]
    df.to_csv(config['out_files']['TM_align_filtered'],index=False, sep='\t')
    updated_ppi_list_db = list(df['PPI_subject'].unique())
    updated_ppi_list_target = list(df['PPI_target'].unique())
    return updated_ppi_list_db, updated_ppi_list_target

# def filter_rmsd(config):
#
#     def get_rmsd_dict(p_val):
#         filter_dict = config['rmsd_threshold_dict']
#         rmsd_dict = {}
#         for i in range(len(filter_dict['RMSD_threshold_{}'.format(p_val)])):
#             rmsd_dict[filter_dict['Motif_length'][i]] = filter_dict['RMSD_threshold_{}'.format(p_val)][i]
#         return rmsd_dict
#
#     rmsd_dict_high = get_rmsd_dict('high') # high confidence p_val<0.05
#     rmsd_dict_medium = get_rmsd_dict('medium') # medium confidence p_val<0.15
#
#     def filter_one(rmsd, motif_len, rmsd_dict):
#         if rmsd <= rmsd_dict[motif_len]:
#             return True
#         else:
#             return False
#
#     df = pd.read_csv(config['out_files']['TM_align_all'], sep='\t')
#     df['RMSD_filter_high'] = df.apply(lambda row: filter_one(row['RMSD_motifs'], row['len_exact_match'], rmsd_dict_high), axis=1)
#     df['RMSD_filter_medium'] = df.apply(lambda row: filter_one(row['RMSD_motifs'], row['len_exact_match'], rmsd_dict_medium), axis=1)
#
#     df = df[df['RMSD_filter_high'] | df['RMSD_filter_medium']]
#     df.to_csv(config['out_files']['TM_align_filtered'], index=False, sep='\t')
#
#     df['PPI_target'] = df.apply(lambda row: get_ppi(row['PDB_target'], row['target_resid_ag'].split(':')[0], config), axis=1)
#     df['PPI_subject'] = df.apply(lambda row: get_ppi(row['PDB_subject'], row['contact_subject_ag'].split(':')[0], config), axis=1)
#
#     updated_ppi_list_db = list(df['PPI_subject'].unique())
#     updated_ppi_list_target = list(df['PPI_target'].unique())
#     return updated_ppi_list_db, updated_ppi_list_target

def get_ppi(pid, ch1, config):
    list_db = [x.strip('\n') for x in open(config['db_list'], 'r').readlines()]
    list_target = [x.strip('\n') for x in open(config['target_list'], 'r').readlines()]
    full_list = list_db+list_target
    for ppi in full_list:
        if pid in ppi and ch1 in ppi:
            return ppi
    return None

def structural_alignment(config):
    # Order - 1st is Target sequence and second is the Subject sequence

    mimicry_file = config['out_files']['filtered_blast_hits']

    df = pd.read_csv(mimicry_file, sep='\t')
    all_rmsds = []
    all_TM_scores = []

    df['exact_match_i'] = df.apply(lambda row: compute_exact_match_res(row['match'], row['hit_i']), axis=1)
    df['len_exact_match'] = df['exact_match_i'].apply(lambda row: int(row.split(',')[1]) - int(row.split(',')[0]) + 1)
    df['PDB_subject'] = df['PDB_subject'].apply(lambda x: x.split('_')[0]) # to remove chain ID
    df['hit_i'] = df['hit_i'].apply(lambda x: x.split('.')[0])

    for i, row in tqdm(df.iterrows()):
        pid_target = row['PDB_target']
        pid_subject = row['PDB_subject']

        exact_match_start, exact_match_end = int(row['exact_match_i'].split(',')[0]), int(row['exact_match_i'].split(',')[1])
        target_ch, contact_res_target = row['target_resid_ag'].split(',')[0].split(':')
        start_res_target = int(contact_res_target) - int(row['hit_i'].split(',')[0]) + 1 + exact_match_start
        end_res_target = int(contact_res_target) - int(row['hit_i'].split(',')[0]) + 1 + exact_match_end + 1

        #pdb.set_trace()
        motif_file_target, seq_target = extract_motif_pdb(pid_target, target_ch, start_res_target, end_res_target, row['query_target'], config)

        subject_ch, contact_res_subject = row['contact_subject_ag'].split(',')[0].split(':')
        start_res_subject = int(contact_res_subject) - int(row['hit_i'].split(',')[0]) + 1 + exact_match_start
        end_res_subject = int(contact_res_subject) - int(row['hit_i'].split(',')[0]) + 1 + exact_match_end + 1

        motif_file_subject, seq_subject = extract_motif_pdb(pid_subject, subject_ch, start_res_subject, end_res_subject, row['subject'], config)

        fasta_file = config['dirs']['motifs_pdb'] + '{}_{}.fasta'.format(seq_target, seq_subject)
        save_fasta_file(seq_target, seq_subject, fasta_file)

        process = subprocess.Popen(
            ['TMalign', motif_file_target, motif_file_subject,
             '-I', fasta_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        rmsd, TM_score = parse_TM_align(stdout)
        all_rmsds.append(rmsd)
        all_TM_scores.append(TM_score)

    df['RMSD_motifs'] = all_rmsds
    df['TM_score'] = all_TM_scores

    df.to_csv(config['out_files']['TM_align_all'] , sep='\t', index=False)
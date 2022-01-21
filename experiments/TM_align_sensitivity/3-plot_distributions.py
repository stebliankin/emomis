from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import zscore
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

def closest(lst, K):
    # find the closest number in a list
    return min(range(len(lst)), key=lambda i: abs(lst[i] - K))

sns.set_style("whitegrid")

rmsd_table = 'out_rmsd/rmsd_distribution.csv'
rmsd_thresholds = 'out_rmsd/thresholds.csv'
zscore_threshold_high = 1.645 # p-value is 5% - high confidence threshold
zscore_threshold_medium =  1.282 # p-value is 10%
#zscore_threshold_medium =  1.037 # p-value is 15% - medium confidence threshold


rmsd_df = pd.read_csv(rmsd_table)

# Plot RMSD distribution
fig, ax = plt.subplots()
fig.set_size_inches(10, 5)

sns.boxplot(y='RMSD', x='Len', data=rmsd_df)
fig.savefig('figures/rmsd_score_violin.jpg', format='jpg', dpi=1200)

fig.clear()

# for i in range(3, 11):`
#     tmp_df = rmsd_df[rmsd_df['Len']==i]
#     percentile = np.percentile(tmp_df['RMSD'], 5)
#     print("RMSD: 5% percentile for len {} is {}".format(i, percentile))

def get_rmsd_threshold(zscore_threshold):
    rmsd_thresholds = []
    for i in range(3, 33):
        tmp_df = rmsd_df[rmsd_df['Len'] == i]
        tm_scores = list(tmp_df['RMSD'])
        tm_scores = sorted(tm_scores)
        tmp_zscores = zscore(tm_scores)
        threshold_i = closest(tmp_zscores, -zscore_threshold)
        rmsd_thresholds.append(tm_scores[threshold_i+1])
        print("len {}: {}".format(i, tm_scores[threshold_i+1]))
    return rmsd_thresholds

rmsd_thresholds_high = get_rmsd_threshold(zscore_threshold_high)
rmsd_thresholds_medium = get_rmsd_threshold(zscore_threshold_medium)
with open(rmsd_thresholds, 'w') as out:
    out.write('Len,RMSD_high_conf,RMSD_medium_conf\n')
    for len_i in range(3, 33):
        out.write('{},{},{}\n'.format(len_i, rmsd_thresholds_high[len_i-3], rmsd_thresholds_medium[len_i-3]))

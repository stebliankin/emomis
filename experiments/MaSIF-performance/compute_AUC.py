
import numpy as np
import os
from sklearn import metrics
import pdb
import sys
import importlib
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from default_config.masif_opts import masif_opts

params = masif_opts["ppi_search"]
custom_params_file_list = sys.argv[1].split(',')
model_distances_df = pd.DataFrame({"Distance":[], "label":[], "model":[]})
out_neg_dist_dir = 'distribution/'

for custom_params_file in custom_params_file_list:
    custom_params = importlib.import_module(custom_params_file, package=None)
    custom_params = custom_params.custom_params

    masif_train = 'lists/masif_original_train.txt' # to exclude complexes used by MaSIF training

    masif_train_list = [x.strip('\n').split('_')[0] for x in open(masif_train).readlines()]

    for key in custom_params:
        print('Setting {} to {} '.format(key, custom_params[key]))
        params[key] = custom_params[key]

    def compute_roc_auc(pos, neg):
        labels = np.concatenate([np.ones((len(pos))), np.zeros((len(neg)))])
        dist_pairs = np.concatenate([pos, neg])
        scores = np.array([1/(x+0.001) for x in dist_pairs])

        auc_score = metrics.roc_auc_score(labels, scores)

        fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
        thresholds = 1/thresholds
        precisions = tpr / (tpr + fpr + 0.00001)
        #print(precisions)
        #print(thresholds)

        #optimal_threshold = thresholds[np.argmax(precisions)]
        #print("Optimal threshold: {}".format(optimal_threshold))
        return auc_score

    def compute_precision_recall(pos, neg, cutoff):
        labels = np.concatenate([np.ones((len(pos))), np.zeros((len(neg)))])
        dist_pairs = np.concatenate([pos, neg])
        predicted=[]
        for el in dist_pairs:
            if el<cutoff:
                predicted.append(1)
            else:
                predicted.append(0)
        return metrics.precision_score(labels, predicted),\
               metrics.recall_score(labels, predicted)

    def search_threshold_recall(pos, neg, target_recall):
        threshold_curr = np.max(pos)
        recall_curr = 1
        while recall_curr>target_recall:
            precision_curr, recall_curr = compute_precision_recall(pos, neg, threshold_curr)
            threshold_curr-=0.001
        #pdb.set_trace()
        return threshold_curr


    def compute_distances(p1_desc, p2pos_desc, p2neg_desc):
        pos_distances = []
        neg_distances = []
        for i in range(0, p1_desc.shape[0]):
            pos_distances.append(np.sqrt(np.sum(np.square(p1_desc[i] - p2pos_desc[i]))))
            neg_distances.append(np.sqrt(np.sum(np.square(p1_desc[i] - p2neg_desc[i]))))
        return pos_distances, neg_distances

    if not os.path.exists('metrics/'):
        os.mkdir('metrics/')

    model_name = custom_params_file.replace('.custom_params','').replace('nn_models.', '')
    out_metrics = 'metrics/'+ model_name + '.tsv'

    descriptors_dir = params['desc_dir']
    cache_dir = params['cache_dir']
    try:
        test_indx = np.load(cache_dir+'/pos_test_idx.npy').astype(int)
    except:
        pdb.set_trace()
    val_indx = np.load(cache_dir+'/pos_val_idx.npy').astype(int)

    p1_desc = np.load(descriptors_dir+'/p1_desc_flipped.npy')
    p2pos_desc = np.load(descriptors_dir+'/p2_desc_straight.npy')
    p2neg_desc = np.load(descriptors_dir+'/p2neg_desc_straight.npy')
    pos_names = np.load(cache_dir + '/p1target_names.npy')

    unique_complexes = set()
    # Remove MaSIF training complexes from the test set:
    test_indx_filtered = []
    for i in test_indx:
        patch_i = pos_names[i]
        pid = patch_i.split('_')[0]
        if pid not in masif_train_list:
            test_indx_filtered.append(i)

            unique_complexes.add(pid)

    print("Number of test patches before filtering: {}".format(len(test_indx)))
    print("Number of test patches after excluding masif train set: {}".format(len(test_indx_filtered)))

    print("Total number of complexes: {}".format(len(unique_complexes)))

    p1_desc_val = p1_desc[test_indx_filtered]
    p2pos_desc_val = p2pos_desc[test_indx_filtered]
    p2neg_desc_val = p2neg_desc[test_indx_filtered]

    pos_distances_val, neg_distances_val = compute_distances(p1_desc_val, p2pos_desc_val, p2neg_desc_val)
    mean_neg_val = np.mean(neg_distances_val)
    mean_pos_val = np.mean(pos_distances_val)

    with open(out_neg_dist_dir+ 'MaSIF_neg_scores_' + custom_params_file + '.csv', 'w') as out:
        #out.write('MaSIF_neg_score\n')
        for dist in neg_distances_val:
            out.write(str(dist)+'\n')

    with open(out_neg_dist_dir+ 'MaSIF_pos_scores_' + custom_params_file + '.csv', 'w') as out:
        for dist in pos_distances_val:
            out.write(str(dist)+'\n')
    #optimal_threshold = (mean_neg_val + mean_pos_val)/2
    #optimal_threshold=2.964738130569458

    #print("Mean pos val: {}; Mean neg val: {}; Optimal cutoff: {}".format(mean_pos_val, mean_neg_val, optimal_threshold))

    p1_desc_test = p1_desc[test_indx_filtered]
    p2pos_desc_test = p2pos_desc[test_indx_filtered]
    p2neg_desc_test = p2neg_desc[test_indx_filtered]
    pos_distances_test, neg_distances_test = compute_distances(p1_desc_test, p2pos_desc_test, p2neg_desc_test)

    #optimal_threshold = search_threshold_recall(pos_distances_test, neg_distances_test, 0.95)
    optimal_threshold = 3
    print("Max descriptors distance: {}".format(np.max(pos_distances_test)))
    print("Threshold when recall is 95%: {}".format(optimal_threshold))

    roc = compute_roc_auc(pos_distances_test, neg_distances_test)
    precision, recall = compute_precision_recall(pos_distances_test, neg_distances_test, optimal_threshold)
    print("ROC AUC: {}; Precision: {}; Recall: {};".format(roc, precision, recall))

    with open(out_metrics, 'w') as out:
        out.write('Model\tROC\tPrecision\tRecall\n')
        out.write(('{}\t{}\t{}\t{}\n'.format(model_name, roc, precision, recall)))

    model_name = custom_params_file.split('.')[1]
    pos_df = pd.DataFrame({'Distance':pos_distances_test})
    pos_df['label'] = 'Positive'
    pos_df['model'] = model_name
    neg_df = pd.DataFrame({'Distance': neg_distances_test})
    neg_df['label'] = 'Negative'
    neg_df['model'] = model_name

    model_distances_df = pd.concat([model_distances_df, neg_df, pos_df])


if not os.path.exists('figures'):
    os.mkdir('figures')

sns.set_theme(style="whitegrid", palette="pastel")

#
ax1 = sns.violinplot(data=model_distances_df, y='Distance', x='model', hue='label', palette=sns.color_palette('bright')[:2], scale="width")
ax1.set_ylim(0,8)
#ax1.axhline(optimal_threshold, ls='--')
#ax1.axhline(2.8, ls='--')
#ax1.axhline(0.98, ls='--',color='r')
#ax1.axhline(1.77, ls='--',color='r')

plt.savefig('figures/distance_violin.png', dpi=300)

#ax1 = sns.violinplot(data=pos_df, y='Distance', x='model', palette=sns.color_palette('bright')[1:2], scale="width")

#ax1.set_ylim(0,8)
#ax1.axhline(optimal_threshold, ls='--')
# ax1.axhline(7.51, ls='--', color='black')
# ax1.axhline(2.78, ls='--',color='g')
# ax1.axhline(7.25, ls='--',color='black')
# ax1.axhline(4.68, ls='--',color='g')

#plt.savefig('figures/distance_violin_pos.png', dpi=300)


model_name = "standard-sc05-min3-nomax"
custom_params = {}
custom_params['cache_dir'] = 'nn_models/{}/cache_sars2/'.format(model_name)
custom_params['model_dir'] = 'nn_models/{}/all_feat/model_data/'.format(model_name)
custom_params['desc_dir'] = 'descriptors/{}/all_feat/sars2/'.format(model_name)
custom_params['cached_list'] = "lists/cached_{}_sars2.txt".format(model_name)

custom_params['training_list'] = 'lists/train_small.txt'
custom_params['testing_list'] = 'lists/sars2_nonredundant.txt'

custom_params['min_per_complex'] = 3
custom_params['max_per_complex'] = float('inf')

custom_params['min_sc_filt'] = 0.5
custom_params['max_sc_filt'] = 1.0
custom_params['gpu_id'] = "/gpu:2"
#custom_params["result_second_stage"] = "result_second_stage/results_masif_sc02-datlsa_{}-model.txt".format(model_name)

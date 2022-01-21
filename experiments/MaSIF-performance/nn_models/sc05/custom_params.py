custom_params = {}
model_name = 'sc05'
cache_model = 'standard-sc05-min3-nomax'
print('Here')
custom_params['cache_dir'] = 'nn_models/{}/cache/'.format(cache_model)
custom_params['model_dir'] = 'nn_models/sc05/all_feat/model_data/'
custom_params['desc_dir'] = 'descriptors/{}/all_feat/all_sabdab/'.format(model_name)
custom_params['gif_eval_out'] = 'nn_models/sc05/gif_eval/'
custom_params['min_sc_filt'] = 0.5
custom_params['max_sc_filt'] = 1.0
custom_params['pos_surf_accept_probability'] = 1.0



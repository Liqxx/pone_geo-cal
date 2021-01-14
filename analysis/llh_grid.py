# usual
import os
import sys
import copy
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, interp2d

# llh stuff
import llh_main

### user parameters
# inputs
direc          = sys.argv[1]
paramfile      = os.path.join(direc, 'out_params.npy')

out_truth      = os.path.join(direc, 'out_truth.npy')
out_data       = os.path.join(direc, 'out_data.npy')
out_data_nn    = os.path.join(direc, 'out_data-nn.npy')

df_truth       = os.path.join(direc, 'out_df-truth.pkl')
df_data        = os.path.join(direc, 'out_df-data.pkl')
df_data_nn     = os.path.join(direc, 'out_df-data-nn.pkl')

# outputs
out_indiv      = os.path.join(direc, 'out_llh-array-indiv')
out_sum        = os.path.join(direc, 'out_llh-array-sum')

# parameters
n_data         = 20 # number of bins to use if smart_binning is off
use_df         = False # whether to use dataframes or dict
store_indiv    = True # whether to store a dict for per-DOM per-POCAM LLH
use_nn_data    = True # whether to use neural net prediction as data or simulation 

# get data
print('Loading data...')

# stdout information
if use_df:
    print('\t- using dataframes.')
else:
    print('\t- using binned numpy arrays.')
if use_nn_data:
    print('\t- using neural net generated data.')
    out_indiv += '-nn-data'
    out_sum   += '-nn-data'    

# now load data
# raw
df_t  = pd.read_pickle(df_truth) if use_df else []
df_d  = pd.read_pickle(df_data) if use_df else []
out_t = np.load(out_truth) if not use_df else []
out_d = np.load(out_data) if not use_df else []
# neural net
df_d  = df_d if not use_nn_data else pd.read_pickle(df_data_nn) if use_df else []
out_d = out_d if not use_nn_data else np.load(out_data_nn, allow_pickle=True)

print('... done!')

### read in
# get simulation set parameters
params = np.load(paramfile, allow_pickle=True, encoding='latin1').item()

# get keys
pocams       = params['pocams']
all_dom_keys = params['oms']

# get binnning
binsize       = params['binsize']
smart_binning = params['smart_binning']
n_bins_new    = params['n_bins_new']

if smart_binning:
    n_data = n_bins_new

# get 5-dimensional llh array summed over DOMs
if not store_indiv:
    llh_arr5d = llh_main.llh_grid5d(df_t, df_d, 
                                    out_t, out_d,
                                    params,
                                    pocams, all_dom_keys,
                                    df=use_df, 
                                    store_indiv=False,
                                    smart_binning=smart_binning, n_data=n_data)

    # save
    np.save(out_sum, llh_arr5d)

if store_indiv:
    llh_arr5d, llh_arr5d_iv = llh_main.llh_grid5d(df_t, df_d, 
                                                  out_t, out_d,
                                                  params,
                                                  pocams, all_dom_keys,
                                                  df=use_df, 
                                                  store_indiv=True,
                                                  smart_binning=smart_binning, n_data=n_data)

    # save
    np.save(out_sum, llh_arr5d)
    np.save(out_indiv, llh_arr5d_iv)



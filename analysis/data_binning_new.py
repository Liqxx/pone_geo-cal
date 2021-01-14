# icecube 
from icecube.dataclasses import *
from icecube.icetray import OMKey
from icecube import dataio, dataclasses, simclasses

# usual
import os
import sys
import copy
import joblib
import numpy as np
from scipy.interpolate import interp1d, interp2d

# storing dataframe
import pandas as pd

# parameter string for files
def param_string(abso, sca, domeff, p0, p1, Nph, N, pocam):
    string = 'ABS-%.5f_SCA-%.5f_DOME-%.5f_P0-%.5f_P1-%.5f_NPH-%.3e_N-%i_POCAM-%s' %(abso, sca, domeff, p0, p1, Nph, N, pocam)
    return string

########################################################################
### user parameters
########################################################################
direc         = sys.argv[1]
file_tag      = direc.rstrip('/').split('/')[-1].split('_')[1]
out_dir       = '/data/user/fhenningsen/deepcore_data/binned_%s' %file_tag
binsize       = 20 # new bin size in ns if not smart_binning
smart_binning = True # enables exponential binning of time profiles
store_df      = True # stores output as dataframe

print('Reading files from: %s' %direc)

# create output folder
if not os.path.exists(out_dir):
    print('Creating output directory: %s\n' %(out_dir))
    os.makedirs(out_dir)

########################################################################
### read in of meta data
########################################################################
# get simulation set parameters
params = np.load(os.path.join(direc, 'PARAMS.npy')).item()

# get global dicts
pocams       = params['pocams']
all_dom_keys = [OMKey(int(i.split('-')[0]), int(i.split('-')[1]), 0) for i in params['oms']]
truth_arr    = params['truth_arr']
truth_str    = params['truth_string']

# get geometry
geofile  = params['geo_file']
geometry = dataio.I3File(geofile)
gframe   = geometry.pop_frame()  
geo      = gframe["I3Geometry"] #access geo file via key

# generate event dict
event = {}
for i in all_dom_keys:
    om = '%i-%i' %(i.string, i.om)
    event[om] = []

# get paramater values
grid = params['grid']

# data / scan
scan_dict  = params['scan_dict']
scan_N     = scan_dict['n']
scan_Nph   = scan_dict['nph']

# truth
truth_dict = params['truth_dict']
t_N     = truth_dict['n']
t_Nph   = truth_dict['nph']
t_abs   = truth_dict['abs']
t_sca   = truth_dict['sca']
t_dome  = truth_dict['domeff']
t_p0    = truth_dict['p0']
t_p1    = truth_dict['p1']

# number of combinations
combinations = len(grid)

# get number of read-in bins
n_bins_old    = params['time_bins'] # in 1ns bins

########################################################################
### use meta data to adjust binning and save params
########################################################################
# smart binning config
if not smart_binning:
    n_bins_new = 5000/binsize
    time_bins  = np.linspace(0, n_bins_old, (n_bins_new+1))
    print 'Using binsize of', binsize, 'ns\n'
else:
    n_bins_new = 30
    print 'Using smart binning.\n'

# record binning in param dict    
params['binsize']       = binsize
params['smart_binning'] = smart_binning
params['n_bins_new']    = n_bins_new

########################################################################
### read-in truth
########################################################################
# number of rows for truth
nrows   = len(pocams) * len(all_dom_keys)

# arrays for meta data
dists    = np.zeros((nrows))
zeniths  = np.zeros((nrows))
t_geos   = np.zeros((nrows))
all_bins = np.zeros((nrows, n_bins_new + 1))

# arrays for data
x     = np.zeros((nrows, n_bins_new)) # data bins
theta = np.zeros((nrows, 12)) # use [pocam_index, pocam_string, pocam_om, string, dom, distance, zenith, abs, sca, dome, p0, p1]

cntr = 0
for pc, pocam in enumerate(pocams):
    par_t = param_string(t_abs, t_sca, t_dome, t_p0, t_p1, t_Nph, t_N, pocam)
    data_tmp_1ns_binned = np.load(os.path.expanduser('%s/TRUTH_%s.npy'%(direc, par_t)))
    
    # get pocam meta
    pc_s = int(pocam.split('-')[0])
    pc_o = int(pocam.split('-')[1])
    
    for key in all_dom_keys:
        omkey = '%i-%i' %(key.string, key.om)
        d  = data_tmp_1ns_binned.item().get(key)
        
        # create dictionary structure
        pocam_om = '%s_%s' %(pocam, omkey)
        
        if not smart_binning:
            # re-sum bins to get new bins of size <binsize>
            data_tmp_binned = [sum(d[i*binsize : (i+1)*binsize]) for i in range(n_bins_new)]
        else:
            # get POCAM coordinates
            p   = pocam.split('-')
            pk  = OMKey(int(p[0]), int(p[1]), 0)
            p_x = geo.omgeo[pk].position.x
            p_y = geo.omgeo[pk].position.y
            p_z = geo.omgeo[pk].position.z

            # get DOM coordinates
            d_x = geo.omgeo[key].position.x
            d_y = geo.omgeo[key].position.y
            d_z = geo.omgeo[key].position.z

            # get distance POCAM/DOM
            distance = np.sqrt((d_x-p_x)**2 + (d_y-p_y)**2 + (d_z-p_z)**2)
            dists[cntr] = distance
            
            # get zenith
            dz     = d_z - p_z
            zenith = np.arccos(abs(dz)/distance) if dz > 0 else np.pi-np.arccos(abs(dz)/distance)
            zeniths[cntr] = zenith
            
            # get geometrical min time of photons
            c_ice = 3e8 / 1.32
            t_geo = int(np.floor(distance / c_ice * 1e9)) # [ns]
            idx   = next((i for i, x in enumerate(d) if x), None)
            t_geos[cntr] = t_geo
            
            # get new smart bins
            smart_exp  = np.log(n_bins_old - t_geo) / np.log(n_bins_new)
            smart_bins = [i**smart_exp for i in range(n_bins_new + 1)]
            smart_bins = np.array([t_geo + int(round(i)) for i in smart_bins])
            all_bins[cntr] = smart_bins
            
            # re-sum data with smart bins
            data_tmp_binned = [sum(d[smart_bins[i]:smart_bins[i+1]]) for i in range(len(smart_bins) - 1)]

        # eventually put to array
        x[cntr]     = data_tmp_binned
        theta[cntr] = [pc+1, pc_s, pc_o, key.string, key.om, distance, zenith, t_abs, t_sca, t_dome, t_p0, t_p1]
            
        # increment nn input counter
        cntr += 1
        
# save
data = np.hstack([x, theta])
np.save(os.path.join(out_dir, 'out_truth.npy'), data)

# update and save params
params['distance']   = dists 
params['zenith']     = zeniths
params['t_geo']      = t_geos
params['smart_bins'] = all_bins
np.save(os.path.join(out_dir, 'out_params.npy'), params)

# in case store df
if store_df:
    print('Saving NN inputs...')
    df           = pd.DataFrame(data)
    meta_columns = ['pocam_index', 'pocam_string', 'pocam_om', 'string', 'om', 'dist', 'zenith','abs', 'sca', 'domeff', 'p0', 'p1']
    df.columns   = ['x%s'%i for i in range(n_bins_new)] + meta_columns
    df.to_pickle(os.path.join(out_dir, 'out_df-truth.pkl')) 
    print('... done!\n')

print "Truth finished. Starting data..."

########################################################################
### read-in data
########################################################################
    
# number of rows for data
nrows = combinations * len(pocams) * len(all_dom_keys)

# arrays for data
x     = np.zeros((nrows, n_bins_new)) # data bins
theta = np.zeros((nrows, 12)) # use [pocam_index, pocam_string, pocam_om, string, dom, distance, zenith, abs, sca, dome, p0, p1]

cntr       = 0
combo_cntr = 1
for gp in grid:
    
    # disentangle grid point
    abso, sca, dome, p0, p1 = gp
        
    # print progress        
    print("%i / %i (parameter combo: %.5f %.5f %.5f %.5f %.5f)" %(combo_cntr, combinations, abso, sca, dome, p0, p1))
            
    pco_cntr = 0        
    for pc, pocam in enumerate(pocams): 
        par = param_string(abso, sca, dome, p0, p1, scan_Nph, scan_N, pocam)
        data_tmp_1ns_binned = np.load(os.path.expanduser('%s/%s.npy'%(direc, par)))
        
        # print progress
        print("\tPOCAM %s" %pocam)
        
        # get pocam meta
        pc_s = int(pocam.split('-')[0])
        pc_o = int(pocam.split('-')[1])

        for key in all_dom_keys:
            omkey = '%i-%i' %(key.string, key.om)
            d = data_tmp_1ns_binned.item().get(key)
            
            if not smart_binning:
                # re-sum bins to new binsize
                data_tmp_binned = [sum(d[i*binsize : (i+1)*binsize]) for i in range(n_bins_new)]
            else:
                smart_bins = [int(i) for i in all_bins[pco_cntr]]
                distance   = dists[pco_cntr]
                zenith     = zeniths[pco_cntr]
                t_geo      = t_geos[pco_cntr]
                # re-sum with new smart bins
                data_tmp_binned = [sum(d[smart_bins[i]:smart_bins[i+1]]) for i in range(len(smart_bins) - 1)]

            # store data
            x[cntr]     = data_tmp_binned
            theta[cntr] = [pc+1, pc_s, pc_o, key.string, key.om, distance, zenith, abso, sca, dome, p0, p1]
            
            # increment pocam-om counter
            pco_cntr += 1
            
            # increment global counter
            cntr += 1
        
    # increment combinations counter
    combo_cntr += 1
                    
print "... done!\n"

# save
data = np.hstack([x, theta])
np.save(os.path.join(out_dir, 'out_data.npy'), data)

# in case store nn input vectors
if store_df:
    print('Saving NN inputs...')
    df          = pd.DataFrame(data)
    meta_columns = ['pocam_index', 'pocam_string', 'pocam_om', 'string', 'om', 'dist', 'zenith','abs', 'sca', 'domeff', 'p0', 'p1']
    df.columns   = ['x%s'%i for i in range(n_bins_new)] + meta_columns
    joblib.dump(df, os.path.join(out_dir, 'out_df-data.pkl'), compress=True) 
    print('... done!\n')

#This script reads in data from i3 files and bins it into 1ns to reduce file size and save memory for later read ins.

from icecube import dataio, dataclasses, simclasses
from icecube.icetray import OMKey
from icecube.dataclasses import *

import os
import gc
import sys
import copy
import joblib
import numpy as np
import pandas as pd
from os.path import join, exists, expanduser

### Define paths ###
data_dir  = sys.argv[1].rstrip('/')
tag       = data_dir.split('/')[-1].split('_')[1]
out_dir   = '/data/user/fhenningsen/pone_data/read_%s' %tag

# create output folder
if not exists(out_dir):
    print('Creating output directory in {}'.format(out_dir))
    os.makedirs(out_dir)

### Get parameters ###
read_in_data  = True

n          = 100  # data runs to be read in
bins       = 5000 # time bins for 5000ns
bin_length = bins
range_bins = (0, 5000)

alldata  = np.array([np.array(i[:-7].split('_')) for i in os.listdir(data_dir) if ".i3.bz2" in i and not "job" in i])

# access parameters via slicing
NPH  = float(alldata[0][7][4:])
N    = float(alldata[0][8][2:])
POC  = np.unique([str(alldata[i][9][6:]) for i in range(len(alldata))]) 

# gcd file 
gcd = '_'.join(alldata[0][1:4])

# get grid
grid_file = join(data_dir, 'GRID.npy')
if not exists(grid_file):
    print('Error: No parameter grid %s found.' %(grid_file))
    sys.exit()
grid = np.load(grid_file)

if n > N:
    print('Error: Cannot read in more runs than generated.')
    sys.exit()

### Define empty event dicionary for time stamps ###
geofile  = "/home/fhenningsen/gcd/%s" %(gcd)
geometry = dataio.I3File(geofile)
gframe = geometry.pop_frame()
geo = gframe["I3Geometry"]

all_dom_keys = []
for i in geo.omgeo.keys():
    if i.pmt == 0 and i.string < 87 and i.om <= 60: # ignore upgrade modules and icetop
        all_dom_keys.append(i)

event = {}
for i in all_dom_keys:
    event[i] = []

### Read in data ###

# number of parameters
combinations = len(grid)

# number of rows for data
nrows = combinations * len(POC) * len(all_dom_keys)

# arrays for data
x        = np.zeros((nrows, bins)) # data bins
theta    = np.zeros((nrows, 10)) # use [pocam_index, pocam_string, pocam_om, string, dom, distance, zenith, abs, sca, dome]
dists    = np.zeros((nrows))
zeniths  = np.zeros((nrows))

if read_in_data:

    combinations = len(grid)

    cntr = 0
    combo_cntr = 1
    for gp in grid:
        
        abso, sca, dome = gp
        
        # print progress
        print("%i / %i (parameters: %.3f, %.3f, %.3f)" %(combo_cntr, combinations, abso, sca, dome))
        
        for pc, pocam in enumerate(POC):
            
            # get pocam meta
            pc_s = int(pocam.split('-')[0])
            pc_o = int(pocam.split('-')[1])
            
            # define parameter string
            param_str = 'GCD_%s_ABS-%.3f_SCA-%.3f_DOME-%.3f_NPH-%.3e' %(gcd, abso, sca, dome, NPH)
            data_file = expanduser('%s/%s_N-%i_POCAM-%s.i3.bz2'%(data_dir, param_str, N, pocam))
            out_file  = expanduser('%s/%s_N-%i_POCAM-%s.npy'%(out_dir, param_str, N, pocam))
                        
            if exists(data_file) and not exists(out_file):
                
                # print progress
                print("\tPOCAM %s" %pocam)
                
                data_tmp = copy.deepcopy(event)
                
                # first we read in the i3 file
                data_tmp_i3 = dataio.I3File(data_file)
                for i in range(n): # iterate over frames
                    data_tmp_fr = data_tmp_i3.pop_daq()
                    data_tmp_k = data_tmp_fr['MCPESeriesMap']
                    for j in data_tmp_k.iteritems():
                        dom_key=j[0]
                        if dom_key.string < 87: # ignore upgrade modules
                            data_tmp[dom_key] += [j.time for j in j[1]] # append time stamps

                # .. now bin into 1ns bins
                data_binned_tmp = copy.deepcopy(event)
                for key in all_dom_keys:  # go through all DOMs
                    
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
                    
                    if len(data_tmp[key]) > 0: # only look at hit DOMs
                        # doing the actual binning:
                        hist_tmp, edges  = np.histogram(data_tmp[key], bins=bins, range=range_bins)
                    else:
                        hist_tmp = np.zeros(bin_length)
                    
                    # store data
                    x[cntr]     = hist_tmp
                    theta[cntr] = [pc+1, pc_s, pc_o, key.string, key.om, distance, zenith, abso, sca, dome]
                    
                    # clean up RAM
                    gc.collect()

                    # advance
                    cntr += 1
                    
        # increment combinations counter
        combo_cntr += 1

# save numpy
print('Storing npy file...')
data = np.hstack([x, theta])
np.save(os.path.join(out_dir, 'out_data.npy'), data)

# save df
print('Storing dataframe file...')
df           = pd.DataFrame(data)
meta_columns = ['pocam_index', 'pocam_string', 'pocam_om', 'string', 'om', 'dist', 'zenith','abs', 'sca', 'domeff']
df.columns   = ['x%s'%i for i in range(bins)] + meta_columns
joblib.dump(df, os.path.join(out_dir, 'out_df-data.pkl'), compress=True) 
        
### Store parameters ###
params = {'pocams'       : POC,
          'oms'          : ['%i-%i' %(i.string, i.om) for i in all_dom_keys],
          'geo_file'     : geofile,
          #'pocam_keys'   : [OMKey(int(p.split('-')[0]), int(p.split('-')[1]), 0) for p in pocams],
          #'om_keys'      : [i for i in all_dom_keys],
          'scan_dict'    : {'n':N,   'nph':NPH},
          'time_bins'    : bins,
          'grid'         : grid,
          'n_data_read'  : n,
          'distance'     : dists,
          'zenith'       : zeniths,
          }
np.save(expanduser('%s/PARAMS.npy' %(out_dir)), params)

print('\nFinished!')







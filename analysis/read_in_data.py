#This script reads in data from i3 files and bins it into 1ns to reduce file size and save memory for later read ins.

from icecube import dataio, dataclasses, simclasses
from icecube.icetray import OMKey
from icecube.dataclasses import *

import os
import gc
import sys
import copy
import numpy as np
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
range_bins = (0, bins)

alldata  = np.array([np.array(i[:-7].split('_')) for i in os.listdir(data_dir) if ".i3.bz2" in i and not "TRUTH" in i])

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
          }
np.save(expanduser('%s/PARAMS.npy' %(out_dir)), params)

### Read in data ###

if read_in_data:

    combinations = len(grid)

    cntr = 1
    for gp in grid:
        
        abso, sca, dome = gp
        
        # print progress
        print("%i / %i (parameters: %.3f, %.3f, %.3f)" %(cntr, combinations, abso, sca, dome))
        
        for pocam in POC:
            
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
                    
                    if len(data_tmp[key]) > 0: # only look at hit DOMs
                        # doing the actual binning:
                        hist_tmp, edges  = np.histogram(data_tmp[key], bins=bins, range=range_bins)
                    else:
                        hist_tmp = np.zeros(bin_length)
                    
                    data_binned_tmp[key] = hist_tmp

                # save the dictionary:
                np.save(out_file, data_binned_tmp)
                
                # clean up RAM
                gc.collect()

        # increment combinations counter
        cntr += 1

print('\nFinished!')







#!/usr/bin/env python

import os
import sys
import time
import numpy as np
import datetime
from os.path import join, exists

# who
datauser = '/data/user/fhenningsen'

# executable
executable = '/home/fhenningsen/pone_geo-cal/data_generation/sim_executable.sh'
mem        = 10 #GB

# tag for simulation
sim_tag = 'water_test'

# create time-sensitive folder for output
t_str = datetime.datetime.now().isoformat('_')[:-7].replace(':', '-')
out_folder = join(datauser, 'deepcore_data' , 'sim_' + sim_tag + '_' + t_str)
if not exists(out_folder):
    print('Creating output directory in {}'.format(out_folder))
    os.makedirs(out_folder)

# create time-sensitive folder for log files
out_log = join('/scratch/fhenningsen', 'logs', t_str)
if not exists(out_log):
    print('Creating log directory in {}'.format(out_log))
    os.makedirs(out_log)

# create time-sensitive folder for log files
out_sub = join('/scratch/fhenningsen', 'subs', t_str)
if not exists(out_sub):
    print('Creating log directory in {}'.format(out_sub))
    os.makedirs(out_sub)

# fixed inputs
N        = 100
N_truth  = N * 1
N_photon = 1e8

# submit 
submit_data  = True
submit_truth = True

# running config parameters
n_par = 3
ABS   = np.linspace(0.9, 1.1,    n_par)
SCA   = np.linspace(0.9, 1.1,    n_par)
DOME  = np.linspace(0.9, 1.1,    n_par)
P0    = np.linspace(-1.0, 1.0,   n_par)
P1    = np.linspace(-0.20, 0.30, n_par)

# fill parameter grid
grid_filled = np.zeros((len(ABS) * len(SCA) * len(DOME) * len(P0) * len(P1)))

cntr = 0
for abso in ABS:
    for sca in SCA:
        for dome in DOME:
            for p0 in P):
                for p1 in P1:
                        grid_filled[cntr] = [abso, sca, dome, p0, p1]
                        cntr += 1

grid_filled = np.round(grid_filled, 3)
np.save(join(out_folder, 'GRID.py'), grid_filled)

# truth
tabs  = 1.00
tsca  = 1.00
tdome = 1.01
tp0   = -0.03
tp1   = 0.05

print('Number of jobs: %i' %(1 + len(ABS) * len(SCA) * len(DOME) * len(P0) * len(P1)))

# simulation sets
for _abs in ABS:

    for sca in SCA:

        for dome in DOME:

            for p0 in P0:

                for p1 in P1:

                    out_file    = 'ABS-%.3f_SCA-%.3f_DOME-%.3f_P0-%.3f_P1-%.3f_NPH-%.3e_N-%i' %(_abs, sca, dome, p0, p1, N_photon, N)
                    args        = "%i %s %.3f %.3f %.3f %.3f %.3f %i" %(N, join(out_folder, out_file), _abs, sca, dome, p0, p1, N_photon)
                    env         = "HDF5_USE_FILE_LOCKING='FALSE'"
                    log_str     = 'job_deepcore_%s' %(out_file)
                    submit_info = 'executable  = {script} \n\
                        universe    = vanilla \n\
                        initialdir = /home/fhenningsen \n\
                        request_gpus = 1 \n\
                        request_memory = {mem}GB \n\
                        log         = {outl}/{logs}.log \n\
                        output      = {out}/{logs}.out \n\
                        error       = {out}/{logs}.err \n\
                        arguments   = "{args}" \n\
                        environment = "{env}" \n\
                        transfer_executable = True \n\
                        queue 1 \n'.format(script = executable,
                                    mem    = mem,
                                    logs   = log_str,
                                    outl   = out_log,
                                    out    = out_folder,
                                    args   = args,
                                    env    = env,
                                    )

                    if submit_data:

                        # write submit file
                        sub_file = '%s/%s.sub' %(out_sub, log_str)
                        with open(sub_file, 'w') as f:
                            f.write(submit_info)

                        # submit them
                        os.system("condor_submit {}".format(sub_file))

# run truth
out_file    = 'TRUTH_ABS-%.3f_SCA-%.3f_DOME-%.3f_P0-%.3f_P1-%.3f_NPH-%.3e_N-%i' %(tabs, tsca, tdome, tp0, tp1, N_photon, N_truth)
args        = "%i %s %.3f %.3f %.3f %.3f %.3f %i" %(N_truth, join(out_folder, out_file),tabs, tsca, tdome, tp0, tp1, N_photon)
env         = "HDF5_USE_FILE_LOCKING='FALSE'"
log_str     = '%s' %(out_file)
submit_info = 'executable  = {script} \n\
    universe    = vanilla \n\
    initialdir = /home/fhenningsen \n\
    request_gpus = 1 \n\
    request_memory = {mem}GB \n\
    log         = {outl}/{logs}.log \n\
    output      = {out}/{logs}.out \n\
    error       = {out}/{logs}.err \n\
    arguments   = "{args}" \n\
    environment = "{env}" \n\
    transfer_executable = True \n\
    queue 1 \n'.format(script = executable,
                mem    = mem,
                logs   = log_str,
                outl   = out_log,
                out    = out_folder,
                args   = args,
                env    = env,
                )

if submit_truth:

    # write submit file
    sub_file = '%s/%s.sub' %(out_sub, log_str)
    with open(sub_file, 'w') as f:
        f.write(submit_info)

    # submit them
    os.system("condor_submit {}".format(sub_file))

print('Number of jobs: %i' %(1 + len(ABS) * len(SCA) * len(DOME) * len(P0) * len(P1)))

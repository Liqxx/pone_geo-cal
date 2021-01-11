#!/usr/bin/env python

import os
import sys
import time
import datetime
from os.path import join, exists

import numpy as np
from scipy import stats

# quasi-random
import quasi_random_rng

# who
datauser = '/data/user/fhenningsen'

# executable
executable = '/home/fhenningsen/osc/deepcore_systematics/data_generation/sim_executable.sh'
mem        = 10 #GB

# tag for simulation
sim_tag = '5d-nn-n-0-2000'

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
    print('Creating sub directory in {}'.format(out_sub))
    os.makedirs(out_sub)

# fixed inputs
N        = 100   # runs
N_truth  = N * 1 # truth runs
N_photon = 1e8   # number of photons per run
dim      = 5     # parameter dimensions
n_params = 2000  # number of (additional) points in n-d grid
n_grid   = 0     # use grid from this index on, if <n_grid> simulated already

# submit 
submit_data  = True
submit_truth = True

# truth
tabs  = 1.02
tsca  = 0.95
tdome = 1.05
tp0   = -0.08
tp1   = 0.06

############## RNG ##############

# create quasi-random input for NN
# quasi-random sampling
z = quasi_random_rng.rnew(dim, n_params + n_grid)

# use normal prior for abs, sca, domeff
# centered at 1, sigma 0.3
# truncate on +- 2 sigma
loc    = 1
scale  = 0.3
clip_a = loc - 2 * scale
clip_b = loc + 2 * scale

# get truncated interval
a, b   = (clip_a - loc) / scale, (clip_b - loc) / scale

# get prior distribution from uniform
# by using inverse CDF
# explanation: CDF transforms distribution to [0,1]
# inverse transforms back
transform = stats.truncnorm(a, b, loc=loc, scale=scale).isf

# replace abs, sca, and domeff priors 
z[:,0] = transform(z[:,0])
z[:,1] = transform(z[:,1])
z[:,2] = transform(z[:,2])

# set p0 range to [-1, 1] uniform
a      = - 1.0
b      =   1.0
z[:,3] = z[:,3] * (b - a) + a

# set p1 range to [-0.2, 0.2] uniform
a      = - 0.2
b      =   0.2
z[:,4] = z[:,4] * (b - a) + a

# optionally start grid from index if other points were simulated before
z = z[n_grid:]

# and save
np.save(join(out_folder, 'GRID.npy'), z)

#################################

print('Number of jobs: %i' %(n_params + 1))

for rng in z:
    
    # get grid point
    abso, sca, dome, p0, p1, = rng

    # submit file
    out_file    = 'ABS-%.5f_SCA-%.5f_DOME-%.5f_P0-%.5f_P1-%.5f_NPH-%.3e_N-%i' %(abso, sca, dome, p0, p1, N_photon, N)
    args        = "%i %s %.8f %.8f %.8f %.8f %.8f %i" %(N, join(out_folder, out_file), abso, sca, dome, p0, p1, N_photon)
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
out_file    = 'TRUTH_ABS-%.5f_SCA-%.5f_DOME-%.5f_P0-%.5f_P1-%.5f_NPH-%.3e_N-%i' %(tabs, tsca, tdome, tp0, tp1, N_photon, N_truth)
args        = "%i %s %.8f %.8f %.8f %.8f %.8f %i" %(N_truth, join(out_folder, out_file),tabs, tsca, tdome, tp0, tp1, N_photon)
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

print('Number of jobs: %i' %(n_params + 1))

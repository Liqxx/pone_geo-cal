#!/usr/bin/env bash

#env:
eval $(/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh)

#configure:
metaproject='/home/fschmuckermaier/IceCube/meta-projects/combo/stable/build'
file_path="/data/user/fschmuckermaier/data_raw/ppc_test"
file_name="gcd_V5_test"
string=88
dom=72
n_runs=5
gcd="/home/fschmuckermaier/gcd/physics_volume_GCD.i3.bz2"
#gcd="/data/sim/IceCubeUpgrade/geometries/GCDs/GeoCalibDetectorStatus_ICUpgrade.v55.mixed_mergedGeo.V5.i3.bz2"

#run:
$metaproject/env-shell.sh python /home/fschmuckermaier/DeepCore_systematics/data_generation/simple_POCAM_simulation.py --gcd-file=$gcd  --output-i3-file="${file_path}/${file_name}" --number-of-runs=$n_runs --string=$string --dom=$dom


#Choose flashing POCAMs:  (all true POCAM positions were shifted by 1 downwards (e.g. true position = 87,4), since only OMs can be flashed somehow...)
# 87 4
# 87 84
# 88 2
# 88 72
# 89 2
# 89 10
# 89 13
# 89 38
# 89 107
# 90 12
# 90 14
# 90 100
# 91 15
# 91 50
# 92 6
# 92 18
# 92 28
# 93 8
# 93 17
# 93 64
# 93 113

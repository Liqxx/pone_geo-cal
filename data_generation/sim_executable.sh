#!/usr/bin/env bash

#env:
# eval $(/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh)
# metaproject='/home/fschmuckermaier/IceCube/meta-projects/combo/stable/build'
eval $(/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh)
metaproject='/home/fhenningsen/osc/py2-ext/build'

### 1. DEFINE INPUT PARAMETER ###

# number of runs
N=$1

# output file
file_name=$2

# Here, we define the 5 key parameters or rather define the lists of parameters we want to iterate through
water_shift_abs=$3
water_shift_sca=$4
dom_eff_corr=$5
#p0_ang_acc=$6
#p1_ang_acc=$7
phot=$6

# geometry file
gcd=$7

### 2. APPLYING THE PARAMETERS BY EDITING FILES ###
# Then we apply the parameters by editing the necessary files

# modification scripts
mod_scripts="/home/fhenningsen/pone/pone_geo-cal/modification_scripts"

# default water tables
water_default="/home/fhenningsen/pone/pone_geo-cal/watermodel/default"

# copy all the tables files so that a new submission doesnt change things
# and create a new tmp directory for the current submission
tmp_dir="tmp-a$3-b$4-de$5-p0$6-p1$7"
# make sure that truth is not deleted after normal jobs finished running and parameters were the same
if [[ $file_name == *"TRUTH"* ]]; then
    water_tables="/home/fhenningsen/pone/pone_geo-cal/watermodel/tmp/truth_${tmp_dir}"
else
    water_tables="/home/fhenningsen/pone/pone_geo-cal/watermodel/tmp/${tmp_dir}"
fi
# and copy all default files
mkdir -p $water_tables
cp ${water_default}/* ${water_tables}/.

# Multiplying the offset to the watermodel parameters:
$metaproject/env-shell.sh python ${mod_scripts}/change_water.py --a_corr=$water_shift_abs --b_corr=$water_shift_sca --outfile=${water_tables}/icemodel.dat

# Replace the overall dom efficiency factor in cfg.txt (Note: check if values greater than 1 work!):
sed "3s/.*/$dom_eff_corr/" ${water_default}/cfg.txt > ${water_tables}/cfg.txt

# taking p0 & p1 as input and produce an as.dat file in the water directory
# format: dima_from_unified.py <p0> <p1> <out_file>
# python ${mod_scripts}/dima_from_unified.py "$p0_ang_acc" "$p1_ang_acc" "$water_tables/as.dat"

### 3. RUN THE SIMULATION ###
#Then we run the simulation for all 7 relevant POCAMs with the new files

### Configure:  ###

### Running through all 7 POCAMs: ###
string=87
dom=84
$metaproject/env-shell.sh python /home/fhenningsen/osc/deepcore_systematics/data_generation/simple_POCAM_simulation.py --gcd-file=$gcd  --output-i3-file="${file_name}" --number-of-photons=$phot --number-of-runs=$N --string=$string --dom=$dom --ice-tables=$water_tables

string=88
dom=72
$metaproject/env-shell.sh python /home/fhenningsen/osc/deepcore_systematics/data_generation/simple_POCAM_simulation.py --gcd-file=$gcd  --output-i3-file="${file_name}" --number-of-photons=$phot --number-of-runs=$N --string=$string --dom=$dom --ice-tables=$water_tables

string=89
dom=38
$metaproject/env-shell.sh python /home/fhenningsen/osc/deepcore_systematics/data_generation/simple_POCAM_simulation.py --gcd-file=$gcd  --output-i3-file="${file_name}" --number-of-photons=$phot --number-of-runs=$N --string=$string --dom=$dom --ice-tables=$water_tables

string=90
dom=100
$metaproject/env-shell.sh python /home/fhenningsen/osc/deepcore_systematics/data_generation/simple_POCAM_simulation.py --gcd-file=$gcd  --output-i3-file="${file_name}" --number-of-photons=$phot --number-of-runs=$N --string=$string --dom=$dom --ice-tables=$water_tables

string=91
dom=50
$metaproject/env-shell.sh python /home/fhenningsen/osc/deepcore_systematics/data_generation/simple_POCAM_simulation.py --gcd-file=$gcd  --output-i3-file="${file_name}" --number-of-photons=$phot --number-of-runs=$N --string=$string --dom=$dom --ice-tables=$water_tables

string=92
dom=28
$metaproject/env-shell.sh python /home/fhenningsen/osc/deepcore_systematics/data_generation/simple_POCAM_simulation.py --gcd-file=$gcd  --output-i3-file="${file_name}" --number-of-photons=$phot --number-of-runs=$N --string=$string --dom=$dom --ice-tables=$water_tables

string=93
dom=64
$metaproject/env-shell.sh python /home/fhenningsen/osc/deepcore_systematics/data_generation/simple_POCAM_simulation.py --gcd-file=$gcd  --output-i3-file="${file_name}" --number-of-photons=$phot --number-of-runs=$N --string=$string --dom=$dom --ice-tables=$water_tables

# remove the current water dir
#rm -rf $water_tables

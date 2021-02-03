import numpy as np
import argparse
from optparse import OptionParser

parser = OptionParser(description="This script varies an icemodels abs and sca params by a factor.")
parser.add_option("--a_corr",default=1.0)
parser.add_option("--b_corr",default=1.0)
parser.add_option("--outfile")
(options, args) = parser.parse_args()

depth=[]
b=[]
a=[]
t=[]

f=open('/home/fhenningsen/workplace/pone/pone_geo-cal/watermodel/default/icemodel.dat','r')
lines=f.readlines()
for x in lines:
    depth.append(np.double(x.split(' ')[0]))

    ### Multiply factor to a and b ###
    b.append(float(options.b_corr)* np.double(x.split(' ')[1]))
    a.append(float(options.a_corr)* np.double(x.split(' ')[2]))
    t.append(np.double(x.split(' ')[3]))

f.close()

txt = open(options.outfile, 'w')
for i in range(len(b)):
    txt.write("%f %f %f %f\n" % (depth[i], b[i], a[i], t[i]))
txt.close()


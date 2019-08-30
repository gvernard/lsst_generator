import os
import sys
import numpy as np

# INPUT
# arg 1: the path to the directory containing the output files of the LSST light curve simulator, e.g. /path/to/the/directory/output/
# arg 2: the index of the trajectory from which continuous and sampled light curves have been extracted, starting from 0
# OUTPUT
# A directory named "lc_<index>" that contains the following files:
#    - 6 files with the magnification values of the continuous light curves having 2 columns: time (MJD) and magnitude
#    - 6 files with the magnification values of the sampled light curves having 3 columns: time (MJD), magnitude. and uncertainty


# Input
target_dir = sys.argv[1]
if not target_dir.endswith('/'):
    target_dir += '/'
index      = int(sys.argv[2])


# Initializing variables
outdir = str(index).zfill(6)
if not os.path.isdir(outdir):
    os.mkdir(outdir)
filters = ['u','g','r','i','z','y']


# Get parameters for transforming the x-values of the continuous light curves to time in MJD.
f = open(target_dir+'time_and_length.dat',"r")
lines = f.readlines()
f.close()
pix = float(lines[0].split(':')[1].split('[')[0].strip()) # pixel size [in 10^14 cm]
t0  = float(lines[2].split(':')[1].split('[')[0].strip()) # beginning of the obsercations [in MJD]
vels = np.loadtxt(target_dir+'velocities.dat')
vel  = vels[index,0] # effective velocity along the light curve trajectory in [km/s]
t_interval = 11574*pix/vel # the time interval of the continuous light curves [in days]


# Read the continuous light curve in all filters
Ntheo = np.loadtxt(target_dir+'theo_length.dat',dtype='int')
Npre = 0
for i in range(0,index):
    Npre += Ntheo[i]
Npost = Npre + Ntheo[index]
continuous_all = np.loadtxt(target_dir+'theo_mag.dat')
cont_t = t0 + np.arange(0,Ntheo[index])*t_interval # create the time values for the continuous light curves

for i in range(0,len(filters)):
    cont_m = continuous_all[Npre:Npost,i]
    np.savetxt(outdir+'/'+filters[i]+'_cont.dat',np.c_[cont_t,cont_m])


# Read the corresponding sampled light curve in the u filter
for i in range(0,len(filters)):
    sampled_m = np.loadtxt(target_dir+'filter_'+filters[i]+'_mag.dat')
    samp_m = sampled_m[index,:]

    sampled_dm = np.loadtxt(target_dir+'filter_'+filters[i]+'_dmag.dat')
    samp_dm = sampled_dm[index,:]

    sampled_t = np.loadtxt(target_dir+filters[i]+'_dates.dat')
    samp_t = sampled_t[:,0]

    np.savetxt(outdir+'/'+filters[i]+'_samp.dat',np.c_[samp_t,samp_m,samp_dm])

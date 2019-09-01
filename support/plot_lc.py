import sys
import numpy as np
import matplotlib.pyplot as plt

# INPUT
# arg 1: the path to an output directory of extract_lc.py

target_dir = sys.argv[1]
if not target_dir.endswith('/'):
    target_dir += '/'


filters = ['u','g','r','i','z','y']
mycolors = ['violet','blue','cyan','green','yellow','red']


fig, ax = plt.subplots(figsize=(9,5.5))

# Plotting the continuous light curves
for i in range(0,len(filters)):
    t,m = np.loadtxt(target_dir+filters[i]+'_cont.dat',unpack=True)
    ax.plot(t,-m,color=mycolors[i])


# Plotting the sampled light curves
for i in range(0,len(filters)):
    t,m,dm = np.loadtxt(target_dir+filters[i]+'_samp.dat',unpack=True)
    ax.errorbar(t,-m,yerr=0.02,fmt='none',color=mycolors[i])
    ax.scatter(t,-m,color=mycolors[i],label=filters[i],marker='*')

ax.legend(fontsize=17)
ax.set_xlabel('MJD',fontsize=20)
ax.set_ylabel('Magnitude',fontsize=20)
ax.tick_params(axis="both",labelsize=17)


plt.tight_layout()
plt.savefig('lc.pdf')


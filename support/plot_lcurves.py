import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os


indir  = sys.argv[1]

if not os.path.isdir(indir+'plots/'):
	os.makedirs(indir+'plots/')

m = sorted(glob.glob(indir+'tablet_*.dat'))
u = sorted(glob.glob(indir+'tableu_*.dat'))
g = sorted(glob.glob(indir+'tableg_*.dat'))
r = sorted(glob.glob(indir+'tabler_*.dat'))
i = sorted(glob.glob(indir+'tablei_*.dat'))
z = sorted(glob.glob(indir+'tablez_*.dat'))
y = sorted(glob.glob(indir+'tabley_*.dat'))

filters = [u,g,r,i,z,y]


for n,j in enumerate(m):
	fig,(ax1,ax2)=plt.subplots(1,2)
	a = np.loadtxt(j).T
	ax1.plot(a[0],a[1],'c',a[0],a[2],'b',a[0],a[3],'g',a[0],a[4],'y',a[0],a[5],'m',a[0],a[6],'r')
	ax1.invert_yaxis()
	ax1.set_ylabel(r'$\Delta$'+' [mag]')
	ax1.set_xlabel('pixel')
	ax1.legend([r'$\sigma_u$', r'$\sigma_g$',r'$\sigma_r$',r'$\sigma_i$',r'$\sigma_z$',r'$\sigma_y$'], loc='upper right')
        ymin,ymax = ax1.get_ylim()

	uu = np.loadtxt(u[n]).T
	gg = np.loadtxt(g[n]).T
	rr = np.loadtxt(r[n]).T
	ii = np.loadtxt(i[n]).T
	zz = np.loadtxt(z[n]).T
	yy = np.loadtxt(y[n]).T
	
	
	ax2.errorbar(uu[0],uu[1],yerr=uu[2],fmt='c')
	ax2.errorbar(gg[0],gg[1],yerr=gg[2],fmt='b')
	ax2.errorbar(rr[0],rr[1],yerr=rr[2],fmt='g')
	ax2.errorbar(ii[0],ii[1],yerr=ii[2],fmt='y')
	ax2.errorbar(zz[0],zz[1],yerr=zz[2],fmt='m')
	ax2.errorbar(yy[0],yy[1],yerr=yy[2],fmt='r')
	ax2.invert_yaxis()
        ax2.set_ylim(ymin,ymax)
	ax2.set_ylabel(r'$\Delta$'+' [mag]')
	ax2.set_xlabel('JD')
        for tick in ax2.get_xticklabels():
                tick.set_rotation(45)
	ax2.legend([r'$\sigma_u$', r'$\sigma_g$',r'$\sigma_r$',r'$\sigma_i$',r'$\sigma_z$',r'$\sigma_y$'], loc='upper right')
	
        plt.tight_layout()
	plt.savefig(indir+'plots/plot_'+str(n).zfill(4)+'.png')
	plt.close('all')

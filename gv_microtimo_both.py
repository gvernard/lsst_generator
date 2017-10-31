# -*- coding: utf-8 -*-
#custom version
import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt
#plt.switch_backend('agg')
from math import *
from scipy import ndimage, signal
from scipy.interpolate import interp1d
import time as pytime
from auxfunctions import veco, makeGaussian, m52snr, getKey, effe, readimg, Dis, irot, gv_readimg
import sys
import os,glob
import pyfits
from astropy.io import fits
import outfunctions
import warnings

#Function to plot magnification
#def micro(rc,numlc,years,s00,lam00,vp,inclination,zs,zl,kappa,gamma,s,sgal2,sdisp2,squas2,ra,dec,dra,ddec,tol,errbaseu,errbaseg,errbaser,errbasei,errbasez,errbasey,get_plot=False,write_table=False,write_img=True,mass=False,vcmb=False,nomacro=False):

def micro(outdir,source,magmap,lsst,other,output):
	starttime = pytime.clock()


	if lsst['years'] == 0:
		print "no reason to plot no time"
		sys.exit("no reason to plot zero time")


        #get magnification map
        img,RES,WIDTH = gv_readimg(magmap['path_to_maps'],str(magmap['id']))
        img=img/np.mean(img)
#	if (output['nomacro'] == True):
#		img=img/np.mean(img)
        print 'magnification map has been read'





	#Calculus of Einstein radius
	RE,DL,DS,DLS=Dis(other['zs'],other['zl'],magmap['mass'])




        if source['type'] == 'gaussian':

                # sigma zero
                s0=source['sigma0']*lt
                
                #real value filters
                leffu=3654.9/((10**10)*(1+other['zs']))
                leffg=4800.3/((10**10)*(1+other['zs']))
                leffr=6222.0/((10**10)*(1+other['zs']))
                leffi=7540.6/((10**10)*(1+other['zs']))
                leffz=8682.1/((10**10)*(1+other['zs']))
                leffy=9925.0/((10**10)*(1+other['zs']))
                lam0=source['lambda0']/((10**10))
                #getting sigma values on pixels
                su=s0*((leffu/lam0)**source['nu'])*RES/WIDTH/RE
                sg=s0*((leffg/lam0)**source['nu'])*RES/WIDTH/RE
                sr=s0*((leffr/lam0)**source['nu'])*RES/WIDTH/RE
                si=s0*((leffi/lam0)**source['nu'])*RES/WIDTH/RE
                sz=s0*((leffz/lam0)**source['nu'])*RES/WIDTH/RE
                sy=s0*((leffy/lam0)**source['nu'])*RES/WIDTH/RE
                
                dsu = makeGaussian(int(sy)*2, su, center=None)
                dsg = makeGaussian(int(sy)*2, sg, center=None)
                dsr = makeGaussian(int(sy)*2, sr, center=None)
                dsi = makeGaussian(int(sy)*2, si, center=None)
                dsz = makeGaussian(int(sy)*2, sz, center=None)
                dsy = makeGaussian(int(sy)*2, sy, center=None)
                gausslimit = float(len(dsy))
                print gausslimit,' gausslimit'
                
                if source['theta'] != 0:
                        dsu = irot(source['theta'],dsu)
                        dsg = irot(source['theta'],dsg)
                        dsr = irot(source['theta'],dsr)
                        dsi = irot(source['theta'],dsi)
                        dsz = irot(source['theta'],dsz)
                        dsy = irot(source['theta'],dsy)

        else:
                
                dsu = pyfits.open(outdir+'/custom/u.fits')[0].data
                dsg = pyfits.open(outdir+'/custom/g.fits')[0].data
                dsr = pyfits.open(outdir+'/custom/r.fits')[0].data
                dsi = pyfits.open(outdir+'/custom/i.fits')[0].data
                dsz = pyfits.open(outdir+'/custom/z.fits')[0].data
                dsy = pyfits.open(outdir+'/custom/y.fits')[0].data
                gausslimit = float(dsy.shape[0])
                print dsu.shape,dsg.shape,dsr.shape,dsi.shape,dsz.shape,dsy.shape
                print gausslimit

        print 'source has been set'




        #seed the random number generator
        #np.random.seed(123)
        np.random.seed()
	
	flag = 0
	lim = 100
	errlim = 0.05




	startconv = pytime.clock()
	newdataconvu = signal.fftconvolve(img,dsu)[limitx1:limitx2,limity1:limity2]
	newdataconvg = signal.fftconvolve(img,dsg)[limitx1:limitx2,limity1:limity2]
	newdataconvr = signal.fftconvolve(img,dsr)[limitx1:limitx2,limity1:limity2]
	newdataconvi = signal.fftconvolve(img,dsi)[limitx1:limitx2,limity1:limity2]
	newdataconvz = signal.fftconvolve(img,dsz)[limitx1:limitx2,limity1:limity2]
	newdataconvy = signal.fftconvolve(img,dsy)[limitx1:limitx2,limity1:limity2]
	print "convolution done"
	endconv = pytime.clock()
	

        startcurve = pytime.clock()
	for sample in range(int(output['numlc'])):
		print 'velocity correction prima: ',vcor[0]*DLS/((1+other['zl'])*DL),' Km/s'
		print 'velocity random: ',vel,' Km/s'
			
		daytolen = (vtot[sample])*86400.0*(1000.0)/RE*(RES/WIDTH)
		print 'length of the curve = ',lar1,' pixels'
		print 'angle = ',ang/pi*180,'º'
	
	
		
		#random initial values between actual possible range
		if ranx[sample]>=0.0:
			x0=x00[sample]*(imx-ranx[sample]) #np.random.uniform(gausslimit,10000-gausslimit-ranx)
			print 0,'< x0=',x0,'<',imx-ranx[sample]
		if ranx[sample]<0.0:
			x0=x00[sample]*(imx)-ranx[sample] #np.random.uniform(-ranx+gausslimit,10000-gausslimit)
			print 0-ranx[sample],'< x0=',x0,'<',imx
	
		if rany[sample]>=0.0:
			y0=y00[sample]*(imy-rany[sample]) #np.random.uniform(gausslimit,10000-gausslimit-rany)
			print 0,'< y0=',y0,'<',imy-rany[sample]
		if rany[sample]<0.0:
			y0=y00[sample]*(imy)-rany[sample] #np.random.uniform(-rany+gausslimit,10000-gausslimit)
			print 0-rany[sample],'< y0=',y0,'<',imy
	
		#coordinates of the final point of the curve
		x2 = x0+ranx[sample]
		y2 = y0+rany[sample]

		limitxi = int(np.floor(np.minimum(x0,x2)))
		limitxf = int(np.ceil(np.maximum(x0,x2))+1)
		limityi = int(np.floor(np.minimum(y0,y2)))
		limityf = int(np.ceil(np.maximum(y0,y2))+1)
	
		newdataconvu2 = newdataconvu[limitxi:limitxf,limityi:limityf]
		newdataconvg2 = newdataconvg[limitxi:limitxf,limityi:limityf]
		newdataconvr2 = newdataconvr[limitxi:limitxf,limityi:limityf]
		newdataconvi2 = newdataconvi[limitxi:limitxf,limityi:limityf]
		newdataconvz2 = newdataconvz[limitxi:limitxf,limityi:limityf]
		newdataconvy2 = newdataconvy[limitxi:limitxf,limityi:limityf]
		newdata2      = newdata[limitxi:limitxf,limityi:limityf]

		#inverted occurence, floor gives us back the missed decimal
		if ranx[sample]>=0:
			b=np.linspace(x0-np.floor(x0),x2-np.floor(x0),int(np.ceil(lar1[sample])+1))
		if ranx[sample]<0:
			b=np.linspace(x0-np.floor(x2),x2-np.floor(x2),int(np.ceil(lar1[sample])+1))
		if rany[sample]>=0:
			a=np.linspace(y0-np.floor(y0),y2-np.floor(y0),int(np.ceil(lar1[sample])+1))
		if rany[sample]<0:
			a=np.linspace(y0-np.floor(y2),y2-np.floor(y2),int(np.ceil(lar1[sample])+1))
	
		#interpolating
		int1su = ndimage.map_coordinates(newdataconvu2,np.vstack((b,a)))#,order=1)
		int1sg = ndimage.map_coordinates(newdataconvg2,np.vstack((b,a)))#,order=1)
		int1sr = ndimage.map_coordinates(newdataconvr2,np.vstack((b,a)))#,order=1)
		int1si = ndimage.map_coordinates(newdataconvi2,np.vstack((b,a)))#,order=1)
		int1sz = ndimage.map_coordinates(newdataconvz2,np.vstack((b,a)))#,order=1)
		int1sy = ndimage.map_coordinates(newdataconvy2,np.vstack((b,a)))#,order=1)
		print ('interpolation (done)')
	
		t = np.linspace(0,lar1[sample],int(np.ceil(lar1[sample])+1))
				
		fsu     = interp1d(t,int1su)
		int12su = fsu(ut*daytolen)
		fsg     = interp1d(t,int1sg)
		int12sg = fsg(gt*daytolen)
		fsr     = interp1d(t,int1sr)
		int12sr = fsr(rt*daytolen)
		fsi     = interp1d(t,int1si)
		int12si = fsi(it*daytolen)
		fsz     = interp1d(t,int1sz)
		int12sz = fsz(zt*daytolen)
		fsy     = interp1d(t,int1sy)
		int12sy = fsy(yt*daytolen)
		
#		print x0,y0,x2,y2,ranx[sample],rany[sample],lar1[sample]
#		print b
#		print a



#               ALL THE PLOTS ARE CREATED HERE
#
#		fig,(ax1,ax2,ax3)=plt.subplots(1,3)
#		fig.set_size_inches(18.5, 10.5)
#
#		#magnification image plot
#		ax1.imshow(np.log10(newdata2),cmap='gray',origin='lower')
#		
#		#source plot
#                if source['type'] == 'gaussian':
#                        circle1=plt.Circle(( (a[0]+a[-1])/2.0   ,(b[0]+b[-1])/2.0),sy,color='r')
#                        circle2=plt.Circle(( (a[0]+a[-1])/2.0   ,(b[0]+b[-1])/2.0),sz,color='m')
#                        circle3=plt.Circle(( (a[0]+a[-1])/2.0   ,(b[0]+b[-1])/2.0),si,color='y')
#                        circle4=plt.Circle(( (a[0]+a[-1])/2.0   ,(b[0]+b[-1])/2.0),sr,color='g')
#                        circle5=plt.Circle(( (a[0]+a[-1])/2.0   ,(b[0]+b[-1])/2.0),sg,color='b')
#                        circle6=plt.Circle(( (a[0]+a[-1])/2.0   ,(b[0]+b[-1])/2.0),su,color='c')
#                        ax1.add_patch(circle1)
#                        ax1.add_patch(circle2)
#                        ax1.add_patch(circle3)
#                        ax1.add_patch(circle4)
#                        ax1.add_patch(circle5)
#                        ax1.add_patch(circle6)
#                        ax1.arrow(a[0],b[0],a[-1]-a[0],b[-1]-b[0], head_length=lar1[sample]/10., head_width=lar1[sample]/10, fc='lime', ec='lime',length_includes_head=True)
#                else:
#                        ax1.arrow(a[0],b[0],a[-1]-a[0],b[-1]-b[0], head_length=lar1[sample]/10., head_width=lar1[sample]/10, fc='lime', ec='lime',length_includes_head=True)
#
#		#magnification image plot labels
#		ax1.set_xlabel('pixel')
#		ax1.set_ylabel('pixel')
#		#magnitude on y axis and pixels on x axis plot of the curve shown by the arrow
#		ax2.plot(t,-2.5*np.log10(int1su),'c',t,-2.5*np.log10(int1sg),'b',t,-2.5*np.log10(int1sr),'g',t,-2.5*np.log10(int1si),'y',t,-2.5*np.log10(int1sz),'m',t,-2.5*np.log10(int1sy),'r')
#		ax2.legend([r'$\sigma_u$', r'$\sigma_g$',r'$\sigma_r$',r'$\sigma_i$',r'$\sigma_z$',r'$\sigma_y$'], loc='upper left')
#		ax2.invert_yaxis()
#		ax2.set_ylabel(r'$\Delta$'+' [mag]')
#		ax2.set_xlabel('pixel')
#		#ax2.set_ylim([0.2,-1.4])
#	
#		ax3.errorbar(tmin+ut-tminplot,-2.5*np.log10(int12su),yerr=m52snr(-2.5*np.log10(int12su)+lsst['errbaseu'],um),fmt='c')
#		ax3.errorbar(tmin+gt-tminplot,-2.5*np.log10(int12sg),yerr=m52snr(-2.5*np.log10(int12sg)+lsst['errbaseg'],gm),fmt='b')
#		ax3.errorbar(tmin+rt-tminplot,-2.5*np.log10(int12sr),yerr=m52snr(-2.5*np.log10(int12sr)+lsst['errbaser'],rm),fmt='g')
#		ax3.errorbar(tmin+it-tminplot,-2.5*np.log10(int12si),yerr=m52snr(-2.5*np.log10(int12si)+lsst['errbasei'],im),fmt='y')
#		ax3.errorbar(tmin+zt-tminplot,-2.5*np.log10(int12sz),yerr=m52snr(-2.5*np.log10(int12sz)+lsst['errbasez'],zm),fmt='m')
#		ax3.errorbar(tmin+yt-tminplot,-2.5*np.log10(int12sy),yerr=m52snr(-2.5*np.log10(int12sy)+lsst['errbasey'],ym),fmt='r')
#		ax3.legend([r'$\sigma_u$', r'$\sigma_g$',r'$\sigma_r$',r'$\sigma_i$',r'$\sigma_z$',r'$\sigma_y$'], loc='upper left')
#		ax3.invert_yaxis()
#		ax3.set_ylabel(r'$\Delta$'+' [mag]')
#		ax3.set_xlabel('JD - 24'+str(tminplot))
#		ax3.set_ylim(ax2.get_ylim())
#		ax3.set_xlim(ax2.get_xlim()/daytolen+tmin-tminplot)	
#
#		if (output['write_image'] == True):
#			plt.savefig(outdir+'/output/plot_'+str(sample).zfill(4)+'.png')
#		if (output['get_plot'] == True):
#			plt.show()
#
#		plt.close('all')	

	
		
		#output tables with the light curve data (theoretical+filters, values, errors, etc)
		if output['full_data']==True:
                        outfunctions.original_output(outdir,sample,lsst,t,tmin,ut,gt,rt,it,zt,yt,int1su,int1sg,int1sr,int1si,int1sz,int1sy,int12su,int12sg,int12sr,int12si,int12sz,int12sy,um,gm,rm,im,zm,ym)
                if output['degraded_data']==True:
                        outfunctions.compressed_output_1(outdir,sample,lsst,daytolen,tmin,ut,gt,rt,it,zt,yt,int1su,int1sg,int1sr,int1si,int1sz,int1sy,int12su,int12sg,int12sr,int12si,int12sz,int12sy,um,gm,rm,im,zm,ym)



        endcurve = pytime.clock()
        endtime = pytime.clock()

        print 'Curve Time =', endcurve - startcurve
        print 'Convolution Time =', endconv - startconv
        print 'Total Time =', endtime - starttime

	return



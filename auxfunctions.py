import pyfits
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from math import *

import matplotlib.mlab as mlab
from scipy import signal
from scipy.fftpack import fft, fftshift
import scipy.misc
from astropy.io import fits
import sys , os
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
from scipy import ndimage

#some constants
C=2.99792458*(10**8)
Ms=1.9891*(10**30)
G=6.67384/(10**11)
Parsec=3.08*(10**16)

def getKey(item):
	return item[0]

def effe(kappa,gamma,s):
	keffe=(1-s)*kappa/(1-s*kappa)
	geffe=gamma/(1-kappa*s)
	return keffe,geffe

def gv_readimg(path_to_maps,mapid):
        f = open(path_to_maps+mapid+'/mapmeta.dat')
        l = f.readlines()
        dum = np.fromstring(l[0],dtype=float, sep=' ')
        AVGMAG = dum[0]
        NRAY   = dum[1]
        RES    = int(l[1])
        WIDTH  = float(l[2])
        fac    = abs(AVGMAG/NRAY)

        f   = open(path_to_maps+mapid+'/map.bin','rb')
        map = np.fromfile(f,'i',-1,"")

        twomap = np.ndarray(shape=(RES,RES),dtype=float)
        for i in range(0,RES):
                for j in range(0,RES):
                        twomap[i,j] = map[i*RES+j]*fac


        return twomap,RES,WIDTH


def readimg(kappa,gamma,s,tol):
	data = np.genfromtxt("/data3/GERLUMPH/NEW/FINAL/table.dat",dtype=None)
	print 'there are ',len(data[:]),' images to use'
	diffk = np.zeros(len(data[:]))
	diffg = np.zeros(len(data[:]))
	diffs = np.zeros(len(data[:]))
	for i in range(len(data[:])):
		diffk[i] = (kappa - data[i][0])**2
	for i in range(len(data[:])):
		diffg[i] = (gamma - data[i][1])**2
	for i in range(len(data[:])):
		diffs[i] = (s - data[i][2])**2
		diff = diffk + diffg + diffs
		diff = diff.tolist()
		im = diff.index(min(diff))
	print data
	if diffk[im]+diffg[im]+diffs[im]> tol:
		keffe,geffe=effe(kappa,gamma,s)
		print "kappa = ",kappa,"gamma = ",gamma
		print "kappaeff = ",keffe,"gammaeff = ",geffe
		print 'using effective values'
		diffk = np.zeros(len(data[:]))
		diffg = np.zeros(len(data[:]))
		diffs = np.zeros(len(data[:]))
		for i in range(len(data[:])):
			diffk[i] = (keffe - data[i][3])**2
		for i in range(len(data[:])):
			diffg[i] = (geffe - data[i][4])**2
		for i in range(len(data[:])):
			diffs[i] = (s - data[i][2])**2
		diff = diffk + diffg + diffs
		diff = diff.tolist()
		im = diff.index(min(diff))
		print 'using image ',data[im][5]
		print 'kappaeff=',data[im][3]
		print 'gammaeff=',data[im][4]
		img1 = pyfits.open("/data3/GERLUMPH/NEW/FINAL/"+data[im][5])
		head1= img1[0].header
		RES=head1['RES']
		WIDTH=head1['WIDTH']
		AVGMAG=1/((1-kappa)*(1-kappa)-(gamma*gamma))
		NRAY=head1['NRAY']
		#print AVGMAG,1/((1-data[im][3])**2-(data[im][4])**2)
		return img1,RES,WIDTH,AVGMAG,NRAY
	else:
		print 'using image ',data[im][5]
		img1 = pyfits.open("/data3/GERLUMPH/NEW/FINAL/"+data[im][5])
		head1= img1[0].header
		RES=head1['RES']
		WIDTH=head1['WIDTH']
		AVGMAG=head1['AVGMAG']
		NRAY=head1['NRAY']
		return img1,RES,WIDTH,AVGMAG,NRAY

#Functions of calculus luminar Distance, angular distance, Einstein radius

def Dis(zs,zl,mass):
	cosmo = FlatLambdaCDM(H0=67.7 * u.km / u.s / u.Mpc, Om0=0.309)
	Dl=cosmo.luminosity_distance(zl).value
	Ds=cosmo.luminosity_distance(zs).value
	DS=Parsec*(10**6)*cosmo.angular_diameter_distance(zs).value
	DL=Parsec*(10**6)*cosmo.angular_diameter_distance(zl).value
	DLS=Parsec*(10**6)*cosmo.angular_diameter_distance_z1z2(zl,zs).value
	RE=DS*sqrt((4*G*Ms*mass/(C*C))*(DLS/(DL*DS)))
	return RE,DL,DS,DLS

def veco(ra,dec):
	ra = ra*15.0
	eq2ga = SkyCoord(ra=ra,dec=dec,unit='degree').galactic
	b,l=np.radians(np.array(eq2ga.b*u.deg)),np.radians(np.array(eq2ga.l*u.deg))
	Vapex=387.
	bapex=(90-48.26)*pi/180.0
	lapex=264.14*pi/180.0
	b=90.0*pi/180.0-b
	#return Vapex*(np.sin(b)*np.sin(bapex*pi/180.)+np.cos(b)*np.cos(bapex*pi/180.)*np.cos(l-(lapex*pi/180.)))
	cor = np.array([np.cos(l)*np.sin(b),np.sin(l)*np.sin(b),np.cos(b)])
	print cor.shape,'cor'
	vcmb = np.array([Vapex*np.cos(lapex)*np.sin(bapex),Vapex*np.sin(lapex)*np.sin(bapex),Vapex*np.cos(bapex)])
	print vcmb.shape,'vcmb'
	alpha = vcmb.dot(cor)
	alpha2 = alpha*cor
	print alpha2.shape,'alpha2'
	r = - alpha2.T + vcmb 	
	print r#.shape,'r'
	if r.size >3:
		norm = np.linalg.norm(r,axis=1)
	else:
		norm = np.linalg.norm(r)
	print norm, norm.shape,'norm'
	if r.size >3:
		b, l = np.arccos(r[:,2]/norm), np.arctan(r[:,1]/r[:,0])
	else:
		b, l = np.arccos(r[2]/norm), np.arctan(r[1]/r[0])
	ga2eq = SkyCoord(frame='galactic',l=l,b=b,unit=(u.degree,u.degree)).icrs
	b = b + 90
	dec2,ra2=np.array(ga2eq.dec*u.rad),np.array(ga2eq.ra*u.rad)
	return np.array(norm),np.array(dec2),np.array(ra2)
	


def vecosimp(ra,dec):
	ra = ra*15.0
	eq2ga = SkyCoord(ra=ra,dec=dec,unit='degree').galactic
	b,l=np.radians(eq2ga.b).value,np.radians(eq2ga.l).value
	Vapex=369.
	bapex=(90-48.4)*pi/180.0
	lapex=264.4*pi/180.0
	print b,l
	b=90.0*pi/180.0-b
	#return Vapex*(np.sin(b)*np.sin(bapex*pi/180.)+np.cos(b)*np.cos(bapex*pi/180.)*np.cos(l-(lapex*pi/180.)))
	cor = np.array([np.cos(l)*np.sin(b),np.sin(l)*np.sin(b),np.cos(b)])
	print cor.shape,'cor'
	vcmb = np.array([Vapex*np.cos(lapex)*np.sin(bapex),Vapex*np.sin(lapex)*np.sin(bapex),Vapex*np.cos(bapex)])
	print vcmb.shape,'vcmb'
	alpha = vcmb.dot(cor)
	alpha2 = alpha*cor
	print alpha2.shape,'alpha2'
	r = - alpha2.T + vcmb 	
	print r#.shape,'r'
	norm = np.linalg.norm(r)
	print norm, norm.shape,'norm'
	r, phi, theta = norm, np.arccos(r[2]/norm), np.arctan(r[1]/r[0])
	print r,phi*180/pi,theta*180/pi
	print sqrt(phi**2+theta**2)
	b = 90 - b
	ga2eq = SkyCoord(frame='galactic',l=l,b=b,unit=(u.degree,u.degree)).icrs
	dec2,ra2=(ga2eq.dec*u.rad).value,(ga2eq.ra*u.rad).value
	return norm,dec2,ra2

def makeGaussian(size, fwhm, center=None):
	"""Make a square gaussian kernel.
	size is the length of a side of the square
	fwhm is full-width-half-maximum, which
	can be thought of as an effective radius."""

	x = np.arange(0, size, 1, float)
	y = x[:,np.newaxis]
	if center is None:
		x0 = y0 = size // 2
	else:
		x0 = center[0]
		y0 = center[1]
	gauss=np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)
	return gauss/np.sum(gauss)

def m52snr(m, m5):
	"""find the SNR for a star of magnitude m obsreved
	under conditions of 5-sigma limiting depth m5.  This assumes
	Gaussianity and might not be strictly true in bluer filters.
	See table 2 and eq 5 in astroph/0805.2366 """
	snr = 5.*10.**(-0.4*(m-m5))
	lc_err = 2.5*np.log10(1.+1./snr)
	return lc_err

def irot(angulo,ds):

	#datos imagenes
	dat = ds
	fil = dat.shape[0]
	#angulo inclinacion
	ang=angulo
	#angulo rotacion
	# determina si los filtros son pares

	#valor de que se va a cortar + imagen amplificada 10 veces(por el angulo)

	v=int(round(fil*((1/abs(cos(ang*pi/180.0)))-1)*abs(cos(ang*pi/180)/2)))*10
	v1=v
	if ang==90.0:
		v=(fil)*10/2
		v1=fil*10/2-1

	#creacion matriz nueva, imagen amplificada
	cubo=np.zeros((abs(fil*10-2*v),fil*10),np.int64)
	#creacion matriz para ang 90 o 270
	cuboa=np.zeros((fil*10,fil*10),np.int64)
	#definicion de las posiciones de la matriz para la interpolcion
	x=np.linspace(0,fil-1,fil)
	y=np.linspace(0,fil-1,fil)
	z=dat

	#definicion de las posiciones para la matriz ampliada sin rotar
	xa=np.linspace(x.min(),x.max(),fil*10)
	ya=np.linspace(x.min(),x.max(),fil*10)
	#definicion de las posiciones para la matriz rotada
	yy=np.linspace(x.min(),x.max(),abs(fil*10-2*v))
	xx=np.linspace(x.min(),x.max(),fil*10)
	#suma de las cuentas matriz original
	sumaz=dat[:,:].sum()
	#definicion de la funcion interpolacion
	newcubo = interpolate.interp2d(x,y,z, kind='linear')
	#matriz interpolada matriz ampliada no rotada
	cuboa=abs(newcubo(xa,ya))
	rubixa=np.array(cuboa)
	zz=rubixa.sum(axis=0)
	zzz=dat.sum(axis=0)


	#rotacion para la matriz cuando es 90 o 270, es la suma por columnas
	if ang==90 or ang==270:
		cubo=rubixa.sum(axis=1)
	#matriz interpolada , ampliada y rotada
	else:
		cubo = abs(newcubo(xx,yy)/cos(ang*pi/180.0))
	#definicion de matriz 0 ampliada y rotada

	cubofinal=[]

	#completacion matriz, matriz ampliada y rotada

	for i in range (v1):
		cuborow=[]
		for l in range (fil*10):
			cuborow.append(0)
		cubofinal.append(cuborow)
	if cubo.ndim==1:
		for e in range(1):
			cuborowf=[]
			for h in range(fil*10):
				cuborowf.append(cubo[h])
			cubofinal.append(cuborowf)
	else:
		for m in range (abs(fil*10-2*v)):
			cuborowf=[]
			for n in range(fil*10):
				cuborowf.append(cubo[m][n])
			cubofinal.append(cuborowf)

	for o in range(v):
		cuborow2=[]
		for p in range(fil*10):
			cuborow2.append(0)
		cubofinal.append(cuborow2)

	#matriz final definida como array
	rubix=np.array(cubofinal)

	if cubo.ndim==1:
		sumar=rubix.sum()

	else:
		sumar=rubix[:,:].sum()




	#proceso para volver la matriz a su tamano original

	cuboo=np.zeros((fil,fil),np.int64)

	xo=np.linspace(0,fil*10.0-1,fil*10.0)
	yo=np.linspace(0,fil*10.0-1,fil*10.0)

	yyo=np.linspace(yo.min(),yo.max(),fil)
	xxo=np.linspace(xo.min(),xo.max(),fil)



	newcuboo = interpolate.interp2d(xo,yo,rubix, kind='linear')



	cuboo=abs(newcuboo(xxo,yyo))

	cubooo=cuboo/cuboo[:,:].sum()

	return cubooo 

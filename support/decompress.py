import numpy as np
import os
import sys
import glob


def get_values(fname,dtype,offset,size,mmin,mmax):
    d = (mmax-mmin)/255.0
    def f(x):
        return (d*x+mmin)
    g = np.vectorize(f)    

    fh = open(fname,'rb')
    fh.seek(offset,os.SEEK_SET)
    a = np.fromfile(fh,np.dtype(dtype),count=size)
    return g(a)


def get_days(fname,dtype,offset,size,daytolen,tmin):
    def f(x):
        return (tmin+x/daytolen)
    g = np.vectorize(f)

    fh = open(fname,'rb')
    fh.seek(offset,os.SEEK_SET)
    a = np.fromfile(fh,np.dtype(dtype),count=size)
    return g(a)


# Beginning of the program

indir  = sys.argv[1]
if len(sys.argv) > 2:
    outdir = sys.argv[2]
else:
    outdir = indir
    

indices = []
dum = sorted(glob.glob(indir+'param_*.dat'))
for d in dum:
    head,tail = os.path.split(d)
    name,ext  = os.path.splitext(tail)
    tmp       = name.split('_')
    indices.append(tmp[1])




for ind in indices:
    pars  = np.loadtxt(indir+'param_'+ind+'.dat',skiprows=1)
    sizes = pars[0].astype(int)
    mins  = pars[1].astype(float)
    maxs  = pars[2].astype(float)

    with open(indir+'param_'+ind+'.dat','r') as f:
        line = f.readline()
    dum = line[2:].split(',')
    daytolen = float(dum[0])
    tmin     = float(dum[1])


    Nt   = sizes[6]
    mint = mins[6]
    maxt = maxs[6]
    u = get_values(indir+'compt_'+ind+'.bin','uint8',0*Nt,Nt,mint,maxt)
    g = get_values(indir+'compt_'+ind+'.bin','uint8',1*Nt,Nt,mint,maxt)
    r = get_values(indir+'compt_'+ind+'.bin','uint8',2*Nt,Nt,mint,maxt)
    i = get_values(indir+'compt_'+ind+'.bin','uint8',3*Nt,Nt,mint,maxt)
    z = get_values(indir+'compt_'+ind+'.bin','uint8',4*Nt,Nt,mint,maxt)
    y = get_values(indir+'compt_'+ind+'.bin','uint8',5*Nt,Nt,mint,maxt)

    theo = np.array([range(sizes[6]),u,g,r,i,z,y]).T
    np.savetxt(outdir+'tablet_'+ind+'.dat',theo,delimiter="	",header="t	u	g	r	i	z	y")



    uval = get_values(indir+'compf_'+ind+'.bin','uint8',0,sizes[0],mint,maxt)
    gval = get_values(indir+'compf_'+ind+'.bin','uint8',sizes[0],sizes[1],mint,maxt)
    rval = get_values(indir+'compf_'+ind+'.bin','uint8',sizes[0]+sizes[1],sizes[2],mint,maxt)
    ival = get_values(indir+'compf_'+ind+'.bin','uint8',sizes[0]+sizes[1]+sizes[2],sizes[3],mint,maxt)
    zval = get_values(indir+'compf_'+ind+'.bin','uint8',sizes[0]+sizes[1]+sizes[2]+sizes[3],sizes[4],mint,maxt)
    yval = get_values(indir+'compf_'+ind+'.bin','uint8',sizes[0]+sizes[1]+sizes[2]+sizes[3]+sizes[4],sizes[5],mint,maxt)

    uerr = get_values(indir+'compe_'+ind+'.bin','uint8',0,sizes[0],mins[0],maxs[0])
    gerr = get_values(indir+'compe_'+ind+'.bin','uint8',sizes[0],sizes[1],mins[1],maxs[1])
    rerr = get_values(indir+'compe_'+ind+'.bin','uint8',sizes[0]+sizes[1],sizes[2],mins[2],maxs[2])
    ierr = get_values(indir+'compe_'+ind+'.bin','uint8',sizes[0]+sizes[1]+sizes[2],sizes[3],mins[3],maxs[3])
    zerr = get_values(indir+'compe_'+ind+'.bin','uint8',sizes[0]+sizes[1]+sizes[2]+sizes[3],sizes[4],mins[4],maxs[4])
    yerr = get_values(indir+'compe_'+ind+'.bin','uint8',sizes[0]+sizes[1]+sizes[2]+sizes[3]+sizes[4],sizes[5],mins[5],maxs[5])
    
    uday = get_days(indir+'compl_'+ind+'.bin','uint16',0,sizes[0],daytolen,tmin)
    gday = get_days(indir+'compl_'+ind+'.bin','uint16',2*sizes[0],sizes[1],daytolen,tmin)
    rday = get_days(indir+'compl_'+ind+'.bin','uint16',2*(sizes[0]+sizes[1]),sizes[2],daytolen,tmin)
    iday = get_days(indir+'compl_'+ind+'.bin','uint16',2*(sizes[0]+sizes[1]+sizes[2]),sizes[3],daytolen,tmin)
    zday = get_days(indir+'compl_'+ind+'.bin','uint16',2*(sizes[0]+sizes[1]+sizes[2]+sizes[3]),sizes[4],daytolen,tmin)
    yday = get_days(indir+'compl_'+ind+'.bin','uint16',2*(sizes[0]+sizes[1]+sizes[2]+sizes[3]+sizes[4]),sizes[5],daytolen,tmin)
    
    ufilter = np.array([uday,uval,uerr]).T
    gfilter = np.array([gday,gval,gerr]).T
    rfilter = np.array([rday,rval,rerr]).T
    ifilter = np.array([iday,ival,ierr]).T
    zfilter = np.array([zday,zval,zerr]).T
    yfilter = np.array([yday,yval,yerr]).T

    np.savetxt(outdir+'tableu_'+ind+'.dat',ufilter,delimiter="	",header="JD	u	uerr",fmt='%10.3f')
    np.savetxt(outdir+'tableg_'+ind+'.dat',gfilter,delimiter="	",header="JD	g	gerr",fmt='%10.3f')
    np.savetxt(outdir+'tabler_'+ind+'.dat',rfilter,delimiter="	",header="JD	r	rerr",fmt='%10.3f')
    np.savetxt(outdir+'tablei_'+ind+'.dat',ifilter,delimiter="	",header="JD	i	ierr",fmt='%10.3f')
    np.savetxt(outdir+'tablez_'+ind+'.dat',zfilter,delimiter="	",header="JD	z	zerr",fmt='%10.3f')
    np.savetxt(outdir+'tabley_'+ind+'.dat',yfilter,delimiter="	",header="JD	y	yerr",fmt='%10.3f')





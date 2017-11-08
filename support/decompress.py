import numpy as np
import os
import sys
import glob


def get_values(fname,dtype,offset,size,mmin,mmax):
    d = (mmax-mmin)/float(np.iinfo(np.dtype(dtype)).max)
    def f(x):
        return (d*x+mmin)
    g = np.vectorize(f)    

    fh = open(fname,'rb')
    fh.seek(np.dtype(dtype).itemsize*offset,0)
    a = np.fromfile(fh,np.dtype(dtype),count=size)
    fh.close()
    return g(a)


# Beginning of the program
indir  = sys.argv[1]
if len(sys.argv) > 2:
    outdir = sys.argv[2]
else:
    outdir = indir

indices = []
dum = sorted(glob.glob(indir+'comp_p_*.dat'))
for d in dum:
    head,tail = os.path.split(d)
    name,ext  = os.path.splitext(tail)
    tmp       = name.split('_')
    indices.append(tmp[2])



for ind in indices:
    with open(indir+'comp_p_'+ind+'.dat') as f:
        content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
        content = [x.strip() for x in content] 
    sizes = map(int,content[0].split(' '))

    dum = content[1].split(' ')
    mytypes = []
    for item in dum:
        item = item.strip('\'')
        if item == 't':
            mytypes.append('uint16')
        elif item == 'h':
            mytypes.append('uint8')
    dtypem = mytypes[0]
    dtypet = mytypes[1]
    dtypee = mytypes[2]
    
    dum = content[2].split(' ')
    minm = float(dum[0])
    maxm = float(dum[1])
    dum = content[3].split(' ')
    mint = float(dum[0])
    maxt = float(dum[1])
    dum = content[4].split(' ')
    mine = float(dum[0])
    maxe = float(dum[1])


    Nt   = sizes[6]/6
    u = get_values(indir+'comp_f_'+ind+'.bin',dtypem,0*Nt,Nt,minm,maxm)
    g = get_values(indir+'comp_f_'+ind+'.bin',dtypem,1*Nt,Nt,minm,maxm)
    r = get_values(indir+'comp_f_'+ind+'.bin',dtypem,2*Nt,Nt,minm,maxm)
    i = get_values(indir+'comp_f_'+ind+'.bin',dtypem,3*Nt,Nt,minm,maxm)
    z = get_values(indir+'comp_f_'+ind+'.bin',dtypem,4*Nt,Nt,minm,maxm)
    y = get_values(indir+'comp_f_'+ind+'.bin',dtypem,5*Nt,Nt,minm,maxm)

    theo = np.array([range(Nt),u,g,r,i,z,y]).T
    #np.savetxt(outdir+'tablet_'+ind+'.dat',theo,delimiter="	",header="t	u	g	r	i	z	y")
    np.savetxt(outdir+'tablet_'+ind+'.dat',theo,delimiter=" ")




    uval = get_values(indir+'comp_m_'+ind+'.bin',dtypem,0,sizes[0],minm,maxm)
    gval = get_values(indir+'comp_m_'+ind+'.bin',dtypem,sizes[0],sizes[1],minm,maxm)
    rval = get_values(indir+'comp_m_'+ind+'.bin',dtypem,sizes[0]+sizes[1],sizes[2],minm,maxm)
    ival = get_values(indir+'comp_m_'+ind+'.bin',dtypem,sizes[0]+sizes[1]+sizes[2],sizes[3],minm,maxm)
    zval = get_values(indir+'comp_m_'+ind+'.bin',dtypem,sizes[0]+sizes[1]+sizes[2]+sizes[3],sizes[4],minm,maxm)
    yval = get_values(indir+'comp_m_'+ind+'.bin',dtypem,sizes[0]+sizes[1]+sizes[2]+sizes[3]+sizes[4],sizes[5],minm,maxm)

    uerr = get_values(indir+'comp_e_'+ind+'.bin',dtypee,0,sizes[0],mine,maxe)
    gerr = get_values(indir+'comp_e_'+ind+'.bin',dtypee,sizes[0],sizes[1],mine,maxe)
    rerr = get_values(indir+'comp_e_'+ind+'.bin',dtypee,sizes[0]+sizes[1],sizes[2],mine,maxe)
    ierr = get_values(indir+'comp_e_'+ind+'.bin',dtypee,sizes[0]+sizes[1]+sizes[2],sizes[3],mine,maxe)
    zerr = get_values(indir+'comp_e_'+ind+'.bin',dtypee,sizes[0]+sizes[1]+sizes[2]+sizes[3],sizes[4],mine,maxe)
    yerr = get_values(indir+'comp_e_'+ind+'.bin',dtypee,sizes[0]+sizes[1]+sizes[2]+sizes[3]+sizes[4],sizes[5],mine,maxe)
    
    uday = get_values(indir+'comp_t_'+ind+'.bin',dtypet,0,sizes[0],mint,maxt)
    gday = get_values(indir+'comp_t_'+ind+'.bin',dtypet,sizes[0],sizes[1],mint,maxt)
    rday = get_values(indir+'comp_t_'+ind+'.bin',dtypet,sizes[0]+sizes[1],sizes[2],mint,maxt)
    iday = get_values(indir+'comp_t_'+ind+'.bin',dtypet,sizes[0]+sizes[1]+sizes[2],sizes[3],mint,maxt)
    zday = get_values(indir+'comp_t_'+ind+'.bin',dtypet,sizes[0]+sizes[1]+sizes[2]+sizes[3],sizes[4],mint,maxt)
    yday = get_values(indir+'comp_t_'+ind+'.bin',dtypet,sizes[0]+sizes[1]+sizes[2]+sizes[3]+sizes[4],sizes[5],mint,maxt)
    
    ufilter = np.array([uday,uval,uerr]).T
    gfilter = np.array([gday,gval,gerr]).T
    rfilter = np.array([rday,rval,rerr]).T
    ifilter = np.array([iday,ival,ierr]).T
    zfilter = np.array([zday,zval,zerr]).T
    yfilter = np.array([yday,yval,yerr]).T

    np.savetxt(outdir+'tableu_'+ind+'.dat',ufilter,delimiter=" ",fmt='%12.6e')
    np.savetxt(outdir+'tableg_'+ind+'.dat',gfilter,delimiter=" ",fmt='%12.6e')
    np.savetxt(outdir+'tabler_'+ind+'.dat',rfilter,delimiter=" ",fmt='%12.6e')
    np.savetxt(outdir+'tablei_'+ind+'.dat',ifilter,delimiter=" ",fmt='%12.6e')
    np.savetxt(outdir+'tablez_'+ind+'.dat',zfilter,delimiter=" ",fmt='%12.6e')
    np.savetxt(outdir+'tabley_'+ind+'.dat',yfilter,delimiter=" ",fmt='%12.6e')




import numpy as np
import array
import math
from auxfunctions import m52snr

def original_output(outdir,sample,lsst,t,tmin,ut,gt,rt,it,zt,yt,int1su,int1sg,int1sr,int1si,int1sz,int1sy,int12su,int12sg,int12sr,int12si,int12sz,int12sy,um,gm,rm,im,zm,ym):
    lc_count = str(sample).zfill(5)

    tablet = np.array([t,lsst['errbaseu']-2.5*np.log10(int1su),lsst['errbaseg']-2.5*np.log10(int1sg),lsst['errbaser']-2.5*np.log10(int1sr),lsst['errbasei']-2.5*np.log10(int1si),lsst['errbasez']-2.5*np.log10(int1sz),lsst['errbasey']-2.5*np.log10(int1sy)]).T
    tableu = np.array([tmin+ut,lsst['errbaseu']-2.5*np.log10(int12su),m52snr(-2.5*np.log10(int12su)+lsst['errbaseu'],um)]).T
    tableg = np.array([tmin+gt,lsst['errbaseg']-2.5*np.log10(int12sg),m52snr(-2.5*np.log10(int12sg)+lsst['errbaseg'],gm)]).T
    tabler = np.array([tmin+rt,lsst['errbaser']-2.5*np.log10(int12sr),m52snr(-2.5*np.log10(int12sr)+lsst['errbaser'],rm)]).T
    tablei = np.array([tmin+it,lsst['errbasei']-2.5*np.log10(int12si),m52snr(-2.5*np.log10(int12si)+lsst['errbasei'],im)]).T
    tablez = np.array([tmin+zt,lsst['errbasez']-2.5*np.log10(int12sz),m52snr(-2.5*np.log10(int12sz)+lsst['errbasez'],zm)]).T
    tabley = np.array([tmin+yt,lsst['errbasey']-2.5*np.log10(int12sy),m52snr(-2.5*np.log10(int12sy)+lsst['errbasey'],ym)]).T
	
    np.savetxt(outdir+'/output/tablet_'+lc_count+'.dat',tablet,delimiter="	",header="t	u	g	r	i	z	y")
    np.savetxt(outdir+'/output/tableu_'+lc_count+'.dat',tableu,delimiter="	",header="JD	u	uerr",fmt='%10.3f')
    np.savetxt(outdir+'/output/tableg_'+lc_count+'.dat',tableg,delimiter="	",header="JD	g	gerr",fmt='%10.3f')
    np.savetxt(outdir+'/output/tabler_'+lc_count+'.dat',tabler,delimiter="	",header="JD	r	rerr",fmt='%10.3f')
    np.savetxt(outdir+'/output/tablei_'+lc_count+'.dat',tablei,delimiter="	",header="JD	i	ierr",fmt='%10.3f')
    np.savetxt(outdir+'/output/tablez_'+lc_count+'.dat',tablez,delimiter="	",header="JD	z	zerr",fmt='%10.3f')
    np.savetxt(outdir+'/output/tabley_'+lc_count+'.dat',tabley,delimiter="	",header="JD	y	yerr",fmt='%10.3f')

    return






def get_errors(values,lsst_error,m):
    err     = m52snr(-2.5*np.log10(values)+lsst_error,m)
    err_min = np.amin(err)
    err_max = np.amax(err)
    err_d   = (err_max-err_min)/255.0
    def ee(x):
        return int(math.floor((x-err_min)/err_d))
    escale = np.vectorize(ee)
    byte_err = array.array('B',escale(err)).tostring()

    return err_min,err_max,byte_err



def compressed_output_1(outdir,sample,lsst,daytolen,tmin,ut,gt,rt,it,zt,yt,int1su,int1sg,int1sr,int1si,int1sz,int1sy,int12su,int12sg,int12sr,int12si,int12sz,int12sy,um,gm,rm,im,zm,ym):
    lc_count = str(sample).zfill(5)



    # Theoretical: get scaled values and convert to bytes
    theo = np.concatenate((lsst['errbaseu']-2.5*np.log10(int1su),lsst['errbaseg']-2.5*np.log10(int1sg),lsst['errbaser']-2.5*np.log10(int1sr),lsst['errbasei']-2.5*np.log10(int1si),lsst['errbasez']-2.5*np.log10(int1sz),lsst['errbasey']-2.5*np.log10(int1sy)),axis=0)
    theo_n   = len(int1su)
    theo_min = np.amin(theo)
    theo_max = np.amax(theo)
    theo_d   = (theo_max-theo_min)/255.0
    def f(x):
        return int(math.floor((x-theo_min)/theo_d))
    scale = np.vectorize(f)
    byte_theo = array.array('B',scale(theo)).tostring()
    with open(outdir+'/output/compt_'+lc_count+'.bin','wb') as f:
        f.write(byte_theo)



    # Filters: get scaled values and convert to bytes
    filters_n = [len(int12su),len(int12sg),len(int12sr),len(int12si),len(int12sz),len(int12sy)]
    filters = np.concatenate((lsst['errbaseu']-2.5*np.log10(int12su),lsst['errbaseg']-2.5*np.log10(int12sg),lsst['errbaser']-2.5*np.log10(int12sr),lsst['errbasei']-2.5*np.log10(int12si),lsst['errbasez']-2.5*np.log10(int12sz),lsst['errbasey']-2.5*np.log10(int12sy)),axis=0)
    byte_filters =  array.array('B',scale(filters)).tostring()
    with open(outdir+'/output/compf_'+lc_count+'.bin','wb') as f:
        f.write(byte_filters)



    # Locations: 
    locations = np.concatenate(((ut*daytolen).astype(np.dtype('uint16')),(gt*daytolen).astype(np.dtype('uint16')),(rt*daytolen).astype(np.dtype('uint16')),(it*daytolen).astype(np.dtype('uint16')),(zt*daytolen).astype(np.dtype('uint16')),(yt*daytolen).astype(np.dtype('uint16'))),axis=0)
    byte_locations =  array.array('H',locations).tostring()
    with open(outdir+'/output/compl_'+lc_count+'.bin','wb') as f:
        f.write(byte_locations)



    # Errors:
    err_min  = np.zeros(6)
    err_max  = np.zeros(6)
    err_min[0],err_max[0],u_byte_err = get_errors(int12su,lsst['errbaseu'],um)
    err_min[1],err_max[1],g_byte_err = get_errors(int12sg,lsst['errbaseg'],gm)
    err_min[2],err_max[2],r_byte_err = get_errors(int12sr,lsst['errbaser'],rm)
    err_min[3],err_max[3],i_byte_err = get_errors(int12si,lsst['errbasei'],im)
    err_min[4],err_max[4],z_byte_err = get_errors(int12sz,lsst['errbasez'],zm)
    err_min[5],err_max[5],y_byte_err = get_errors(int12sy,lsst['errbasey'],ym)
    with open(outdir+'/output/compe_'+lc_count+'.bin','wb') as f:
        f.write(u_byte_err)
        f.write(g_byte_err)
        f.write(r_byte_err)
        f.write(i_byte_err)
        f.write(z_byte_err)
        f.write(y_byte_err)



    # Parameters:
    table_dum = np.array([np.append(filters_n,theo_n),np.append(err_min,theo_min),np.append(err_max,theo_max)])
    np.savetxt(outdir+'/output/param_'+lc_count+'.dat',table_dum,delimiter="	",header=str(daytolen)+','+str(tmin))

    
    return

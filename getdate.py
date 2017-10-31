# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import lsst.sims.maf.db as db
import lsst.sims.maf.utils as utils
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
import asciitable

def getdate(outdir,path_to_dbfile,runName,ra,dec):
  outDir = outdir+'/dates'
  #runName = 'minion_1016'
  #dbFile = 'enigma_1189_sqlite.db'
  dbFile = path_to_dbfile + runName + '_sqlite.db'
  opsimdb = utils.connectOpsimDb(dbFile)
  resultsDb = db.ResultsDb(outDir=outDir)

  filters = ['u','g','r','i','z','y']
  #colors={'u':'cyan','g':'g','r':'y','i':'r','z':'m', 'y':'k'}
  ra=ra*15.0
  #ra = np.radians(ra)
  #dec = np.radians(dec)
  # SNR limit (Don't use points below this limit)
  snrLimit = 5.
  # Demand this many points above SNR limit before plotting LC
  nPtsLimit = 6
  # The pass metric just passes data straight through.
  metric = metrics.PassMetric(cols=['filter','fiveSigmaDepth','expMJD'])
  slicer = slicers.UserPointsSlicer(ra,dec,lonCol='ditheredRA',latCol='ditheredDec')
  sql = ''
  bundle = metricBundles.MetricBundle(metric,slicer,sql)
  bg =  metricBundles.MetricBundleGroup({0:bundle}, opsimdb, outDir=outDir, resultsDb=resultsDb)
  bg.runAll()
  bundle.metricValues.data[0]['filter']


  print '%i Observations total at this point (All SNR levels)' % bundle.metricValues.data[0].size
  for fname in filters:
    good = np.where(bundle.metricValues.data[0]['filter'] == fname)
    print '%i Observations in %s' % (good[0].size, fname)
    mags=bundle.metricValues.data[0]['expMJD'][good[0]],bundle.metricValues.data[0]['fiveSigmaDepth'][good[0]]
    mags=np.array(mags)
    mags=np.transpose(mags)
    mags=mags[mags[:,0].argsort()]
    asciitable.write({'JD': mags[:,0],'5sigmadepth': mags[:,1]},outDir+'/'+fname+".dat",names=['JD','5sigmadepth'])

    
# return mags#'Fin'



#read dates
  u = np.loadtxt(outdir+'/dates/u.dat',skiprows=1)
  g = np.loadtxt(outdir+'/dates/g.dat',skiprows=1)
  r = np.loadtxt(outdir+'/dates/r.dat',skiprows=1)
  i = np.loadtxt(outdir+'/dates/i.dat',skiprows=1)
  z = np.loadtxt(outdir+'/dates/z.dat',skiprows=1)
  y = np.loadtxt(outdir+'/dates/y.dat',skiprows=1)
  zz = u.tolist()+g.tolist()+r.tolist()+i.tolist()+z.tolist()+y.tolist()
  zz = sorted(zz, key=getKey)
  
  if len(u) == 0:
    u = np.empty((1,2))
    u[:] = np.NAN
  if len(g) == 0:
    g = np.empty((1,2))
    g[:] = np.NAN
  if len(r) == 0:
    r = np.empty((1,2))
    r[:] = np.NAN
  if len(i) == 0:
    i = np.empty((1,2))
    i[:] = np.NAN
  if len(z) == 0:
    z = np.empty((1,2))
    z[:] = np.NAN
  if len(y) == 0:
    y = np.empty((1,2))
    y[:] = np.NAN
	
  days = lsst['years']*365.25

  t1 = np.asarray(zz)[:,0]
  tmin = min(t1)
  tminplot = np.floor(tmin/1000.0)*1000.0

  ut,gt,rt,it,zt,yt = (u[:,0]-tmin),g[:,0]-tmin,r[:,0]-tmin,i[:,0]-tmin,z[:,0]-tmin,y[:,0]-tmin
  um,gm,rm,im,zm,ym = u[:,1],g[:,1],r[:,1],i[:,1],z[:,1],y[:,1]
  tmax = max(t1)
  
  conditionu = ut <= days
  conditiong = gt <= days
  conditionr = rt <= days
  conditioni = it <= days
  conditionz = zt <= days
  conditiony = yt <= days

  ut,gt,rt,it,zt,yt = np.extract(conditionu,ut),np.extract(conditiong,gt),np.extract(conditionr,rt),np.extract(conditioni,it),np.extract(conditionz,zt),np.extract(conditiony,yt)
  um,gm,rm,im,zm,ym = um[:len(ut)],gm[:len(gt)],rm[:len(rt)],im[:len(it)],zm[:len(zt)],ym[:len(yt)]
	
  if len(ut) == 0:
    ut = np.empty(1)
    ut[:] = np.NAN
    um = np.empty(1)
    um[:] = np.NAN
  if len(gt) == 0:
    gt = np.empty(1)
    gt[:] = np.NAN
    gm = np.empty(1)
    gm[:] = np.NAN
  if len(rt) == 0:
    rt = np.empty(1)
    rt[:] = np.NAN
    rm = np.empty(1)
    rm[:] = np.NAN
  if len(it) == 0:
    it = np.empty(1)
    it[:] = np.NAN
    im = np.empty(1)
    im[:] = np.NAN
  if len(zt) == 0:
    zt = np.empty(1)
    zt[:] = np.NAN
    zm = np.empty(1)
    zm[:] = np.NAN
  if len(yt) == 0:
    yt = np.empty(1)
    yt[:] = np.NAN
    ym = np.empty(1)
    ym[:] = np.NAN
	



#velocity model

#velocity correction
  if other['vcmb'] == False: vcor = 0
  if other['vcmb'] != False: vcor = veco(lsst['ra'],lsst['dec'])
  print 'velocity correction: ',vcor[0],' Km/s'
  print 'Einstein radius = ',RE,' meters'
  lar = days #tmax - tmin
  print lar,"LARGO"
  time=lar*86400
  lt=2.5*(10**13)


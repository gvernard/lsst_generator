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

    
  return mags#'Fin'


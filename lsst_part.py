import json
import sys
import numpy
import random
import os

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


# This code gets the ra and dec as an input, together with some LSST parameters (years, opsims file).
# The output is 6 files, one for each LSST filter, having the Modified Julian Date and depth in each column.
# These files can be empty if no LSST observation is planned on the given time scale.



if len(sys.argv) > 1:
    inp_file = sys.argv[1]
else:
    inp_file = "input.json"



print(">>>>>>>>>>>>>>>>>> Reading input <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
with open(inp_file) as json_input:    
    myinput = json.load(json_input)
    
outdir = myinput['path_to_out'] + str(myinput['outdir'])
path_to_dbfile = myinput['lsst']['path_to_dbfile']
runName = myinput['lsst']['survey']
ra = myinput['system']['ra']
dec = myinput['system']['dec']
days = myinput['lsst']['years']*365.25 # in days



print(">>>>>>>>>>>>>>>>>> Calculating angular diameter distances <<<<<<<<<<<<")
H0 = 67.7 # in km s^-1 Mpc^-1
Om0 = 0.309 # Omega matter at t=t_0 (now)
cosmo = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om0)
Ds  = cosmo.angular_diameter_distance(myinput['system']['zs']).value # in Mpc
Dl  = cosmo.angular_diameter_distance(myinput['system']['zl']).value # in Mpc
Dls = cosmo.angular_diameter_distance_z1z2(myinput['system']['zl'],myinput['system']['zs']).value # in Mpc



print(">>>>>>>>>>>>>>>>>> Reading LSST dates and depths <<<<<<<<<<<<<<<<<<<<<")
outDir = outdir + '/output'
dbFile = path_to_dbfile + runName + '.db'
opsimdb = db.opsimDatabase.OpsimDatabase(dbFile)
resultsDb = db.ResultsDb(outDir=outDir)

#ra=ra*15.0
# SNR limit (Don't use points below this limit)
snrLimit = 5.0
# The pass metric just passes data straight through.
metric = metrics.PassMetric(cols=['filter','fiveSigmaDepth','observationStartMJD'])
slicer = slicers.UserPointsSlicer(ra,dec,lonCol='fieldRA',latCol='fieldDec')
sql = ''
bundle = metricBundles.MetricBundle(metric,slicer,sql)
bg =  metricBundles.MetricBundleGroup({0:bundle},opsimdb,outDir=outDir,resultsDb=resultsDb)
bg.runAll()
bundle.metricValues.data[0]['filter']

filters = myinput['filters']
print('%i Observations total at this point (All SNR levels)' % bundle.metricValues.data[0].size)
for fname in filters:
    good = numpy.where(bundle.metricValues.data[0]['filter'] == fname)
    print('%i Observations in %s' % (good[0].size, fname))
    mags=bundle.metricValues.data[0]['observationStartMJD'][good[0]],bundle.metricValues.data[0]['fiveSigmaDepth'][good[0]]
    mags=numpy.array(mags)
    mags=numpy.transpose(mags)
    mags=mags[mags[:,0].argsort()]
    numpy.savetxt(outDir+'/'+fname+"_dates.dat",numpy.c_[mags[:,0],mags[:,1]],header="JD   5-sigma-depth")


# These two files are generated somewhere above but it is hard to find where, they are not used in any way though.
os.remove(outDir+"/resultsDb_sqlite.db")
os.remove(outDir+"/opsim_Pass_filter_fiveSigmaDepth_observationStartMJD_USER.npz")


print(">>>>>>>>>>>>>>>>>> Write json input for the GERLUMPH part <<<<<<<<<<<<")
out = myinput

out["system"]["Ds"]  = Ds
out["system"]["Dl"]  = Dl
out["system"]["Dls"] = Dls
out["path_2_dates"]  = myinput["path_to_out"] + str(myinput["outdir"]) + "/output/"
out["path_2_custom"] = myinput["path_to_out"] + str(myinput["outdir"]) + "/custom/"
out["path_2_output"] = myinput["path_to_out"] + str(myinput["outdir"]) + "/output/"

with open(myinput["path_to_out"] + str(myinput["outdir"]) + '/input_cpp.json','w') as outfile:
    json.dump(out,outfile,indent=4,separators=(',',':'))

print('Done.')


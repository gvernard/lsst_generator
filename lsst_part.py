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
# The output is 6 files, one for each LSST filter, having the Julian date and depth in each column.
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
ra = myinput['lsst']['ra']
dec = myinput['lsst']['dec']
days = myinput['lsst']['years']*365.25 # in days



print(">>>>>>>>>>>>>>>>>> Calculating angular diameter distances <<<<<<<<<<<<")
H0 = 67.7 # in km s^-1 Mpc^-1
Om0 = 0.309 # Omega matter at t=t_0 (now)
cosmo = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om0)
Ds  = cosmo.angular_diameter_distance(myinput['other']['zs']).value # in Mpc
Dl  = cosmo.angular_diameter_distance(myinput['other']['zl']).value # in Mpc
Dls = cosmo.angular_diameter_distance_z1z2(myinput['other']['zl'],myinput['other']['zs']).value # in Mpc



print(">>>>>>>>>>>>>>>>>> Reading LSST dates and depths <<<<<<<<<<<<<<<<<<<<<")
outDir = outdir + '/output'
dbFile = path_to_dbfile + runName + '.db'
opsimdb = db.opsimDatabase.OpsimDatabase(dbFile)
resultsDb = db.ResultsDb(outDir=outDir)

#ra=ra*15.0
# SNR limit (Don't use points below this limit)
snrLimit = 5.
# The pass metric just passes data straight through.
metric = metrics.PassMetric(cols=['filter','fiveSigmaDepth','observationStartMJD'])
slicer = slicers.UserPointsSlicer(ra,dec,lonCol='fieldRA',latCol='fieldDec')
sql = ''
bundle = metricBundles.MetricBundle(metric,slicer,sql)
bg =  metricBundles.MetricBundleGroup({0:bundle},opsimdb,outDir=outDir,resultsDb=resultsDb)
bg.runAll()
bundle.metricValues.data[0]['filter']

filters = ['u','g','r','i','z','y']
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
out = {}
out["profile"] = myinput["profile"]

out["lrest"]   = myinput["lrest"]

out["output"]  = {}
out["output"]["Nlc"]  = myinput["output"]["numlc"]
out["output"]["seed"] = random.randint(1,999)
out["output"]["full_data"] = myinput["output"]["full_data"]
out["output"]["degraded_data"] = myinput["output"]["degraded_data"]
out["output"]["velocities"] = myinput["output"]["velocities"]

out["maps"] = myinput["maps"]

out["vel"] = {}
out["vel"]["years"] = myinput["lsst"]["years"]
out["vel"]["ra"]    = myinput["lsst"]["ra"]
out["vel"]["dec"]   = myinput["lsst"]["dec"]
out["vel"]["sigma_l"]    = myinput["other"]["sgal"]
out["vel"]["sigma_s"]    = myinput["other"]["squas"]
out["vel"]["sigma_disp"] = myinput["other"]["sdisp"]
out["vel"]["epsi"]       = myinput["other"]["epsi"]
out["vel"]["zs"]  = myinput["other"]["zs"]
out["vel"]["zl"]  = myinput["other"]["zl"]
out["vel"]["Ds"]  = Ds
out["vel"]["Dl"]  = Dl
out["vel"]["Dls"] = Dls

out["filters"] = filters

out["path_2_dates"]  = myinput["path_to_out"] + str(myinput["outdir"]) + "/output/"
out["path_2_custom"] = myinput["path_to_out"] + str(myinput["outdir"]) + "/custom/"
out["path_2_output"] = myinput["path_to_out"] + str(myinput["outdir"]) + "/output/"

errbase = [
    myinput["lsst"]["errbaseu"],
    myinput["lsst"]["errbaseg"],
    myinput["lsst"]["errbaser"],
    myinput["lsst"]["errbasei"],
    myinput["lsst"]["errbasez"],
    myinput["lsst"]["errbasey"]
]
out["errbase"] = errbase

with open(myinput["path_to_out"] + str(myinput["outdir"]) + '/input_cpp.json','w') as outfile:
    json.dump(out,outfile)

print('Done.')


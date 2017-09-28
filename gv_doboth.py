# -*- coding: utf-8 -*-
#from microtimo33 import *
from getdate import *
from gv_microtimo_both import *
import numpy as np
import os
import time as pytime
import json
import sys

def gv_doboth(inp_file):
	startall=pytime.clock()
        

        ##############################################################################################################################
        # Get input
        ##############################################################################################################################
        print ">>>>>>>>>>>>>>>>>> Reading input <<<<<<<<<<<<<<<<<<<<<<<<<<<"
        with open(inp_file) as json_input:    
                myinput = json.load(json_input)

        outdir = myinput['path_to_out'] + str(myinput['outdir'])

        ##############################################################################################################################
        # Setup LSST things
        ##############################################################################################################################
        print ">>>>>>>>>>>>>>>>>> Doing LSST things <<<<<<<<<<<<<<<<<<<<<<<"
	a=getdate(outdir,myinput['lsst']['path_to_dbfile'],myinput['lsst']['survey'],myinput['lsst']['ra'],myinput['lsst']['dec'])
	endlsst=pytime.clock()
	


        ##############################################################################################################################
        # Calculate light curves
        ##############################################################################################################################
        print ">>>>>>>>>>>>>>>>>> Calculating light curves <<<<<<<<<<<<<<<<"
	startcurves=pytime.clock()
#	micro(rc, values[0], values[1], values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9],values[10],values[11],values[12],values[13],values[14],values[15],values[16],values[17],values[18],values[19],values[20],values[21],values[22],values[23],values[24], get_plot=trufal[values[25]], write_table=trufal[values[26]]  , write_img=trufal[values[27]]  ,mass=trufal[values[28]], vcmb=trufal[values[29]],nomacro=trufal[values[30]])



        micro(outdir,myinput['source'],myinput['magmap'],myinput['lsst'],myinput['other'],myinput['output'])
	
	endall=pytime.clock()
	print 'Time for dates',endlsst-startall
	print 'Time for curves',endall-startcurves

        return
        


#print str(sys.argv)
if len(sys.argv) > 1:
        inp_file = sys.argv[1]
else:
        inp_file = "input.json"

gv_doboth(inp_file)
print 'Done.'        
        
        
        

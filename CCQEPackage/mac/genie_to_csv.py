'''
    usage:
    ------
    Afro files:
    python mac/genie_to_csv.py -mA=0.99 -evf=0.01
    python mac/genie_to_csv.py --DataType=argon40_numu_mA_0_99 -evf=0.01
    make && python mac/genie_to_csv.py --DataType=uboone_CCinclMEC_argon40_numu_mA_0_99 -evf=1
    
    My files:
    python mac/genie_to_csv.py --DataType=40Ar_spline_CCinclMEC_muons_mA0.40 -evf=0.01
    
'''

import ROOT , time , os, sys , math
sys.path.insert(0, '/Users/erezcohen/larlite/UserDev/mySoftware/MySoftwarePackage/mac')
import matplotlib as mpl, pandas as pd , numpy as np
import input_flags
from ROOT import GenieFile
from matplotlib import pyplot as plt
from my_tools import *
from calc_tools import *
flags = input_flags.get_args()
print flags


#infilename = flags.DataType

gf = GenieFile( "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/mA/afro_genie_samples/" # path
               ,"CC_100k_mA_%.2f"%float(flags.mA) # RootFileName
               ,"gst" # RootFileName
               ,flags.verbose )
gf.HeaderCSV()


Nevents = gf.GetNevents()
ctr_CC1p0pi = 0
for i_event in range(int(flags.evnts_frac*Nevents)): #{
    if i_event%(Nevents/10)==0: print 'reading event',i_event,'(%.0f'%(100*float(i_event)/Nevents),'%)'
    gf.ReadEvent(i_event)
    gf.SetTopology()
    #    gf.Print()
    if gf.GetCC_1p_200MeVc_0pi() is True:#{
        #        print 'CC1p0pi event!'
        ctr_CC1p0pi += 1
    #}
    gf.StreamToCSV()
#}

gf.EndJob()

print 'done genie_to_csv. wrote %d CC1p0pi events'%ctr_CC1p0pi

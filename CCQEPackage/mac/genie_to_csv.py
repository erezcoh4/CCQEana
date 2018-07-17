'''
    usage:
    ------
    Afro files, June 2018:
    make && python mac/genie_to_csv.py -mA=0.99 -evf=0.001
    make && python mac/genie_to_csv.py -mA=0.80 -evf=1 --nskip=80000 -v6
    
    
    and to run on all files:
    python mac/genie_to_csv.py -mA=0 -evf=1
    
    My files, May 2018:
    python mac/genie_to_csv.py --DataType=40Ar_spline_CCinclMEC_muons_mA0.40 -evf=0.01
    
'''

import ROOT , time , os, sys , math
sys.path.insert(0, '/Users/erezcohen/larlite/UserDev/mySoftware/MySoftwarePackage/mac')
import matplotlib as mpl, pandas as pd , numpy as np, input_flags
from ROOT import GenieFile
from matplotlib import pyplot as plt
from my_tools import *
from calc_tools import *
flags = input_flags.get_args()
print flags


#infilename = flags.DataType
acceptance_map_path = "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/acceptance_maps/"
Pmu_theta_acceptance_map_name = "Pmu_theta_12x12_bins"
Pp_theta_acceptance_map_name = "Pp_theta_12x12_bins"
Q2_acceptance_map_name = "Q2_11_bins"

# run on a single files
if flags.mA>0:#{
    print "flags.mA:",flags.mA
    gf = GenieFile( "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/mA/afro_genie_samples/" # path
                   ,"CC_100k_mA_%.2f"%float(flags.mA) # RootFileName
                   ,"gst" # RootFileName
                   ,flags.verbose
                   ,acceptance_map_path
                   ,Pmu_theta_acceptance_map_name
                   ,Pp_theta_acceptance_map_name
                   ,Q2_acceptance_map_name)
    gf.HeaderCSV()

    Nevents = gf.GetNevents()
    ctr_CC1p0pi = 0
    print "stepping through events ",flags.nskip,"to",int(flags.evnts_frac*Nevents)
    for i_event in range(flags.nskip,int(flags.evnts_frac*Nevents)): #{
        if i_event%(Nevents/10)==0: print 'reading event',i_event,'(%.0f'%(100*float(i_event)/Nevents),'%)'
        
        gf.ReadEvent(i_event)
        gf.SetTopology()
        gf.MimicDetectorVolume()
        gf.SetMicroBooNEWeights()
        
        if flags.verbose>2:#{
            print 'XXXXXXXX \n event %d \nXXXXXXXX'%i_event
            gf.Print()
        #}
        if gf.GetCC_1p_200MeVc_0pi() is True:#{
        #        print 'CC1p0pi event!'
            ctr_CC1p0pi += 1
        #}
        gf.StreamToCSV()
    #}
    gf.EndJob()
    print 'done stepping through %d CC1p0pi events with mA=%.2f'%(ctr_CC1p0pi,flags.mA)
#}


# run on all files
elif flags.mA==0:#{
    for mA in [0.6,0.7,0.8,0.9,0.95,0.99,1.05,1.1,1.2,1.3,1.4]:#{
        gf = GenieFile( "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/mA/afro_genie_samples/" # path
                       ,"CC_100k_mA_%.2f"%float(mA) # RootFileName
                       ,"gst" # RootFileName
                       ,flags.verbose
                       ,acceptance_map_path
                       ,Pmu_theta_acceptance_map_name
                       ,Pp_theta_acceptance_map_name
                       ,Q2_acceptance_map_name)
        gf.HeaderCSV()
        Nevents = gf.GetNevents()
        ctr_CC1p0pi = 0
        for i_event in range(flags.nskip,int(flags.evnts_frac*Nevents)): #{
            if i_event%(Nevents/10)==0: print 'reading event',i_event,'(%.0f'%(100*float(i_event)/Nevents),'%)'
            
            gf.ReadEvent(i_event)
            gf.SetTopology()
            gf.MimicDetectorVolume()
            gf.SetMicroBooNEWeights()
            
            if gf.GetCC_1p_200MeVc_0pi() is True: ctr_CC1p0pi += 1
            gf.StreamToCSV()
        #}
        gf.EndJob()
        print 'done stepping through %d CC1p0pi events with mA=%.2f'%(ctr_CC1p0pi,mA)
    #}
#}


print 'done genie_to_csv'

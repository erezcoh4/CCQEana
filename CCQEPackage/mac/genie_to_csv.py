'''
    usage:
    ------
    Afro files, June 2018:
    make && python mac/genie_to_csv.py -evf=0.001
    python mac/genie_to_csv.py --option=CC_4_9E20_POT_mA_0_99_nominal -evf=1
    python mac/genie_to_csv.py --option=CC_4_9E20_POT_mA_0_99_hA2015 -evf=1
    python mac/genie_to_csv.py --option=CC_4_9E20_POT_mA_0_99_hN2015 -evf=1
    python mac/genie_to_csv.py --option=CC_4_9E20_POT_mA_0_99_hA_SRC -evf=1
    python mac/genie_to_csv.py --option=CC_4_9E20_POT_mA_0_99_hA_Tune3 -evf=1
    
    make && python mac/genie_to_csv.py -evf=1 --nskip=80000 -v6
    
    
'''

import ROOT , time , os ,  sys
from ROOT import GenieFile
sys.path.insert(0, '/Users/erezcohen/larlite/UserDev/mySoftware/MySoftwarePackage/mac')
import input_flags
from my_tools import *
from calc_tools import *
flags = input_flags.get_args()
print flags


genie_path = "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/afro_genie_samples/" # path
filename = flags.option
if filename == '': filename = "CC_4_9E20_POT_mA_0_99"

acceptance_map_path = "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/acceptance_maps/"
Pmu_theta_acceptance_map_name = "Pmu_theta_12x12_bins"
Pp_theta_acceptance_map_name = "Pp_theta_12x12_bins"
Q2_acceptance_map_name = "Q2_11_bins"
Q2_gen_rec_map_name = "Q2_gen_rec"


# run on a single files
gf = GenieFile( genie_path
                   ,filename
                   ,"gst" # RootFileName
                   ,flags.verbose
                   ,acceptance_map_path
                   ,Pmu_theta_acceptance_map_name
                   ,Pp_theta_acceptance_map_name
                   ,Q2_acceptance_map_name)
gf.SetQ2_gen_rec_map( acceptance_map_path, Q2_gen_rec_map_name )
gf.HeaderCSV()

Nevents = gf.GetNevents()
ctr_CC1p = 0
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
    if gf.GetCC_1p_200MeVc() is True:#{
        ctr_CC1p += 1
    #}
    gf.StreamToCSV()
#}
gf.EndJob()
print 'done stepping through %d CC1p events with'%(ctr_CC1p)
print 'done genie_to_csv'

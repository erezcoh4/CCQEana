'''
    The sam list Joel supplied for files without PandoraCosmic removal pass 
    had the following weird structure:
    
    fnal-dcache:/pnfs/uboone/scratch/users/joelam/mcc8/v06_26_01_11/reco_noCosmic/prodcosmics_corsika_cmc_uboone/4182171_1241/prodcosmics_corsika_cmc_uboone_2119_20180214T123435_gen2_3163479b-edf4-4453-8e47-451ee5837091_20180226T221126_reco1_20180227T001006_reco2.root
    
    This code peals-off the "fnal-dcache:"
    '''

filename = "/Users/erezcohen/Desktop/prodcosmics_reco_withoutPandoraCosmicRemovalPass_joelam.list"
string = ''
with open(filename,'r') as file:
    lines = file.readlines()
    for line in lines:
        list=line.split("fnal-dcache:")
        fixed_filename = list[1]
        # print 'adding',fixed_filename
        string = string+fixed_filename
print '---------'
print string
'''
    usage:
    -----
    * for MC-BNB/DATA-cosmic:
    
    python scripts/download_files.py --name=prodgenie_bnb_nu_uboone_overlay_mcc8_v9
    python scripts/download_files.py --name=prodgenie_bnb_nu_uboone_overlay_mcc8_v9 --option=makeup --continue_makeup=7216920_63 --ctr=761


    python scripts/download_files.py --name=prodgenie_bnb_nu_uboone_overlay_mcc8_v4
    
    python scripts/download_files.py --name=ecohen_physical_files_adi_prodgenie_bnb_nu_uboone_overlay_cosmic_data_100K_reco2
    python scripts/download_files.py --name=ecohen_physical_files_adi_prodgenie_bnb_nu_uboone_overlay_cosmic_data_100K_reco2 --option=makeup --continue_makeup=3547582_592 --ctr=765

    * for MC-BNB/MC-cosmic:

    python scripts/download_files.py --name=prodgenie_bnb_nu_cosmic_uboone_mcc8.7_reco2_dev


    * for On Beam:
    
    python scripts/download_files.py --name=prod_reco_optfilter_bnb_v12_unblind_mcc8_05
    python scripts/download_files.py --name=prod_reco_optfilter_bnb_v12_unblind_mcc8_04
    
    * for Off Beam:
    
    python scripts/download_files.py --name=prod_reco_optfilter_extbnb_v12_mcc8_dev_05
    python scripts/download_files.py --name=prod_reco_optfilter_extbnb_v12_mcc8_dev_04
    
    * for cosmic-pairs analysis:
    
    python scripts/download_files.py --name=prodcosmics_corsika_cmc_uboone_mcc8.7_reco2 --tag=withPandoraCosmic_pass
    
'''

import sys, os, time, argparse
from prompter import yesno
time_name = "%4d_%02d_%02d" % time.localtime()[0:3]


parser = argparse.ArgumentParser()
parser.add_argument('-n','--name', default='ccqe_ana_MCBNBCosmicDATA', type=str )
parser.add_argument('-o','--option', default='download', type=str )
parser.add_argument('--continue_makeup', default='3215875_563', type=str )
parser.add_argument('--ctr', default=0, type=int )
parser.add_argument('--tag', default='withPandoraCosmic_pass', type=str )
name = parser.parse_args().name
option = parser.parse_args().option
continue_makeup = parser.parse_args().continue_makeup
ctr = parser.parse_args().ctr
tag = parser.parse_args().tag

uboone = 'ecohen@uboonegpvm04.fnal.gov'
csv_path = '/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/'

# for the new overlay: AdiReco27List_18Jan2018
if name=="ccqe_ana_MCBNBCosmicDATA":
    default_pnfsjob = "/pnfs/uboone/scratch/users/ecohen/mcc8/06_26_01_09/ccqe_ana/"+name
    default_outdirname = csv_path+"ccqe_candidates/"

# new overlay using SAM definition / MCC8.7 MC-BNB + MC-Cosmic
if "prodgenie_bnb_nu_uboone_overlay_cosmic_data" in name or "prod_reco" or "prodgenie" in name:
    default_book = "/uboone/data/users/ecohen/book/mcc8/06_26_01_09/ccqe_ana/"+name+"_CCQE"
    default_pnfsjob = "/pnfs/uboone/scratch/users/ecohen/mcc8/06_26_01_09/ccqe_ana/"+name+"_CCQE"
    default_outdirname = csv_path+"ccqe_candidates/"


# pandora pairs
elif "prodcosmics_corsika" in name:
    default_pnfsjob = "/pnfs/uboone/scratch/users/ecohen/mcc8/06_26_01_09/cosmic_pairs/" + name + "_" + tag + "/"
    default_outdirname = csv_path+"pandora_pairs/"

default_indirname = csv_path+"from_grid/"+name+"/"


book = raw_input("enter book:...<"+default_book+">") or default_book
print 'book: ',book
#pnfsjob = raw_input("enter pnfs:...<"+default_pnfsjob+">") or default_pnfsjob
#print 'pnfsjob: ',pnfsjob
indirname = raw_input("enter indir:...<"+default_indirname+">") or default_indirname
print 'indirname: ',indirname
outdirname = raw_input("enter outdir:...<"+default_outdirname+">") or default_outdirname
print 'outdirname: ',outdirname
print


# for the new overlay: AdiReco27List_18Jan2018
if name=="ccqe_ana_MCBNBCosmicDATA":
    outfilename = outdirname+'/'+name+'_'+time_name+'_vertices.csv'
    summaryfilename = csv_path+'summary/'+name+'_'+time_name+'_summary.csv'

# new overlay using SAM definition
if "prodgenie_bnb_nu_uboone_overlay_cosmic_data" in name or "prod_reco" in name or "prodgenie" in name:
    outfilename = outdirname+'/'+name+'_'+time_name+'_vertices.csv'
    summaryfilename = csv_path+'summary/'+name+'_'+time_name+'_summary.csv'


# pandora pairs
elif "prodcosmics_corsika" in name:
    outfilename = outdirname+'/'+name+ "_" + tag +'_'+time_name+'_vertices.csv'
    summaryfilename = csv_path+'summary/'+name+ "_" + tag + '_'+time_name+'_summary.csv'


# step 1: create a list of files to download
if yesno('# step 1: create a list of files to download?'):#{
    print 'created the <indirname> directory (if exists already, remove the existing one):'
    print indirname
    print
    os.system("rm -fr "+indirname)
    os.system("mkdir "+indirname)
    print 'creating a list of files to download from',book
    print 'into '+indirname+'/files_to_download.list'
    print "ssh "+uboone+" ls "+book+"/*/*.csv > "+indirname+"/files_to_download.list"
    print 'loading...'
    os.system("ssh "+uboone+" ls "+book+"/*/*.csv > "+indirname+"/files_to_download.list")
#}
print
# step 2: grab this list
if yesno('# step 2: grab this list?'):#{
    files_to_download_name = "files_to_download"
    if option=="makeup": files_to_download_name = files_to_download_name + "_since_" + continue_makeup
    with open(indirname+'/'+files_to_download_name+".list") as f:
        files = f.read().splitlines()
#}
print

# step 3: iterate over the list and download the files
if yesno('# step 3: iterate over the list and download the files?'):#{
    for file in files:#{
        try:
            print 'downloading',file
            destination = file.rsplit('/', 1)[1] + '_%d.csv'%ctr
            print 'into destination',destination
            os.system('scp '+uboone+':'+file+' '+indirname+'/'+destination)
            ctr = ctr + 1
        except IOError:
            pass
    #}
#}
print

# step 4: merge the csv files into a large dataframe
print '# step 4: merge the csv files into a large dataframe'
print 'creating ',name+'_vertices.csv file'
import os, pandas as pd
list_files = os.listdir(indirname)
print "list of files to combine from:",list_files
print
summary_df_array,vertices_df_array,tracks_df_array ,genie_df_array,events_df_array = [],[],[],[],[]
list_files = os.listdir(indirname)
for file in list_files:#{
    if os.stat(indirname+'/'+file).st_size == 0: continue
    try:
        if "vertices.csv" in file:#{
            print 'adding ',indirname+'/'+file
            vertices_df_array.append(pd.read_csv( ( indirname+'/'+file )))
        #}
        if "summary.csv" in file:#{
            print 'adding ',indirname+'/'+file
            summary_df_array.append(pd.read_csv( ( indirname+'/'+file )))
        #}
        if "tracks.csv" in file:#{
            print 'adding ',indirname+'/'+file
            tracks_df_array.append(pd.read_csv( ( indirname+'/'+file )))
        #}
        if "genie.csv" in file:#{
            print 'adding ',indirname+'/'+file
            genie_df_array.append(pd.read_csv( ( indirname+'/'+file )))
        #}
        if "events.csv" in file:#{
            print 'adding ',indirname+'/'+file
            events_df_array.append(pd.read_csv( ( indirname+'/'+file )))
        #}
    except IOError:
        pass
#}
vertices_all = pd.concat(vertices_df_array)
print len(vertices_all),'total pair-vertices'
vertices_all.to_csv(outfilename)
print 'concatenated %d'%len(vertices_all),'vertices to\n',outfilename,'\ndone.'
print
summary_all = pd.concat(summary_df_array)
print len(summary_all),'is the length summary of all events'
summary_all.to_csv(summaryfilename)
print 'concatenated %d'%len(summary_all),'summary files to\n',summaryfilename,'\ndone.'
print
if (len(tracks_df_array)>0):#{
    tracks_all = pd.concat(tracks_df_array)
    print len(tracks_all),'total tracks'
    outfilename = csv_path+'tracks/'+name+'_'+time_name+'_tracks.csv'
    tracks_all.to_csv(outfilename)
    print 'concatenated %d'%len(tracks_all),'tracks to\n',outfilename,'\ndone.'
#}
print
if (len(genie_df_array)>0):#{
    genie_all = pd.concat(genie_df_array)
    print len(genie_all),'total genie'
    outfilename = csv_path+'genie/'+name+'_'+time_name+'_genie.csv'
    genie_all.to_csv(outfilename)
    print 'concatenated %d'%len(genie_all),'genie to\n',outfilename,'\ndone.'
#}
print
if (len(events_df_array)>0):#{
    events_all = pd.concat(events_df_array)
    print len(events_all),'events in total'
    outfilename = csv_path+'events/'+name+'_'+time_name+'_events.csv'
    events_all.to_csv(outfilename)
    print 'concatenated %d'%len(events_all),'events to\n',outfilename,'\ndone.'
#}
print

os.system('say "download files from grid has completed"')
os.system('say "concatenated %d vertices to outputfile from %d events".'%(len(vertices_all),len(events_all)))


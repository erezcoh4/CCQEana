'''
    python scripts/download_csv_from_grid.csv
    '''
# combine *vertices.csv and *summary.csv files from a grid job
import os, shutil

# for the previous overlay: prod_reco2_extbnb_v8_mcc8
# default_name = "prod_reco2_extbnb_v8_mcc8"
# default_indirname = "/uboone/data/users/ecohen/book/mcc8/v06_42_00/ccqe_ana/"+name+"_CCQE/"
# default_outdirname = "/uboone/data/users/ecohen/CCQEanalysis/csvFiles/ccqe_candidates/"+name+"/"


# for the new overlay: AdiReco27List_18Jan2018
default_name = "ccqe_ana_MCBNBCosmicDATA"
name = raw_input("enter list name:...<"+default_name+">") or default_name

default_indirname = "/pnfs/uboone/scratch/users/ecohen/mcc8/06_26_01_09/ccqe_ana/"+name+"/"
indirname = raw_input("enter indir:...<"+default_indirname+">") or default_indirname

default_outdirname = "/uboone/data/users/ecohen/CCQEanalysis/csvFiles/ccqe_candidates/"+name+"/"
outdirname = raw_input("enter outdir:...<"+default_outdirname+">") or default_outdirname

if not os.path.exists(outdirname):
    os.makedirs(outdirname)
    print 'created a directory:\n',outdirname

print 'copying ',name,'files to',outdirname
list_dir = os.listdir(indirname)
print "list of directories to copy from:",list_dir

ctr = 0
for dir_name in list_dir:#{
    if os.path.isdir(indirname+'/'+dir_name): #{
        print 'in directory',dir_name
        list_files = os.listdir(indirname+'/'+dir_name)
        print 'list_files:\n',list_files
        for file in list_files:#{
            if "vertices.csv" in file or "summary.csv" in file:#{
                source_to_cp = indirname+'/'+dir_name+'/'+file
                destination = outdirname+file+'_%d'%ctr+'.csv'
                print 'copying ',source_to_cp,'to',destination
                shutil.copyfile(source_to_cp, destination)
                ctr = ctr + 1
            #}
        #}
    #}
#}
print 'done, moved ',ctr,'files. see result in'
print outdirname

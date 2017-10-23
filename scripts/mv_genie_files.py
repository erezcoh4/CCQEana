# combine *vertices.csv files from a grid job
default_name = "prod_reco2_extbnb_v8_mcc8"
name = raw_input("enter list name:...<"+default_name+">") or default_name

defauld_indirname = "/uboone/data/users/ecohen/book/mcc8/v06_42_00/cc_genie/"+name+"_CCQE/"
indirname = raw_input("enter indir:...<"+defauld_indirname+">") or defauld_indirname

defauls_outdirname = "/uboone/data/users/ecohen/CCQEanalysis/csvFiles/genie/"+name+"/"
outdirname = raw_input("enter outdir:...<"+defauls_outdirname+">") or defauls_outdirname

print 'copying ',name,'files to',outdirname
import os, shutil
list_dir = os.listdir(indirname)
print "list of directories to copy from:",list_dir

ctr = 0
for dir_name in list_dir:#{
    if os.path.isdir(indirname+'/'+dir_name): #{
        print 'in directory',dir_name
        list_files = os.listdir(indirname+'/'+dir_name)
        for file in list_files:#{
            if "genie.csv" in file:#{
                source_to_cp = indirname+'/'+dir_name+'/'+file
                destination = outdirname+file+'_%d'%ctr+'.csv'
                print 'copying ',source_to_cp,'to',destination
                shutil.copyfile(source_to_cp, destination)
                ctr = ctr + 1
            #}
        #}
    #}
#}
print 'done, see:\n',outdirname

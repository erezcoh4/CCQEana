# combine *vertices.csv and *summary.csv files from a grid job


# for the previous overlay: prod_reco2_extbnb_v8_mcc8
# default_name = "prod_reco2_extbnb_v8_mcc8"
# default_indirname = "/uboone/data/users/ecohen/book/mcc8/v06_42_00/ccqe_ana/"+name+"_CCQE/"
# default_outdirname = "/uboone/data/users/ecohen/CCQEanalysis/csvFiles/ccqe_candidates/"+name+"/"
import sys , os , time
time_name = "%4d_%02d_%02d" % time.localtime()[0:3]


# for the new overlay: AdiReco27List_18Jan2018
default_name = "ccqe_ana_MCBNBCosmicDATA"
name = raw_input("enter list name:...<"+default_name+">") or default_name

default_indirname = "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/from_grid/"+name+"/"
indirname = raw_input("enter indir:...<"+default_indirname+">") or default_indirname

default_outdirname = "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/ccqe_candidates/"
outdirname = raw_input("enter outdir:...<"+default_outdirname+">") or default_outdirname


print 'creating ',name+'_vertices.csv file'
import os, pandas as pd
list_files = os.listdir(indirname)
print "list of files to combine from:",list_files

vertices_df_array = []
summary_df_array = []
list_files = os.listdir(indirname)
for file in list_files:#{
    if "vertices.csv" in file:#{
        print 'adding ',indirname+'/'+file
        vertices_df_array.append(pd.read_csv( ( indirname+'/'+file )))
    #}
    if "summary.csv" in file:#{
        print 'adding ',indirname+'/'+file
        summary_df_array.append(pd.read_csv( ( indirname+'/'+file )))
    #}
#}
vertices_all = pd.concat(vertices_df_array)
print len(vertices_all),'total pair-vertices'
outfilename = outdirname+'/'+name+'_'+time_name+'_vertices.csv'
vertices_all.to_csv(outfilename)
print 'concatenated %d'%len(vertices_all),'vertices to\n',outfilename,'\ndone.'

summary_all = pd.concat(summary_df_array)
print len(summary_all),'is the length summary of all events'
outfilename = outdirname+'/'+name+'_'+time_name+'_summary.csv'
summary_all.to_csv(outfilename)
print 'concatenated %d'%len(summary_all),'summary files to\n',outfilename,'\ndone.'

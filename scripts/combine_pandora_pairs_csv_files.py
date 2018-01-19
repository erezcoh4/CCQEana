# combine *vertices.csv files from a grid job
default_name = "prodcosmics_corsika_cmc_uboone_mcc8.4_reco2"
name = raw_input("enter list name:...<"+default_name+">") or default_name

defauld_indirname = "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/from_grid/"+name+"/"
indirname = raw_input("enter indir:...<"+defauld_indirname+">") or defauld_indirname

defauls_outdirname = "/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/pandora_pairs/"
outdirname = raw_input("enter outdir:...<"+defauls_outdirname+">") or defauls_outdirname


print 'creating ',name+'_vertices.csv file'
import os, pandas as pd
list_files = os.listdir(indirname)
print "list of files to combine from:",list_files

vertices_df_array = []
list_files = os.listdir(indirname)
for file in list_files:#{
    if "vertices.csv" in file:#{
        print 'adding ',indirname+'/'+file
        vertices_df_array.append(pd.read_csv( ( indirname+'/'+file )))
    #}
#}
vertices_all = pd.concat(vertices_df_array)
print len(vertices_all),'total pair-vertices'
outfilename = outdirname+'/'+name+'_vertices.csv'
vertices_all.to_csv(outfilename)
print 'concatenated %d'%len(vertices_all),'vertices to\n',outfilename,'\ndone.'

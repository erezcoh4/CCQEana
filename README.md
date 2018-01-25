# CCQEana in Yun-Tse branch of LArSoft 
rsync -avz LArSoftCodes/ ecohen@uboonegpvm04.fnal.gov:/uboone/app/users/ecohen/LArSoft/srcs/uboonecode/uboone/ErezCCQEana/

# CCQEana
rsync -avz LArSoftCodes/ ecohen@uboonegpvm04.fnal.gov:/uboone/app/users/ecohen/ErezCCQE/srcs/uboonecode/uboone/ErezCCQEana/ 

# to sync scripts
rsync -avz xml_files/  $uboone:/uboone/app/users/ecohen/ErezCCQE/xml_files/
rsync -avz fcl_files/  $uboone:/uboone/app/users/ecohen/ErezCCQE/fcl_files/
rsync -avz scripts/  $uboone:/uboone/app/users/ecohen/ErezCCQE/scripts/


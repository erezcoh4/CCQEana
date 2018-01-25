#
# Add any personal extra databases here:
#
#UPS_EXTRA_DIR=$HOME/p/upsdb; export UPS_EXTRA_DIR
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# get ups environment, and then setup the login product
#

pa=/grid/fermiapp/products/common/etc 
if [  -r "$pa/setups.sh" ]
then
  . "$pa/setups.sh"
  if ups exist login
  then
      setup login
  fi
fi

#
# make sure our .shrc gets run...
#
ENV=$HOME/.shrc
export ENV 
if [ "`basename $SHELL`" != ksh -a -r $ENV ]
then
    . $ENV
fi
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The default umask setting is 022, disabling writes except by
# owner.  Uncomment the following line to allow group writes.
# umask 002

# set prompt
PS1="<`hostname`> "; export PS1

# set default printer, etc.
# FLPQUE=fcc2w_ps;	export FLPQUE;
# FLPHOST=fnprt;		export FLPHOST;

# set timezone, esp if you don't want central time
# export TZ;		TZ=CST6CDT


# ---------------------------------------------------
# Erez
# Jan-25, 2018
# ---------------------------------------------------
echo 'now adding my own personal settings...'

# -------- exports
export data_uboone=/uboone/data/users/ecohen
export example_art_files=/uboone/data/users/ecohen/example_art_files
export app_uboone=/uboone/app/users/ecohen
export myJobs=$data_uboone/jobs
export fcl_files=/uboone/app/users/ecohen/ErezCCQE/fcl_files
export xml_files=/uboone/app/users/ecohen/ErezCCQE/xml_files
export ArtFilesMCC8=/uboone/data/users/ecohen/Lists/ArtFilesLists/MCC8
export CCQEcandidates=/uboone/data/users/ecohen/CCQEanalysis/csvFiles/ccqe_candidates/


# -------- aliases
alias ll='ls -lh --color=auto'
alias enw='emacs -nw'

# -------- GoTo
alias GoDataUboone='cd $data_uboone'
alias GoAppUboone='cd $app_uboone'
alias GoJobs='cd $myJobs'
alias GoErezCCQE='cd /uboone/app/users/$USER/ErezCCQE'
alias GoMakeErezCCQE='cd $app_uboone/ErezCCQE/build_slf6.x86_64/uboonecode/uboone/ErezCCQEana && mrbsetenv && make install -j4'


# -------- LArSoft shell
source /grid/fermiapp/products/uboone/setup_uboone.sh
setup uboonecode v06_26_01_10 -q e10:prof
GoErezCCQE
source localProducts_larsoft_v06_26_01_09_e10_prof/setup
mrbslp
echo "set my LArSoft preferences (v06_26_01_10 e10:prof)"



echo "setup Erez personal preferences"
echo "at $(pwd)"
echo "start working..."



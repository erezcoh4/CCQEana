import sys; sys.path.insert(0, '../..');
from ccqe_notebook_tools import *
from mupClassification_notebooks import *
from onbeam_offbeam_notebooks import *


'''
    Setup for CCQE-like cross-section extraction
    We extract CC1p cross section
    
    Signal:
    -------
    a CC event with a single proton above 200 MeV/c 
    and no charged pions above 70 MeV/c
    
    Background:
    -----------
    1mu-1p which are not CC1p
    other pairs
    cosmic in the overlay (which are correlated with the beam)
    cosmic not in the overlay (off-beam data)
'''
flux = 3.601e10 # cm^-2
flux_err = 0
Ntargets = 5.00e31
Ntargets_err = 0

Bins = dict({
            'Pmu':linspace(0.1,1.5,7)
            ,'theta(mu)':linspace(0,120,7)
            ,'cos(theta(mu))':linspace(-0.5,1,7)
            ,'Pp':linspace(0.2,1.2,7)
            ,'theta(p)':linspace(0,90,7)
            ,'cos(theta(p))':linspace(-0.2,1,7)
            })





vlabels = dict({
               'Pmu':"p_{\mu}"
               ,'theta(mu)':r"\theta_{\mu}"
               ,'cos(theta(mu))':r"\cos(\theta_{\mu})"
               ,'Pp':"p_{p}"
               ,'theta(p)':r"\theta_{p}"
               ,'cos(theta(p))':r"\cos(\theta_{p})"
               })


Vlabels = dict({
              'Pmu':r"$p_{\mu}$"
               ,'theta(mu)':r"$\theta_{\mu}$"
               ,'cos(theta(mu))':r"$\cos(\theta_{\mu})$"
               ,'Pp':r"$p_{p}$"
               ,'theta(p)':r"$\theta_{p}$"
               ,'cos(theta(p))':r"$\cos(\theta_{p})$"
              })

Units = dict({
             'Pmu':r"GeV/c"
             ,'theta(mu)':r"deg."
             ,'cos(theta(mu))':None
             ,'Pp':r"GeV/c"
             ,'theta(p)':r"deg."
             ,'cos(theta(p))':None
               })



Colors = dict({
              'overlay':'forestgreen'
              ,'CC 1p':'blue'
              ,'beam off':'black'
              ,'beam on':'tomato'
              })


Xsec_path = '/Users/erezcohen/Desktop/uBoone/CCQEanalysis/Xsec/'
Paths = dict({'migration maps':Xsec_path+'migration_maps/'
             ,'background maps':Xsec_path+'background_maps/'
             ,'efficiency maps':Xsec_path+'efficeincy_maps/'
              ,'1d Xsec':Xsec_path+'1d_Xsec/'})






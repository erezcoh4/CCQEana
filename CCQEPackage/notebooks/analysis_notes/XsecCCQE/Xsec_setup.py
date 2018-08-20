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


Bins = dict({
            'Pmu':linspace(0.1,1.5,7)
            ,'theta(mu)':linspace(0,120,7)
            ,'Pp':linspace(0.2,1.2,7)
            ,'theta(p)':linspace(0,90,7)
            })





Labels = dict({
            'Pmu':r"$p_{\mu}$ [GeV/c]"
            ,'theta(mu)':r"$\theta_{\mu}$ [deg.]"
            ,'Pp':r"$p_{p}$ [GeV/c]"
            ,'theta(p)':r"$\theta_{p}$ [deg.]"
            })




Colors = dict({
              'overlay':'forestgreen'
              ,'CC 1p':'blue'
              ,'beam off':'black'
              ,'beam on':'tomato'
              })
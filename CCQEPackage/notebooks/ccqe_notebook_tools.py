import sys, pandas as pd, numpy as np, matplotlib.pyplot as plt, ast, matplotlib as mpl
sys.path.insert(0, '/Users/erezcohen/larlite/UserDev/mySoftware/MySoftwarePackage/mac')
from plot_tools import *
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter, MultipleLocator, FormatStrFormatter
from matplotlib import ticker



#---------------------------------------------------------------------------------------------
def print_and_say(string=''):#{
    print string
    os.system('say "%s".'%string)
#}
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
# Aug-4, 2017
pair_types   = ['1mu-1p'   ,'other pairs','cosmic' ,'CC 1p 0pi'    ]
MClabels     = ['1$\\mu$1p','other pairs','cosmic' ,'CC 1p 0$\\pi$']
MCcolors     = ['teal'     ,'red'        ,'orange' ,'blue'         ]
MCcmaps      = ['Greens'   ,'Reds'       ,'Oranges','Blues'        ]
#---------------------------------------------------------------------------------------------


# ------------------------------------------------
# Aug-5, 2017
def set_color_label_cmap(t):
    pdg = t.GetMCpdgCode()
    color, label, cmap = 'purple','%d'%pdg,'Purples'
    if pdg==-13:
        color, label, cmap = 'black','$\\mu^+$','Greys'
    elif pdg==13:
        color, label, cmap = 'black','$\\mu^-$','Greys'
    elif pdg==-11:
        color, label, cmap = 'red','$e^+$','Reds'
    elif pdg==11:
        color, label, cmap = 'red','$e^-$','Reds'
    elif pdg==2212:
        color, label, cmap = 'blue','$p$','Blues'
    elif pdg==2112:
        color, label, cmap = 'salmon','$n$','Oranges'
    elif pdg==-211:
        color, label, cmap = 'teal','$\\pi^-$','Greens'
    elif pdg==211:
        color, label, cmap = 'teal','$\\pi^+$','Greens'
    elif pdg==111:
        color, label, cmap = 'teal','$\\pi^0$','Greens'
    elif pdg==22:
        color, label, cmap = 'cyan','$\\gamma$','cool'
    return color , label , cmap
# ------------------------------------------------

# ------------------------------------------------
def find_a_straight_line( x_array , y_array ):
    [x0,x1] = x_array
    [y0,y1] = y_array
    slope = (y1-y0)/((x1-x0) if np.abs(x1-x0)>0.001 else (y1-y0)/(0.001*np.sign(x1-x0)))
    intercept = y1 - slope*x1
    return slope,intercept
# ------------------------------------------------



# ------------------------------------------------
# April-23, 2018
def sample_in_FV(sample=None, max_FV_y = 110, # 115 in pandoraNu tracks collection
                 min_FV_z = 5, max_FV_z = 1037,
                 min_FV_x = 3, max_FV_x = 250): # 257
    sample_in_FV = sample[
                          (np.abs(sample['starty_muCandidate']) < max_FV_y)
                          & (np.abs(sample['starty_pCandidate']) < max_FV_y)
                          & (np.abs(sample['endy_muCandidate']) < max_FV_y)
                          & (np.abs(sample['endy_pCandidate']) < max_FV_y)
                          
                          & ((sample['startz_muCandidate'] > min_FV_z) & (sample['startz_muCandidate'] < max_FV_z) )
                          & ((sample['startz_pCandidate'] > min_FV_z) & (sample['startz_pCandidate'] < max_FV_z) )
                          & ((sample['endz_muCandidate'] > min_FV_z) & (sample['endz_muCandidate'] < max_FV_z) )
                          & ((sample['endz_pCandidate'] > min_FV_z) & (sample['endz_pCandidate'] < max_FV_z) )
                          
                          & ((sample['startx_muCandidate'] > min_FV_x) & (sample['startx_muCandidate'] < max_FV_x) )
                          & ((sample['startx_pCandidate'] > min_FV_x) & (sample['startx_pCandidate'] < max_FV_x) )
                          & ((sample['endx_muCandidate'] > min_FV_x) & (sample['endx_muCandidate'] < max_FV_x) )
                          & ((sample['endx_pCandidate'] > min_FV_x) & (sample['endx_pCandidate'] < max_FV_x) )
                          ]
    return sample_in_FV
# ------------------------------------------------



#
## ------------------------------------------------
#def sample_in_FV(sample=None, max_FV_y = 110, # 115 in pandoraNu tracks collection
#                 min_FV_z = 5, max_FV_z = 1037,
#                 min_FV_x = 3, max_FV_x = 250): # 257
#    sample_in_FV = sample[
#                          (np.abs(sample['starty_assigned_muon']) < max_FV_y)
#                          & (np.abs(sample['starty_assigned_proton']) < max_FV_y)
#                          & (np.abs(sample['endy_assigned_muon']) < max_FV_y)
#                          & (np.abs(sample['endy_assigned_proton']) < max_FV_y)
#
#                          & ((sample['startz_assigned_muon'] > min_FV_z) & (sample['startz_assigned_muon'] < max_FV_z) )
#                          & ((sample['startz_assigned_proton'] > min_FV_z) & (sample['startz_assigned_proton'] < max_FV_z) )
#                          & ((sample['endz_assigned_muon'] > min_FV_z) & (sample['endz_assigned_muon'] < max_FV_z) )
#                          & ((sample['endz_assigned_proton'] > min_FV_z) & (sample['endz_assigned_proton'] < max_FV_z) )
#
#                          & ((sample['startx_assigned_muon'] > min_FV_x) & (sample['startx_assigned_muon'] < max_FV_x) )
#                          & ((sample['startx_assigned_proton'] > min_FV_x) & (sample['startx_assigned_proton'] < max_FV_x) )
#                          & ((sample['endx_assigned_muon'] > min_FV_x) & (sample['endx_assigned_muon'] < max_FV_x) )
#                          & ((sample['endx_assigned_proton'] > min_FV_x) & (sample['endx_assigned_proton'] < max_FV_x) )
#                          ]
#    return sample_in_FV
## ------------------------------------------------



#---------------------------------------------------------------------------------------------
# April-30, 2017
def get_fraction_in_cut( data=None , cut_var='distance', mul=1 , xmin=0.1, xmax=10 , Nbins=10 ,  cut_type= 'max' ):
    x_array = np.linspace(xmin,xmax,Nbins)
    frac , frac_err = [] , []
    denominator = len(data)
    
    for x in x_array:
        if cut_type is 'max':
            reduced = data[mul*data[cut_var]<x]
        elif cut_type is 'min':
            reduced = data[mul*data[cut_var]>x]
        numerator = float(len(reduced))
        
        frac.append(100 * numerator / denominator)
        frac_err.append( frac[-1] * np.sqrt(1./numerator + 1./denominator) ) if numerator>0 else frac_err.append( frac[-1]/np.sqrt(denominator) )
    
    return np.array(x_array), np.array(frac) , np.array(frac_err)
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
# April-30, 2017
def get_fraction_in_symmetriccut( data=None , cut_var='delta_phi', mul=1,xcenter=0.1, delta_x_min=0, delta_x_max=100 , Nbins=10 ):
    delta_x_array = np.linspace(delta_x_min,delta_x_max,Nbins)
    frac , frac_err = [] , []
    denominator = len(data)
    
    for delta_x in delta_x_array:
        reduced = data[np.abs(mul*data[cut_var]-xcenter)<delta_x]
        numerator = float(len(reduced))
        frac.append(100 * numerator / denominator)
        frac_err.append( frac[-1] * np.sqrt(1./numerator + 1./denominator) ) if numerator>0 else frac_err.append( frac[-1]/np.sqrt(denominator) )
    
    return np.array(delta_x_array), np.array(frac) , np.array(frac_err)
#---------------------------------------------------------------------------------------------

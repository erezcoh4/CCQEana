import sys, pandas as pd, numpy as np, matplotlib.pyplot as plt, ast, matplotlib as mpl
sys.path.insert(0, '/Users/erezcohen/larlite/UserDev/mySoftware/MySoftwarePackage/mac')
from plot_tools import *
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter, MultipleLocator, FormatStrFormatter
from matplotlib import ticker




#---------------------------------------------------------------------------------------------
# Aug-4, 2017
pair_types = ['1mu-1p'   ,'other pairs','cosmic' ,'CC 1p 0pi'    ]
labels     = ['1$\\mu$1p','other pairs','cosmic' ,'CC 1p 0$\\pi$']
colors     = ['teal'     ,'red'        ,'orange' ,'blue'         ]
cmaps      = ['Greens'   ,'Reds'       ,'Oranges','Blues'        ]
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
        color, label, cmap = 'cyal','$\\gamma$','cool'
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
# Aug-4, 2017
def sample_in_FV(sample=None, max_FV_y = 110, # 115 in pandoraNu tracks collection
                 min_FV_z = 5, max_FV_z = 1037,
                 min_FV_x = 3, max_FV_x = 250): # 257
    sample_in_FV = sample[
                          (np.abs(sample['starty_assigned_muon']) < max_FV_y)
                          & (np.abs(sample['starty_assigned_proton']) < max_FV_y)
                          & (np.abs(sample['endy_assigned_muon']) < max_FV_y)
                          & (np.abs(sample['endy_assigned_proton']) < max_FV_y)
                          
                          & ((sample['startz_assigned_muon'] > min_FV_z) & (sample['startz_assigned_muon'] < max_FV_z) )
                          & ((sample['startz_assigned_proton'] > min_FV_z) & (sample['startz_assigned_proton'] < max_FV_z) )
                          & ((sample['endz_assigned_muon'] > min_FV_z) & (sample['endz_assigned_muon'] < max_FV_z) )
                          & ((sample['endz_assigned_proton'] > min_FV_z) & (sample['endz_assigned_proton'] < max_FV_z) )
                          
                          & ((sample['startx_assigned_muon'] > min_FV_x) & (sample['startx_assigned_muon'] < max_FV_x) )
                          & ((sample['startx_assigned_proton'] > min_FV_x) & (sample['startx_assigned_proton'] < max_FV_x) )
                          & ((sample['endx_assigned_muon'] > min_FV_x) & (sample['endx_assigned_muon'] < max_FV_x) )
                          & ((sample['endx_assigned_proton'] > min_FV_x) & (sample['endx_assigned_proton'] < max_FV_x) )
                          ]
    return sample_in_FV
# ------------------------------------------------

#---------------------------------------------------------------------------------------------
# July-11, 2017
def plot_feature_pairs(cut_name='${PID}_A$',
                       var='l_long',x_label='$l_{long}$ [cm]',mul=1,
                       bins=np.linspace(0,300,50),
                       figsize=(12,8),legend_fontsize=25,fontsize=25,
                       do_add_legend=False,legend_loc='upper center',
                       ticks_color='black'):
    fig,ax = plt.subplots(figsize=figsize)
    max_h=0
    text_colors=[]
    for i,(pair_type,label,cmap,color) in enumerate(zip(pair_types[1:],labels[1:],cmaps[1:],colors[1:])):
        sample = reduced_samples[cut_name][pair_type]
        if len(sample) < 10: continue
        h,bins,_=plt.hist(mul*sample[var],normed=1,bins=bins,histtype='step',linewidth=3,color=color)
        text_colors.append(color)
        p = plt.plot([0,0], label=label,linestyle='-',linewidth=6,color=color)
        if np.max(h)>max_h:
            max_h=np.max(h)

    if do_add_legend:
        leg=ax.legend(fontsize=legend_fontsize,loc=legend_loc)
        for text_color,text in zip(text_colors,leg.get_texts()):
            text.set_color(text_color)

    set_axes(ax,x_label=x_label,fontsize=fontsize,ticks_color=ticks_color,do_add_grid=True)
    ax.set_xlim(np.min(bins),np.max(bins))
    ax.set_ylim(0,1.05*max_h)
    ax.xaxis.set_major_locator(LinearLocator(5));ax.yaxis.set_major_locator(LinearLocator(4))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    plt.tight_layout()
    return ax
#---------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------
# July-11, 2017
def plot_cut_samples (reduced_cut_name='${PID}_A$',
                      cut_name='maximal distance between tracks',mul=1,
                      cut_var ='distance',
                      cut_type= 'max',
                      x_label = 'maximal tracks distance [cm]', y_label='% of sample',
                      xcenter=0,figsize=(12,8),fontsize=25,
                      xmin=0.1, xmax=10 , Nbins=10, do_add_legend=True, legend_loc='bbox',legend_fontsize=25,
                      ticks_color='black'):
    fig,ax=plt.subplots(figsize=figsize)
    for i,(pair_type,label,cmap,color) in enumerate(zip(pair_types[1:],labels[1:],cmaps[1:],colors[1:])):
        sample = reduced_samples[reduced_cut_name][pair_type]
        if cut_type=='max' or cut_type=='min':
            x , frac , frac_err = get_fraction_in_cut( data=sample , cut_var=cut_var , mul=mul , cut_type=cut_type , xmin=xmin, xmax=xmax , Nbins=Nbins )
        elif cut_type=='symmetric':
            x , frac , frac_err = get_fraction_in_symmetriccut( data=sample , cut_var=cut_var , mul=mul , xcenter=xcenter, delta_x_min=xmin, delta_x_max=xmax , Nbins=Nbins )
        plt.errorbar(x , y=frac, yerr=frac_err , xerr=0, fmt='o' , markersize=markers_size , label=label, color=color)
    if do_add_legend:
        if 'bbox' not in legend_loc:
            leg=ax.legend(fontsize=legend_fontsize,loc=legend_loc,markerscale=2.)
        else:
            leg=ax.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.,fontsize=legend_fontsize,markerscale=2.)
        for color,text in zip(colors[1:],leg.get_texts()):
            text.set_color(color)
    ax.set_ylim(0,101)
    ax.set_xlim(xmin,xmax)
    set_axes(ax,x_label=x_label,y_label=y_label,fontsize=fontsize,ticks_color=ticks_color,xticks=np.linspace(xmin,xmax,7),yticks=[25,50,75,100])
    ax.grid(linestyle='--',alpha=0.75)
    plt.tight_layout()
    return ax,leg
#---------------------------------------------------------------------------------------------



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

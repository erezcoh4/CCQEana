
import sys; sys.path.insert(0, '../../'); sys.path.insert(0,'../mupClassification/')
from ccqe_notebook_tools import *
import matplotlib.patches as patches


debug = 0
Nevents=dict()
Nevents['OffBeam sof.trig. efficiency'] = 0.04462
Nevents['OnBeam sof.trig. efficiency'] = 0.05135


Nevents['v04 after sof.trig.'] = 378787 # from python scripts/count_events.py on <prod_reco2_extbnb_v8_mcc8_v04_26_04_05_v04>
Nevents['v05 after sof.trig.'] = 1815 # from python scripts/count_events.py on <prod_reco2_extbnb_v8_mcc8_v04_26_04_05_v05>
Nevents['OffBeam after sof.trig.'] = Nevents['v04 after sof.trig.'] + Nevents['v05 after sof.trig.']
Nevents['OffBeam before sof.trig.'] = Nevents['OffBeam after sof.trig.']/Nevents['OffBeam sof.trig. efficiency']


Nevents['OnBeam after sof.trig.'] = 544114 # from python scripts/count_events.py on <prod_reco2_bnb_v8_mcc8>
Nevents['OnBeam before sof.trig.'] = Nevents['OnBeam after sof.trig.']/Nevents['OnBeam sof.trig. efficiency']
Nevents['OnBeam POT'] = 4.93e19

OffBeam_scaling = Nevents['OnBeam before sof.trig.']/Nevents['OffBeam before sof.trig.']
print "OffBeam_scaling:",OffBeam_scaling,"= N(on beam)/N(off beam) before sof. trig."


Nevents['MC-BNB/Cosmic-DATA overlay'] = 96350 # from python scripts/count_events.py on <prodgenie_bnb_nu_uboone_overlay_mcc8_reco2>
Nevents['MC-BNB/Cosmic-DATA overlay POT'] = 9.773e19
MC_scaling = Nevents['OnBeam POT']/Nevents['MC-BNB/Cosmic-DATA overlay POT']
print "MC_scaling:",MC_scaling,"= N(POT on beam)/N(POT MC)"


# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Sep-3,2017
def plot_stacked_MCsamples( ax=None, MCsamples = None
                           , var=None, bins=None , N_OnBeam_minus_OffBeam=1, alpha=0.8):
    Nall_pairs = len(MCsamples['1mu-1p']+MCsamples['cosmic']+MCsamples['other pairs'])
    #     MC_norm_fact = N_OnBeam_minus_OffBeam/Nall_pairs
    global debug
    x_array, weights_array = [] , []
    label_array , color_array = [] , []
    # stack background (cosmic, other-pairs) + 1u1p pairs
    for i_pair_type in [2,1,0]:
        pair_type=pair_types[i_pair_type]
        sample = MCsamples[pair_type];
        label_array.append(MClabels[i_pair_type]);
        color_array.append(MCcolors[i_pair_type]);
        x = sample[var]
        x_array.append(x)
        # normalize the MC to have the same number of events as the total On-Off beam sample
        weights_array.append (MC_scaling * np.ones(len(x)) )
    # -- - - - --------- - - -- ---- -  - --- -- -- -- --
    bin_width = bins[1]-bins[0]
    h,bins_arr,_=ax.hist( x_array , weights=weights_array
                         , bins=bins-0.5*bin_width, width=bin_width
                         , stacked=True
                         , color=color_array
                         , label=label_array
                         , alpha=alpha)
    h_stack = h[0]+h[1]+h[2]
    if np.max(h[2])>np.max(ax.get_ylim()):#{
        ax.set_ylim(np.min(ax.get_ylim()),1.05*np.max(h[2]))
    #}
    # add CC1p0pi as a box inside the 1u1p
    sample = MCsamples[pair_types[3]]
    hCC1p0pi,edges = np.histogram( sample[var] , weights=MC_scaling*np.ones(len(sample)) , bins=bins )
    for bin in range(len(bins[:-2])):#{
        x, dx = bins[bin+1] - 0.4*bin_width, bin_width
        y, dy = 0.99*h[2][bin+1] - hCC1p0pi[bin], hCC1p0pi[bin]
        ax.add_patch( patches.Rectangle( (x, y),dx,dy, facecolor=MCcolors[3],alpha=0.8*alpha,label=MClabels[3] if bin==0 else None))
    #}
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -



# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Sep-3,2017
def OnBeam_minus_OffBeam_1d( OnBeamSample=None , OffBeamSample=None, MCsamples=None
                            , var='PIDa_assigned_proton' , x_label='$PID_a^p$'
                            , bins=np.linspace(0,30,31)
                            , ax=None, figsize=(14,6),fontsize=25
                            , color='purple'
                            , do_add_MCoverlay=True
                            , do_add_legend=True , legend_loc='best', MCalpha=0.5):
    global debug
    bin_width = bins[1]-bins[0]
    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    h_OnBeam,edges = np.histogram( OnBeamSample[var] , bins=bins )
    h_OnBeam_err = np.sqrt(h_OnBeam)
    h_OffBeam,edges = np.histogram( OffBeamSample[var] , bins=bins )
    h_OffBeam_err = np.sqrt(h_OffBeam)
    
    h_OnBeam_minus_OffBeam = h_OnBeam - OffBeam_scaling*h_OffBeam
    h_OnBeam_minus_OffBeam_err = np.sqrt( np.square(h_OnBeam_err) + np.square(OffBeam_scaling*h_OffBeam_err)  )
    
    plt.errorbar( x = bins[:-1], xerr=bin_width/2.
                 , y=h_OnBeam_minus_OffBeam , yerr=h_OnBeam_minus_OffBeam_err
                 , fmt='o', color=color , ecolor='black', label='(On-Off) Beam'
                 )
    ax.set_xlim(np.min(bins)-bin_width,np.max(bins)+bin_width);
    ax.set_ylim(np.min([0,np.min(h_OnBeam_minus_OffBeam-1.1*h_OnBeam_minus_OffBeam_err)])
                             ,np.max(h_OnBeam_minus_OffBeam+1.1*h_OnBeam_minus_OffBeam_err));
                 
    plt.plot(ax.get_xlim(),[0,0],'--',color='black',linewidth=2)
    set_axes(ax,x_label=x_label,y_label='counts',do_add_grid=True,fontsize=fontsize)
                 
    if do_add_MCoverlay:#{
        N_OnBeam_minus_OffBeam = len(OnBeamSample) - OffBeam_scaling*len(OffBeamSample)
        if debug: print 'Number of On-Off:',N_OnBeam_minus_OffBeam
        plot_stacked_MCsamples( ax=ax
                               , MCsamples = MCsamples , var=var, bins=bins
                               , N_OnBeam_minus_OffBeam=N_OnBeam_minus_OffBeam , alpha=MCalpha)
    #}
    if do_add_legend:
        plt.legend(fontsize=fontsize,loc=legend_loc)
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -


# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Sep-3,2017
def OnBeam_minus_OffBeam_2d( OnBeamSample=None , OffBeamSample=None
                            , varx='l_assigned_proton' , x_label='$l^{\\mu}$'
                            , vary='l_assigned_muon' , y_label='$l^p$'
                            , bins=(np.linspace(0,100,51),np.linspace(0,150,51))
                            , ax=None, figsize=(12,8),fontsize=25
                            , cmap='hot_r'):
    global debug
    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    h_OnBeam_T,xedges, yedges = np.histogram2d( OnBeamSample[varx] , OnBeamSample[vary] , bins=bins )
    h_OffBeam_T,xedges, yedges = np.histogram2d( OffBeamSample[varx] , OffBeamSample[vary] , bins=bins )
    
    h_OnBeam_minus_OffBeam_T = h_OnBeam_T - OffBeam_scaling*h_OffBeam_T
    h_OnBeam_minus_OffBeam = h_OnBeam_minus_OffBeam_T.T
    
    X, Y = np.meshgrid(xedges, yedges)
    pcmesh = ax.pcolormesh(X, Y, h_OnBeam_minus_OffBeam ,cmap=cmap)
    cbar = plt.colorbar(pcmesh,label='(On-Off) Beam [cts]')
    cbar.ax.yaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(size=fontsize))
    cbar.ax.tick_params(labelsize=fontsize)
    
    bin_width = [bins[0][1]-bins[0][0],bins[1][1]-bins[1][0]]
    ax.set_xlim(np.min(bins[0])-bin_width[0],np.max(bins[0])+bin_width[0]);
    ax.set_ylim(np.min(bins[1])-bin_width[1],np.max(bins[1])+bin_width[1]);
    set_axes(ax,x_label=x_label,y_label=y_label,do_add_grid=True,fontsize=fontsize)
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -

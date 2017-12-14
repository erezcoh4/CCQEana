
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
MC_scaling_DATAcosmic = Nevents['OnBeam POT']/Nevents['MC-BNB/Cosmic-DATA overlay POT']
print "MC_scaling_DATAcosmic:",MC_scaling_DATAcosmic,"= N(POT on beam)/N(POT MC)"

Nevents['MC-BNB/Cosmic-MC overlay'] = 358800 # from python scripts/count_events.py on <prodgenie_bnb_nu_cosmic_uboone_mcc8.2_reco2>
Nevents['MC-BNB/Cosmic-MC overlay POT'] = 3.61901e20
MC_scaling_MCcosmic = Nevents['OnBeam POT']/Nevents['MC-BNB/Cosmic-MC overlay POT']
print "MC_scaling_MCcosmic:",MC_scaling_MCcosmic,"= N(POT on beam)/N(POT MC)"

OnBeamColor = 'teal'
OffBeamColor = 'purple'


# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Dec-6,2017
def plot_OnBeam(OnBeamSample=None,OnBeamFV=None
                , var='PIDa_assigned_proton' , x_label='$PID_a^p$'                 
                , bins=np.linspace(0,30,31)                 
                , ax=None, figsize=(14,6),fontsize=25                
                , color=OnBeamColor
                , do_add_legend=True , legend_loc='best'):
    bin_width = bins[1]-bins[0]
    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    x = OnBeamSample[var]
    h_OnBeam,edges = np.histogram( x , bins=bins )
    h_OnBeam_err = np.sqrt(h_OnBeam)
    
    plt.errorbar( x = bins[:-1], xerr=bin_width/2., markersize=12
                 , y=h_OnBeam , yerr=h_OnBeam_err
                 , fmt='o', color=color , ecolor=color
                 , label='BNB (%d=%.1f'%(len(OnBeamSample),100*float(len(OnBeamSample))/len(OnBeamFV))+'%)'
                )
    plt.plot([0,0],[0,0],'--',color='black',linewidth=2)
    
    set_axes(ax,x_label=x_label,y_label='counts',do_add_grid=True,fontsize=fontsize
             ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
            )
    if do_add_legend: 
        if legend_loc=='bbox':
            leg=plt.legend(bbox_to_anchor=(1.,1.05),fontsize=fontsize,loc=2)
        else:
            leg=plt.legend(fontsize=fontsize,loc=legend_loc)
    plt.tight_layout()
    return ax,leg
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -

# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Dec-6,2017
def plot_OffBeam(OffBeamSample=None,OffBeamFV=None
                , var='PIDa_assigned_proton' , x_label='$PID_a^p$'                 
                , bins=np.linspace(0,30,31)                 
                , ax=None, figsize=(14,6),fontsize=25                
                , color=OffBeamColor
                , do_add_legend=True , legend_loc='best'
                 , do_OffBeam_scaling=True):
    bin_width = bins[1]-bins[0]
    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    x = OffBeamSample[var]
    
    h_OffBeam,edges = np.histogram( x , bins=bins )
    h_OffBeam_err = np.sqrt(h_OffBeam)
    
    plt.errorbar( x = bins[:-1], xerr=bin_width/2., markersize=12
                 , y=OffBeam_scaling*h_OffBeam if do_OffBeam_scaling else h_OffBeam
                 , yerr=OffBeam_scaling*h_OffBeam_err if do_OffBeam_scaling else h_OffBeam_err
                 , fmt='s', color=color , ecolor=color
                 , label='extBNB (%.1f=%.1f'%(OffBeam_scaling*len(OffBeamSample),100*float(len(OffBeamSample))/len(OffBeamFV))+'%)'
                )
    plt.plot([0,0],[0,0],'--',color='black',linewidth=2)
    
    set_axes(ax,x_label=x_label,y_label='counts',do_add_grid=True,fontsize=fontsize
             ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
            )
    if do_add_legend: 
        if legend_loc=='bbox':
            leg=plt.legend(bbox_to_anchor=(1.,1.05),fontsize=fontsize,loc=2)
        else:
            leg=plt.legend(fontsize=fontsize,loc=legend_loc)
    plt.tight_layout()
    return ax,leg
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -

# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Dec-4,2017
def draw_var(cut_name=None,var=None,x_label=None,bins=None,debug=0
             ,OnBeamFV=None,OffBeamFV=None
             ,do_cosmic_only=True,do_bnb_only=True
             ,chi2_xrange=None,chi2_xy=None
             ,reduced_OnBeam=None,reduced_OffBeam=None
             ,MCbnbDATAcosmicSamples=None, reduced_MCbnbDATAcosmicSamples=None
             ,OriginalOnBeamSample=None , OriginalOffBeamSample=None
             ,do_save_fig=True,figures_path=None):

    Nsubplots = 1
    i_subplot = 0
    if do_cosmic_only: Nsubplots += 1
    if do_bnb_only: Nsubplots += 1
        
    fig = plt.figure(figsize=(16,6*Nsubplots))        
    if do_cosmic_only:
        i_subplot+=1
        ax = fig.add_subplot(Nsubplots,1,i_subplot)
        ax,leg=extBNBvsCosmicOverlay(OffBeamSample=reduced_OffBeam[cut_name]
                                     ,OffBeamFV=OffBeamFV,MCbnbDATAcosmicSamples=MCbnbDATAcosmicSamples
                                     ,cosmic_overlay_sample = reduced_MCbnbDATAcosmicSamples[cut_name]['cosmic']
                                     ,var=var , color='black' ,x_label=x_label, bins=bins , ax=ax, legend_loc='bbox')
    if do_bnb_only:
        i_subplot+=1
        ax = fig.add_subplot(Nsubplots,1,i_subplot)
        plot_OnBeam(OnBeamSample=reduced_OnBeam[cut_name],OnBeamFV=OnBeamFV
                                     ,var=var , x_label=x_label, bins=bins , ax=ax, legend_loc='bbox')
        plot_OffBeam(OffBeamSample=reduced_OffBeam[cut_name],OffBeamFV=OffBeamFV
                                     ,var=var , x_label=x_label, bins=bins , ax=ax, legend_loc='bbox')


    ax = fig.add_subplot(Nsubplots,1,Nsubplots)
    ax,leg=OnBeam_minus_OffBeam_1d(debug=debug
                                   ,OnBeamSample=reduced_OnBeam[cut_name] 
                                   ,OffBeamSample=reduced_OffBeam[cut_name] 
                                   ,MCsamples=reduced_MCbnbDATAcosmicSamples[cut_name]
                                   ,MCbnbDATAcosmicSamples=MCbnbDATAcosmicSamples
                                   ,var=var
                                   ,x_label=x_label 
                                   ,bins=bins 
                                   ,ax=ax
                                   ,legend_loc='bbox'                                       
                                   ,do_add_chi2_MC_data=True , chi2_xrange=chi2_xrange, chi2_xy=chi2_xy
                                   ,OriginalOnBeamSample=reduced_OnBeam['no cut'] , OriginalOffBeamSample=reduced_OffBeam['no cut'])
    if do_save_fig:
        plt.savefig(figures_path+var+'_'+'after_cut_'+cut_name+'.pdf', bbox_inches='tight')    
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -





# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Nov-9,2017
def extBNBvsCosmicOverlay(OffBeamSample=None,OffBeamFV=None
                            , var='PIDa_assigned_proton' , x_label='$PID_a^p$' 
                            , bins=np.linspace(0,30,31) 
                            , ax=None, figsize=(14,6),fontsize=25
                            , color='purple'
                            , do_add_cosmic_overlay=True , cosmic_overlay_sample=None, MCbnbDATAcosmicSamples=None
                            , do_add_legend=True , legend_loc='best', overlay_alpha=0.5):
    bin_width = bins[1]-bins[0]
    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    x = OffBeamSample[var]
    h_OffBeam,edges = np.histogram( x , bins=bins )
    h_OffBeam_err = np.sqrt(h_OffBeam)
    Int_OffBeam = np.sum([h_OffBeam[i]*bin_width for i in range(len(h_OffBeam))])
    
    plt.errorbar( x = bins[:-1], xerr=bin_width/2.
                 , y=h_OffBeam , yerr=h_OffBeam_err
                 , fmt='o', color=color , ecolor='black'
                 , label='extBNB (%d=%.1f'%(len(OffBeamSample),100*float(len(OffBeamSample))/len(OffBeamFV))+'%)'
                )
    plt.plot([0,0],[0,0],'--',color='black',linewidth=2)
    
    if do_add_cosmic_overlay:        
        plot_cosmic_overlay( ax=ax                               
                            , cosmic_overlay_sample = cosmic_overlay_sample , var=var, bins=bins 
                            , Int_OffBeam = Int_OffBeam
                            , alpha=overlay_alpha,MCbnbDATAcosmicSamples=MCbnbDATAcosmicSamples)
    set_axes(ax,x_label=x_label,y_label='counts',do_add_grid=True,fontsize=fontsize
             ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
    #              ,ylim=(0,1.05*np.max(h_OffBeam))
            )
    if do_add_legend: 
        if legend_loc=='bbox':
            leg=plt.legend(bbox_to_anchor=(1.,1.05),fontsize=fontsize,loc=2)
        else:
            leg=plt.legend(fontsize=fontsize,loc=legend_loc)
    plt.tight_layout()
    return ax,leg
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -


# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Nov-9,2017
def plot_cosmic_overlay( ax=None, cosmic_overlay_sample=None, var=None, bins=None
                        , alpha=0.5, Int_OffBeam=1,MCbnbDATAcosmicSamples=None):
    x = cosmic_overlay_sample[var]
    bin_width = bins[1]-bins[0]
        
    h_overlay,edges = np.histogram( x , bins=bins )
    Int_overlay = np.sum([h_overlay[i]*bin_width for i in range(len(h_overlay))])
    SUMoverlay = np.sum(h_overlay)    

    cosmic_overlay_scaling = Int_OffBeam/Int_overlay    
    ax.hist( x , weights=cosmic_overlay_scaling*np.ones(len(x))
                     , bins=bins-0.5*bin_width, width=bin_width
                     , color=MCcolors[2]
                     , label='cosmic overlay (%d=%.1f'%(len(cosmic_overlay_sample),100*float(len(cosmic_overlay_sample))/len(MCbnbDATAcosmicSamples['cosmic']))+'%)'
                     , alpha=alpha)       
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -        



# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# written Nov-9,2017  (last edit Dec-6)
def OnBeam_minus_OffBeam_1d( OnBeamSample=None , OffBeamSample=None , debug=0
                            , var='PIDa_assigned_proton' , x_label='$PID_a^p$' 
                            , bins=np.linspace(0,30,31) 
                            , ax=None, figsize=(14,6),fontsize=25
                            , color='tomato'
                            , do_add_MCoverlay=True , MCsamples=None, MCbnbDATAcosmicSamples=None
                            , MC_scaling=MC_scaling_DATAcosmic
                            , do_add_legend=True , legend_loc='best', MCalpha=0.7
                            , do_add_chi2_MC_data=False , chi2_xrange=None, chi2_xy=(0,0)
                            , OriginalOnBeamSample=None , OriginalOffBeamSample=None):
    bin_width = bins[1]-bins[0]
    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    h_OnBeam,edges = np.histogram( OnBeamSample[var] , bins=bins )
    h_OnBeam_err = np.sqrt(h_OnBeam)
    h_OffBeam,edges = np.histogram( OffBeamSample[var] , bins=bins )
    h_OffBeam_err = np.sqrt(h_OffBeam)
    
    h_OnBeam_minus_OffBeam = h_OnBeam - OffBeam_scaling*h_OffBeam
    h_OnBeam_minus_OffBeam_err = np.sqrt( np.abs(h_OnBeam + OffBeam_scaling*OffBeam_scaling*h_OffBeam) )
    Integral = len(OnBeamSample) - OffBeam_scaling*len(OffBeamSample)
    Integral_Original = len(OriginalOnBeamSample) - OffBeam_scaling*len(OriginalOffBeamSample)

    
    plt.errorbar( x = bins[:-1], xerr=bin_width/2., markersize=12
                 , y=h_OnBeam_minus_OffBeam , yerr=h_OnBeam_minus_OffBeam_err
                 , fmt='o', color=color , ecolor='black', label=r'(On-Off) Beam ($\int=$%.1f=%.1f'%(Integral,100*Integral/Integral_Original)+'%)'
                )
    if debug>1: print "OnBeam-OffBeam (bins[:-1]):\n",bins[:-1]
    plt.plot([np.min(ax.get_xlim()),np.min(ax.get_xlim())],[0,0],'--',color='black',linewidth=2)
    set_axes(ax,x_label=x_label,y_label='counts',do_add_grid=True,fontsize=fontsize
             ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
             ,ylim=(np.min([0,np.min(h_OnBeam_minus_OffBeam-1.1*h_OnBeam_minus_OffBeam_err)])               
                    ,np.max(h_OnBeam_minus_OffBeam+1.1*h_OnBeam_minus_OffBeam_err))
            )
    
    if do_add_MCoverlay:        
        h_MC, bins_MC = plot_stacked_MCsamples( ax=ax, debug=debug
                               , MCsamples = MCsamples , MCbnbDATAcosmicSamples=MCbnbDATAcosmicSamples
                               , var=var, bins=bins , MC_scaling=MC_scaling
                               , alpha=MCalpha)
    if do_add_chi2_MC_data:
        chi2 , ndf = chi2_two_histograms( bins=bins, chi2_xrange=chi2_xrange
                                         , h1=h_MC , h2=h_OnBeam_minus_OffBeam 
                                         , h1err=np.sqrt(h_MC), h2err=h_OnBeam_minus_OffBeam_err
                                         , debug=debug )
        plt.text(chi2_xy[0],chi2_xy[1],r'$\chi^2/ndf=%.1f/%d$'%(chi2,ndf),fontsize=fontsize)
        if chi2_xrange is not None:
            plt.plot([chi2_xrange[0]-0.5*bin_width,chi2_xrange[0]-0.5*bin_width],ax.get_ylim(),'-',alpha=0.3)
            plt.plot([chi2_xrange[1]-0.5*bin_width,chi2_xrange[1]-0.5*bin_width],ax.get_ylim(),'-',alpha=0.3)

    if do_add_legend: 
        if legend_loc=='bbox':
            leg=plt.legend(bbox_to_anchor=(1.,1.05),fontsize=fontsize,loc=2)
        else:
            leg=plt.legend(fontsize=fontsize,loc=legend_loc)
    plt.tight_layout()
    return ax,leg
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -



# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Nov-20,2017 (last editted Dec-6)
def plot_stacked_MCsamples( ax=None, debug=0
                           , MCsamples=None
                           , MC_scaling=MC_scaling_DATAcosmic
                           , MCbnbDATAcosmicSamples=None
                           , var=None, bins=None , alpha=0.8):
    '''
    return: h, bins
            stacked histogram values and bins,
            of the samples from the overlay: 
            (cosmic, other-pairs) + 1mu 1p pairs
    '''
    bin_width = bins[1]-bins[0]    
    h=dict()
    labels = dict()
    colors = dict()
    for i_pair_type,pair_type in enumerate(pair_types):
        sample = MCsamples[pair_type]; 
        labels[pair_type] = MClabels[i_pair_type]+' (%.1f=%.1f'%(MC_scaling*len(sample),100*float(len(sample))/len(MCbnbDATAcosmicSamples[pair_type]))+'%)'
        colors[pair_type] = MCcolors[i_pair_type];
        x = sample[var]
        histo,edges = np.histogram(x,bins=bins)
        h[pair_type] = MC_scaling*histo 
    # -- - - - --------- - - -- ---- -  - --- -- -- -- --
    if debug>1: print "stacked_MCsamples (bins[:-1]):\n",bins[:-1]
    plt.bar(bins[:-1]-0.5*bin_width,h['cosmic']+h['other pairs']+h['1mu-1p'] , width=bin_width
                     ,color=colors['1mu-1p'],alpha=alpha, label=labels['1mu-1p'])
    # CC 1p 0pi
    plt.bar(bins[:-1]-0.5*bin_width,h['cosmic']+h['other pairs']+h['CC 1p 0pi'] , width=bin_width
                     ,color=colors['CC 1p 0pi'],alpha=alpha, label=labels['CC 1p 0pi'])
    plt.bar(bins[:-1]-0.5*bin_width,h['cosmic']+h['other pairs'] , width=bin_width
                     ,color=colors['other pairs'],alpha=alpha , label=labels['other pairs'])
    plt.bar(bins[:-1]-0.5*bin_width,h['cosmic'] , width=bin_width
                     ,color=colors['cosmic'],alpha=alpha, label=labels['cosmic'])
    h_stack = h['cosmic']+h['other pairs']+h['1mu-1p']
    if np.max(h_stack)>np.max(ax.get_ylim()): ax.set_ylim(np.min(ax.get_ylim()),1.05*np.max(h_stack))
    return h_stack , bins
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -



# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Dec-3, 2017 (last edit Dec-9)
def chi2_two_histograms( bins=None, chi2_xrange=None
                        , h1=None , h2=None 
                        , h1err=None, h2err=None 
                        , debug=0):
    '''
    compare the two histograms using a chi2 test.
    return: chi2, ndf
    '''
    chi2 = 0
    Nbins_compares = 0
    for i_bin in range(len(bins)-1):
        if chi2_xrange is None or (bins[i_bin]>=chi2_xrange[0] and bins[i_bin]<=chi2_xrange[1]):
            if debug: 
                print 'comparing in bin:',bins[i_bin],'h1:',h1[i_bin],'+/-',h1err[i_bin],'h2:',h2[i_bin],'+/-',h2err[i_bin]
                print 'chi2+=',(np.square( h1[i_bin] - h2[i_bin] ) / np.max([( np.square(h1err[i_bin]) + np.square(h2err[i_bin]) ),1]))
            chi2 += np.square( h1[i_bin] - h2[i_bin] ) / np.max([( np.square(h1err[i_bin]) + np.square(h2err[i_bin]) ),1])
            Nbins_compares += 1
    ndf = Nbins_compares - 1
    if debug: print 'chi2,ndf:',chi2,ndf
    return chi2,ndf
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -







import matplotlib.colors as colors
# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)
    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    
    

# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Sep-3,2017 (last edit Dec-9)
def OnBeam_minus_OffBeam_2d( OnBeamSample=None , OffBeamSample=None , debug=0
                            , varx='l_assigned_proton' , x_label='$l^{\\mu}$'
                            , vary='l_assigned_muon' , y_label='$l^p$'
                            , bins=(np.linspace(0,100,51),np.linspace(0,150,51))
                            , ax=None, figsize=(12,8),fontsize=25
                            , cmap='hot_r'):
    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    h_OnBeam_T,xedges, yedges = np.histogram2d( OnBeamSample[varx] , OnBeamSample[vary] , bins=bins )
    h_OffBeam_T,xedges, yedges = np.histogram2d( OffBeamSample[varx] , OffBeamSample[vary] , bins=bins )
    
    h_OnBeam_minus_OffBeam_T = h_OnBeam_T - OffBeam_scaling*h_OffBeam_T
    h_OnBeam_minus_OffBeam = h_OnBeam_minus_OffBeam_T.T
    
    X, Y = np.meshgrid(xedges, yedges)
    elev_min, elev_max = np.min(h_OnBeam_minus_OffBeam) , np.max(h_OnBeam_minus_OffBeam)
    pcmesh = ax.pcolormesh(X, Y, h_OnBeam_minus_OffBeam ,cmap=cmap
                           , clim=(elev_min, elev_max), norm=MidpointNormalize(midpoint=0,vmin=elev_min, vmax=elev_max))
    cbar = plt.colorbar(pcmesh,label='(On-Off) Beam [cts]')
    cbar.ax.yaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(size=fontsize))
    cbar.ax.tick_params(labelsize=fontsize)
    
    bin_width = [bins[0][1]-bins[0][0],bins[1][1]-bins[1][0]]
    ax.set_xlim(np.min(bins[0])-bin_width[0],np.max(bins[0])+bin_width[0]);
    ax.set_ylim(np.min(bins[1])-bin_width[1],np.max(bins[1])+bin_width[1]);
    set_axes(ax,x_label=x_label,y_label=y_label,do_add_grid=True,fontsize=fontsize)
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -


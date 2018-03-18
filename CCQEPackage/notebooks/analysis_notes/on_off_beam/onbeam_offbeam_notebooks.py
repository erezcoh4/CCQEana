
import sys; sys.path.insert(0, '../../'); sys.path.insert(0,'../mupClassification/')
from ccqe_notebook_tools import *
from numpy import sqrt,square
import matplotlib.patches as patches


debug = 0
Nevents=dict()
# following docdb 5640 [https://microboone-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=5640&filename=MCC8_Data_factors.pdf&version=5]
# we take the normalization factors
#Nevents['OffBeam sof.trig. efficiency'] = 0.04462
#Nevents['OnBeam sof.trig. efficiency'] = 0.05135
#Nevents['v04 after sof.trig.'] = 378787 # from python scripts/count_events.py on <prod_reco2_extbnb_v8_mcc8_v04_26_04_05_v04>
#Nevents['v05 after sof.trig.'] = 1815 # from python scripts/count_events.py on <prod_reco2_extbnb_v8_mcc8_v04_26_04_05_v05>
#Nevents['OffBeam after sof.trig.'] = Nevents['v04 after sof.trig.'] + Nevents['v05 after sof.trig.']
#Nevents['OffBeam before sof.trig.'] = Nevents['OffBeam after sof.trig.']/Nevents['OffBeam sof.trig. efficiency']
#Nevents['OnBeam after sof.trig.'] = 544114 # from python scripts/count_events.py on <prod_reco2_bnb_v8_mcc8>
#Nevents['OnBeam before sof.trig.'] = Nevents['OnBeam after sof.trig.']/Nevents['OnBeam sof.trig. efficiency']
#OffBeam_scaling = Nevents['OnBeam before sof.trig.']/Nevents['OffBeam before sof.trig.']
# MC-BNB/Cosmic-MC overlay
#Nevents['MC-BNB/Cosmic-MC overlay'] = 358800 # from python scripts/count_events.py on <prodgenie_bnb_nu_cosmic_uboone_mcc8.2_reco2>
#Nevents['MC-BNB/Cosmic-MC overlay POT'] = 3.61901e20
#MC_scaling_MCcosmic = Nevents['OnBeam POT']/Nevents['MC-BNB/Cosmic-MC overlay POT']
#print "MC_scaling_MCcosmic:",MC_scaling_MCcosmic,"= N(POT on beam)/N(POT MC)"

'''
    Off Beam scaling factor (Ariana, Feb 2018)
    -- ---- - ---- - -- -- --
    
    ** do not source localProducts
    ** > setup sam_web_client v2_1
    
    0) Grab getDataInfo.py from here:
    /uboone/app/users/zarko
    
    More info on the tool as of last fall can be found here:
    [https://microboone-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=11078&filename=pot_and_beamq.pdf&version=1]
    
    1) Feed the script a file-list to calculate the trigger content of each sample we run over.  There are some other options other than file-list that you can feed the script if you don't have this info; take a look at the options in the script to understand which is best for the info you have.
    
    2) Trigger counts for various triggers is returned by the script in the following form :
    EXT	Gate2	E1DCNT	tor860	tor875	E1DCNT_wcut	tor860_wcut	tor875_wcut
    
    To calculate the scaling factor of off-to-on beam,
    take the ratio of E1DCNT_wcut (from OnBeam sample) / EXT (from OffBeam sample).
    (
    Zarko:
    for MCC8 you need to add -v2 option.
    Or for any sample that includes runs after ~6300
    those samples have to be processed with new beam data quality cuts v2 otherwise they'll all fail beam cuts
    )
    for example: MCC8.5 reprocessing
    -- ---- - ---- -
    
    > python scripts/getDataInfo.py --defname="prod_reco_optfilter_bnb_v11_unblind_mcc8" -v2
    
    Definition prod_reco_optfilter_bnb_v11_unblind_mcc8 contains 4110 files
           EXT         Gate2        E1DCNT        tor860        tor875   E1DCNT_wcut   tor860_wcut   tor875_wcut
       5888530      11578672      11583581     4.963e+19     4.957e+19      10948876     4.909e+19     4.903e+19

    > python scripts/getDataInfo.py --defname="prod_reco_optfilter_extbnb_v11_mcc8_dev" -v2
    
    Definition prod_reco_optfilter_extbnb_v11_mcc8_dev contains 5789 files
           EXT         Gate2        E1DCNT        tor860        tor875   E1DCNT_wcut   tor860_wcut   tor875_wcut
      15499028      22044264      22056782     9.146e+19     9.135e+19      20128930     9.057e+19     9.045e+19
    Warning!! BNB data for some of the requested runs/subruns is not in the database.
    Runs missing BNB data (number of subruns missing the data): 5762 (1),
    (Zarko, Mar-02, 2018: You can ignore the warning)
    '''
OffBeam_scaling = float(10948876)/15499028
Nevents['OnBeam POT'] = 4.957e+19
print "OffBeam_scaling:",OffBeam_scaling,"= N(on beam)/N(off beam) before sof. trig."

# MC-BNB/Cosmic-DATA overlay
summary = pd.read_csv('/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/summary/ecohen_physical_files_adi_prodgenie_bnb_nu_uboone_overlay_cosmic_data_100K_reco2_2018_02_23_summary.csv')
Nevents['MC-BNB/Cosmic-DATA overlay'] = np.sum(summary.Nevents)
Nevents['MC-BNB/Cosmic-DATA overlay POT'] = np.sum(summary.POT)
MC_scaling_DATAcosmic = Nevents['OnBeam POT']/Nevents['MC-BNB/Cosmic-DATA overlay POT']
print "MC_scaling_DATAcosmic:",MC_scaling_DATAcosmic,"= N(POT on beam)/N(POT MC)"

OnBeamColor = 'teal'
OffBeamColor = 'orange'


# ------------------------------------------------
# March-6, 2018
def get_Nreduced(MCbnbDATAcosmicSamples=None,reduced = dict()):
    Noriginal , Nreduced , freduced = dict() , dict() , dict()
    for pair_type in pair_types:
        sam = MCbnbDATAcosmicSamples[pair_type]
        Noriginal[pair_type] = len(MCbnbDATAcosmicSamples[pair_type])
        Nreduced[pair_type] = float(len(reduced[pair_type]))
        freduced[pair_type] = 100.0 * Nreduced[pair_type]/Noriginal[pair_type]
    return Nreduced , freduced
# ------------------------------------------------

# ------------------------------------------------
# March 5, 2018
def get_pureff_cut(MCbnbDATAcosmicSamples=None,reduced=None,pureff=None, cut_name = 'PIDa'):
    eff,pur = dict(),dict()
    Nreduced , freduced = get_Nreduced(MCbnbDATAcosmicSamples=MCbnbDATAcosmicSamples,reduced=reduced)
    Ntot = (Nreduced['1mu-1p']+Nreduced['cosmic']+Nreduced['other pairs'])
    eff['1mu-1p'] = freduced['1mu-1p']
    pur['1mu-1p'] = 100.*Nreduced['1mu-1p']/Ntot if Ntot>0 else 0
    eff['CC 1p 0pi'] = freduced['CC 1p 0pi']
    pur['CC 1p 0pi'] = 100.*Nreduced['CC 1p 0pi']/Ntot if Ntot>0 else 0
    pureff_cut = pd.DataFrame({'label':cut_name
                              ,'$\mu p$ eff.':'%.1f'%eff['1mu-1p']+'%'
                              ,'$\mu p$ pur.':'%.1f'%pur['1mu-1p']+'%'
                              ,'CC$0\pi 1 p$ eff.':'%.1f'%freduced['CC 1p 0pi']+'%'
                              ,'CC$0\pi 1 p$ pur.':'%.1f'%(100.*Nreduced['CC 1p 0pi']/Ntot if Ntot>0 else 0)+'%'}
                              , index=[cut_name]
                              )

    for pair_type in pair_types: pureff_cut[pair_type] = '%.1f'%freduced[pair_type]+'%' +' (%.0f)'%Nreduced[pair_type]
    pureff = pureff.append(pureff_cut)
    return pureff
# ------------------------------------------------


# ------------------------------------------------
# March-6, 2018
def apply_cuts_to_overlay(MCbnbDATAcosmicSamples=None
                          ,PIDa_p_min=13
                          ,minPEcut = 100
                          ,maxdYZcut = 200
                          ,delta_theta_12=55  # deg.
                          ,opt_box=(50,100) # [Nwires x Nticks]
                          ,r_max_RdQ_CC1p0pi = 0.35 # sphere in U,V,Y space, apply a cut only to CC1p0pi
                          ,delta_Delta_phi=35 # deg.
                          ,Pt_max=0.35        # GeV/c
                          ,cuts_order=['no cut','PIDa','flash','length','non-collinearity','vertex activity','delta phi','soft Pt']
                          ):
    reducedSamples = dict()
    pureffOverlay = pd.DataFrame()
    
    reducedSamples['no cut'] = dict()
    for pair_type in pair_types: reducedSamples['no cut'][pair_type] = MCbnbDATAcosmicSamples[pair_type]
    pureffOverlay = get_pureff_cut(MCbnbDATAcosmicSamples=MCbnbDATAcosmicSamples
                                   ,pureff=pureffOverlay,cut_name='no cut',reduced=reducedSamples['no cut'])
    
    for i_cut,cut in zip(range(1,len(cuts_order)),cuts_order[1:]):#{
        reduced = dict()
        print 'grabbing reduced samples after (',cuts_order[i_cut-1],') and applying cut on (',cuts_order[i_cut],')'
        samples_previous_cut = reducedSamples[cuts_order[i_cut-1]]
        for pair_type in pair_types:#{
            sam = samples_previous_cut[pair_type]
            if cut == 'PIDa':
                reduced[pair_type] = sam[sam['PIDa_assigned_proton']>PIDa_p_min]
            
            elif cut == 'flash':
                reduced[pair_type] = sam[(sam['Nflashes']>0)
                                         &(sam['ClosestFlash_TotalPE'] > minPEcut)
                                         &(sam['ClosestFlash_YZdistance'] < maxdYZcut)]
            
            elif cut == 'length':
                reduced[pair_type] = sam[sam['PIDa_long'] < sam['PIDa_short']]
    
            elif cut == 'non-collinearity':
                reduced[pair_type] = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
            
            elif cut == 'vertex activity':
                R_str = 'RdQaroundVertex'
                box_str='[%d wires x %d ticks]'%(opt_box[0],opt_box[1])
                Ru,Rv,Ry = R_str+'[plane 0]'+box_str,R_str+'[plane 1]'+box_str,R_str+'[plane 2]'+box_str
                reduced[pair_type] = sam[(sam[Ru]==1) | (sam[Rv]==1) | (sam[Ry]==1)
                                         |
                                         (sqrt( np.square(sam[Ru]-1) + square(sam[Rv]-1) + square(sam[Ry]-1) )
                                          <= r_max_RdQ_CC1p0pi) ]
            elif cut == 'delta phi':
                reduced[pair_type] = sam[np.abs(sam['delta_phi']-180.)<delta_Delta_phi]

            elif cut == 'soft Pt':
                reduced[pair_type] = sam[sam['reco_Pt']<Pt_max]
        #}
        reducedSamples[cut] = reduced
        pureffOverlay = get_pureff_cut(MCbnbDATAcosmicSamples=MCbnbDATAcosmicSamples
                                       ,pureff=pureffOverlay,cut_name=cut,reduced=reduced)
    #}
    return reducedSamples,pureffOverlay
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -


# ------------------------------------------------
# March-6, 2018
def apply_cuts_to_data(PIDa_p_min=13
                       ,minPEcut = 100
                       ,maxdYZcut = 200
                       ,delta_theta_12=55  # deg.
                       ,opt_box=(50,100) # [Nwires x Nticks]
                       ,r_max_RdQ_CC1p0pi = 0.35 # sphere in U,V,Y space, apply a cut only to CC1p0pi
                       ,delta_Delta_phi=35 # deg.
                       ,Pt_max=0.35        # GeV/c
                       ,OnBeamFV=None,OffBeamFV=None
                       ,cuts_order=['no cut','PIDa','flash','length','non-collinearity','vertex activity','delta phi','soft Pt']
                       ):
    reducedOnBeam, reducedOffBeam = dict(), dict()
    
    reducedOnBeam['no cut'] = OnBeamFV
    reducedOffBeam['no cut'] = OffBeamFV
    
    for i_cut,cut in zip(range(1,len(cuts_order)),cuts_order[1:]):#{
        print 'grabbing reduced data samples after (',cuts_order[i_cut-1],') and applying cut on (',cuts_order[i_cut],')'
        OnBeam_previous_cut = reducedOnBeam[cuts_order[i_cut-1]]
        OffBeam_previous_cut = reducedOffBeam[cuts_order[i_cut-1]]
        
        for sam,sam_name in zip([OnBeam_previous_cut,OffBeam_previous_cut]
                                ,['OnBeam','OffBeam']):#{
            print 'len('+sam_name+'):',len(sam)
            if cut == 'PIDa':
                sam = sam[sam['PIDa_assigned_proton']>PIDa_p_min]
            
            elif cut == 'flash':
                sam = sam[(sam['Nflashes']>0)
                          &(sam['ClosestFlash_TotalPE'] > minPEcut)
                          &(sam['ClosestFlash_YZdistance'] < maxdYZcut)]
            
            elif cut == 'length':
                sam = sam[sam['PIDa_long'] < sam['PIDa_short']]
    
            elif cut == 'non-collinearity':
                sam = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
                    
            elif cut == 'vertex activity':
                R_str = 'RdQaroundVertex'
                box_str='[%d wires x %d ticks]'%(opt_box[0],opt_box[1])
                Ru,Rv,Ry = R_str+'[plane 0]'+box_str,R_str+'[plane 1]'+box_str,R_str+'[plane 2]'+box_str
                sam = sam[(sam[Ru]==1) | (sam[Rv]==1) | (sam[Ry]==1)
                          |
                          (sqrt( np.square(sam[Ru]-1) + square(sam[Rv]-1) + square(sam[Ry]-1) )
                           <= r_max_RdQ_CC1p0pi) ]
            elif cut == 'delta phi':
                sam = sam[np.abs(sam['delta_phi']-180.)<delta_Delta_phi]
                                      
            elif cut == 'soft Pt':
                sam = sam[sam['reco_Pt']<Pt_max]
                                              
            if sam_name=='OnBeam': reducedOnBeam[cut] = sam
            if sam_name=='OffBeam': reducedOffBeam[cut] = sam
        #}

    #}
    return reducedOnBeam,reducedOffBeam
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -



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
                , do_add_legend=True , legend_loc='best',y_label='counts'
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
    
    set_axes(ax,x_label=x_label,y_label=y_label,do_add_grid=True,fontsize=fontsize
             ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
            )
    if do_add_legend: 
        if legend_loc=='bbox':
            leg=plt.legend(bbox_to_anchor=(1.,1.05),fontsize=fontsize,loc=2)
        else:
            leg=plt.legend(fontsize=fontsize,loc=legend_loc)
    plt.tight_layout()
    if do_add_legend is False: return ax
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
# written Nov-9,2017  (last edit Mar-4, 2018)
def OnBeam_minus_OffBeam_1d( OnBeamSample=None , OffBeamSample=None , debug=0
                            , var='PIDa_assigned_proton' , x_label='$PID_a^p$' , y_label='counts'
                            , bins=np.linspace(0,30,31) 
                            , ax=None, figsize=(14,6),fontsize=25
                            , color=OnBeamColor
                            , do_add_MCoverlay=True , MCsamples=None, MCbnbDATAcosmicSamples=None
                            , MC_scaling=MC_scaling_DATAcosmic
                            , do_add_legend=True , legend_loc='best', MCalpha=0.7
                            , do_add_chi2_MC_data=False , chi2_xrange=None, chi2_xy=(0,0), doOffBeam_scaling=True
                            , OriginalOnBeamSample=None , OriginalOffBeamSample=None):
    bin_width = bins[1]-bins[0]
    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    h_OnBeam,edges = np.histogram( OnBeamSample[var] , bins=bins )
    h_OnBeam_err = np.sqrt(h_OnBeam)
    h_OffBeam,edges = np.histogram( OffBeamSample[var] , bins=bins )
    h_OffBeam_err = np.sqrt(h_OffBeam)
    
    if doOffBeam_scaling==True:#{
        h_OnBeam_minus_OffBeam = h_OnBeam - OffBeam_scaling*h_OffBeam
        h_OnBeam_minus_OffBeam_err = np.sqrt( np.abs(h_OnBeam + OffBeam_scaling*OffBeam_scaling*h_OffBeam) )
        Integral = len(OnBeamSample) - OffBeam_scaling*len(OffBeamSample)
        Integral_Original = len(OriginalOnBeamSample) - OffBeam_scaling*len(OriginalOffBeamSample)
    #}
    else:#{
        h_OnBeam_minus_OffBeam = h_OnBeam - h_OffBeam
        h_OnBeam_minus_OffBeam_err = np.sqrt( np.abs(h_OnBeam + h_OffBeam) )
        Integral = len(OnBeamSample) - len(OffBeamSample)
        Integral_Original = len(OriginalOnBeamSample) - len(OriginalOffBeamSample)
    #}


    plt.errorbar( x=0.5*(edges[1:]+edges[:-1])-bin_width, xerr=bin_width/2., markersize=12
                 , y=h_OnBeam_minus_OffBeam , yerr=h_OnBeam_minus_OffBeam_err
                 , fmt='o', color=color , ecolor='black'
                 , label=r'(On-Off) Beam ($\int=$%.1f=%.1f'%(Integral,100*Integral/Integral_Original)+'%)'
                )
    if debug>1: print "OnBeam-OffBeam (bins[:-1]):\n",bins[:-1]
    plt.plot([np.min(ax.get_xlim()),np.min(ax.get_xlim())],[0,0],'--',color='black',linewidth=2)
    
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
    set_axes(ax,x_label=x_label,y_label=y_label,do_add_grid=True,fontsize=fontsize
             ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
            #              ,ylim=(np.min([0,np.min(h_OnBeam_minus_OffBeam-1.1*h_OnBeam_minus_OffBeam_err)])#                     ,np.max(h_OnBeam_minus_OffBeam+1.1*h_OnBeam_minus_OffBeam_err))
            )
    plt.tight_layout()
    if do_add_legend==False: return ax
    return ax,leg
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -


# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Nov-20,2017 (last editted Dec-6)
def plot_stacked_MCsamples( ax=None, debug=0
                           , MCsamples=None
                           , MC_scaling=MC_scaling_DATAcosmic
                           , MCbnbDATAcosmicSamples=None
                           , var=None, x_label='',y_label='', bins=None , alpha=0.8, fontsize=25):
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
    set_axes(ax,x_label=x_label,y_label=y_label,do_add_grid=True,fontsize=fontsize
                 ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
                 )
    
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


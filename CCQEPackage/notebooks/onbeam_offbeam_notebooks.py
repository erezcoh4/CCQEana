import sys; sys.path.insert(0, '../../'); sys.path.insert(0,'../mupClassification/')
from ccqe_notebook_tools import *
import matplotlib.patches as patches
import pprint; pp = pprint.PrettyPrinter(indent=4)
from matplotlib.pylab import *


# -- - --- - --- - --- - --- - --- - --- - --- - --- - --- - --- - -
versions = dict({'Overlay':'prodgenie_bnb_nu_uboone_overlay_mcc8.11_reco2'
                ,'OnBeam':'prod_reco_optfilter_bnb_v12_unblind_mcc8'
                ,'OffBeam':'prod_reco_optfilter_extbnb_v12_mcc8_dev'
                ,'date':'2018_05_11'
                ,'overlay date':'2018_05_11'
                ,'data date':'2018_05_11'
                })
print 'versions:'
pp.pprint(versions)
# -- - --- - --- - --- - --- - --- - --- - --- - --- - --- - --- - -


debug = 0
Nevents=dict()
OffBeam_scaling = float(10947004)/15499027 # from note "Off Beam scaling factor (Ariana, Feb 2018)"
summary = pd.read_csv('/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/summary/'+versions['date']+'/'+versions['Overlay']+'_'+versions['date']+'_summary.csv')

Nevents['OnBeam POT'] = 4.908e+19
Nevents['overlay']      = np.sum(summary.Nevents)
Nevents['overlay POT']  = np.sum(summary.POT)
Nevents['f(POT)']       = Nevents['OnBeam POT']/Nevents['overlay POT']

OnBeamColor,OffBeamColor = 'teal','orange'

print "OffBeam_scaling:",OffBeam_scaling,"= N(on beam)/N(off beam) before SwT"
print "f(POT):",Nevents['f(POT)'],"= N(POT on beam)/N(POT MC)"




# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# written April-29, 2018
def OnBeam_OffBeam( OnBeamSample=None , OffBeamSample=None , debug=0
                   , var='PIDa_assigned_proton' , x_label=r'$PID_a^p$' , y_label='counts'
                   , bins=np.linspace(0,30,31)
                   , ax=None, figsize=(14,6),fontsize=25
                   , color=OnBeamColor
                   , doOffBeam_scaling=True
                   , OriginalOnBeamSample=None , OriginalOffBeamSample=None
                   , remove_ticks_x=False, remove_ticks_y=False):
    bin_width = bins[1]-bins[0]
    mid = 0.5*(bins[:-1]+bins[1:])
    
    
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
    
    
    plt.errorbar( x=mid, xerr=bin_width/2., markersize=12
                 , y=h_OnBeam_minus_OffBeam , yerr=h_OnBeam_minus_OffBeam_err
                 , fmt='o', color=color , ecolor='black'
                 , label=r'(On-Off) Beam ($\int=$%.1f=%.1f'%(Integral,100*Integral/Integral_Original)+'%)'
                 )
    if debug>1: print "OnBeam-OffBeam (bins[:-1]):\n",bins[:-1]
    plt.plot([np.min(ax.get_xlim()),np.min(ax.get_xlim())],[0,0],'--',color='black',linewidth=2)
    
    set_axes(ax,x_label=x_label,y_label=y_label,do_add_grid=True,fontsize=fontsize
             ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
             ,remove_ticks_x=remove_ticks_x, remove_ticks_y=remove_ticks_y
             )
    plt.tight_layout()
    return h_OnBeam_minus_OffBeam,h_OnBeam_minus_OffBeam_err
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -


# ------------------------------------------------
# April-29, 2018
def plot_OnOff_nocut_finalcut(var='theta_12',x_label= r'$\theta_{1,2}$ [deg.]',bins=linspace(0,180,31)
                              ,scaling_name='only POT'
                              ,color ='black',x_ticks=None
                              ,figsize=(24,8)
                              ,overlay_scalings=None,do_OffBeam=False
                              ,reducedOffBeam=None,reducedOnBeam=None,reducedOverlay=None
                              ,f_Cosmic=None,f_CosmicName=r'$f_{cosmic}$'
                              ,chi2_xrange=None,xlim=None
                              ,last_cut_name='soft Pt',last_cut_label='detection + kinematical cuts'
                              ,debug=0
                              ,do_show_cut=True,x_varcut=(0,np.inf)
                              ,do_only_preselection=False
                              ):
    fig = plt.figure(figsize=figsize)
    
    for i_cut,(cut_name,cut_label) in enumerate(zip(['no cut',last_cut_name]
                                                    ,['preselection',last_cut_label])):#{
        
        ax = fig.add_subplot(1,1 if do_only_preselection else 2 ,i_cut + 1)
        
        h_OnOff,h_OnOff_err = OnBeam_OffBeam( OnBeamSample=reducedOnBeam[cut_name] , OffBeamSample=reducedOffBeam[cut_name] , debug=debug
                                             , var=var , x_label=x_label , bins=bins
                                             , ax=ax, color=color
                                             , doOffBeam_scaling=True
                                             , OriginalOnBeamSample=reducedOnBeam['no cut'] , OriginalOffBeamSample=reducedOffBeam['no cut']
                                             )
            
        h_stack,_=plot_stacked_MCsamples(reducedOverlay = reducedOverlay
                                                                              , debug=0
                                                                              , overlay_scaling=overlay_scalings[scaling_name]
                                                                              , cut_name=cut_name
                                                                              , var=var, x_label=x_label, y_label='overlay prediction',xlim=xlim
                                                                              , bins=bins , alpha=0.8, ax=ax
                                                                              , do_add_legend=False
                                                                              );

        chi2 , ndf = chi2_two_histograms( bins=bins, chi2_xrange=(np.min(bins),np.max(bins))
                                                                              , h1=h_OnOff , h2=h_stack
                                                                              , h1err=h_OnOff_err, h2err=np.sqrt(h_stack)
                                                                              , debug=0 )
        chi2_txt = r'$\chi^2/ndf=%.1f/%d$'%(chi2,ndf)
        ax.set_title(cut_label+(', '+f_CosmicName+'=%.2f'%f_Cosmic if f_Cosmic is not None else '') + ',' + chi2_txt
                                                          ,y=1.02,fontsize=25)

        if do_show_cut:
            plt.plot([x_varcut[0],x_varcut[0]],ax.get_ylim(),'--',[x_varcut[1],x_varcut[1]],ax.get_ylim(),'--',color='black')
        if debug: print cut_label,': sum of h_OnOff:',np.sum(h_OnOff),',sum of h_stack:',np.sum(h_stack)
    
        if do_only_preselection and i_cut==0: return

    #}
    plt.tight_layout()
    return h_OnOff,h_stack
# ------------------------------------------------


# ------------------------------------------------
# April-26
# for overlay_vs_BeamOn-Off
def plot_nocut_finalcut(var='theta_12',x_label= r'$\theta_{1,2}$ [deg.]',bins=linspace(0,180,31)
                        ,scaling_name='N(On) scaling and f(cosmic)=0.92'
                        ,scaling_color ='black',x_ticks=None
                        ,figsize=(24,8)
                        ,overlay_scalings=None,do_OffBeam=False
                        ,reducedOffBeam=None,reducedOnBeam=None,reducedOverlay=None
                        ,f_Cosmic=None
                        ,chi2_xrange=None,xlim=None
                        ,last_cut_name='Pt & delta phi',last_cut_label='detection + kinematical cuts'
                        ,debug=0
                        ,do_show_cut=True,x_varcut=(0,np.inf)
                        ,do_only_preselection=False
                        ):
    fig = plt.figure(figsize=figsize)
    
    for i_cut,(cut_name,cut_label) in enumerate(zip(['no cut',last_cut_name]
                                                  ,['preselection',last_cut_label])):#{
                                                  
        ax = fig.add_subplot(1,1 if do_only_preselection else 2 ,i_cut + 1)
        if do_OffBeam:#{
            _,h_OffBeam = plot_OffBeam(OffBeamSample=reducedOffBeam[cut_name],OffBeamFV=reducedOffBeam['no cut']
                                       ,color='white', var=var, bins=bins, ax=ax , do_OffBeam_scaling=True)
        #}
        _,h_OnBeam = plot_OnBeam(OnBeamSample=reducedOnBeam[cut_name],OnBeamFV=reducedOnBeam['no cut']
                                 ,color='black', var=var, bins=bins, ax=ax)
                                 
        h_stack,_=plot_stacked_MCsamples(reducedOverlay = reducedOverlay
                                         , debug=0
                                         , overlay_scaling=overlay_scalings[scaling_name]
                                         , cut_name=cut_name
                                         , var=var, x_label=x_label, y_label='overlay prediction',xlim=xlim
                                         , bins=bins , alpha=0.8, ax=ax
                                         , do_add_legend=False
                                 );
        chi2 , ndf = chi2_two_histograms( bins=bins, chi2_xrange=(np.min(bins),np.max(bins))
                                     , h1=h_OnBeam , h2=h_stack
                                     , h1err=np.sqrt(h_OnBeam), h2err=np.sqrt(h_stack)
                                     , debug=0 )
        chi2_txt = r'$\chi^2/ndf=%.1f/%d$'%(chi2,ndf)
        ax.set_title(cut_label+(', $f_{cosmic}=%.2f$'%f_Cosmic if f_Cosmic is not None else '') + ',' + chi2_txt
                                                  ,y=1.02,fontsize=25)
        if do_show_cut:
            plt.plot([x_varcut[0],x_varcut[0]],ax.get_ylim(),'--',[x_varcut[1],x_varcut[1]],ax.get_ylim(),'--',color='black')
        if debug: print cut_label,': sum of h_OnBeam:',np.sum(h_OnBeam),',sum of h_stack:',np.sum(h_stack)
        if do_only_preselection and i_cut==0: return
    #}
    plt.tight_layout()
    return h_OnBeam,h_stack
# ------------------------------------------------



# ------------------------------------------------
# April-26
def get_samples_scaling( N_total=1 # total integral of all overlay
                        , f_Cosmic=None  # fraction of cosmic in the overlay, need to be in the range 0-1
                        , f_OverlayCosmic=None # by which factor should we increase the cosmic part in the overlay
                        , OverlaySubsamples=None
                        , debug=0
                        , f_POT=Nevents['f(POT)']):
    # return the total scaling factor of each of the overlay subsamples
    scaling,N,f = dict(),dict(),dict()
    
    for pair_type in pair_types: N[pair_type] = float(len(OverlaySubsamples[pair_type]))
    N['MC'] = N['1mu-1p']+N['other pairs']
    
    # option 1:
    # scale the total integral of the overlay to the number of POT in the BeamOn data
    if N_total=='MC(tot) x POT':#{
        N['total']=(N['MC']+N['cosmic'])*f_POT
    #}
    else:#{
        N['total']=N_total
    #}
    
    # option 2:
    # change the fraction of cosmic pairs in the overlay,
    # simultaneously to changing the beam part,
    # so that the total integral is retained
    if f_Cosmic is None:#{
        f['Cosmic'] = N['cosmic']/(N['MC']+N['cosmic'])
    #}
    else:#{
        f['Cosmic'] = f_Cosmic
    #}
    f['MC'] = 1-f['Cosmic']
    
    f['1mu-1p in MC'] = N['1mu-1p']/N['MC']
    f['CC 1p 0pi in MC'] = N['CC 1p 0pi']/N['MC']
    f['other pairs in MC'] = N['other pairs']/N['MC']
    
    # f_Cosmic = fraction of cosmic in the overlay
    scaling['cosmic']     = N['total']*f['Cosmic']/N['cosmic']
    scaling['other pairs']= N['total']*f['MC']*f['other pairs in MC']/N['other pairs']
    scaling['1mu-1p']     = N['total']*f['MC']*f['1mu-1p in MC']/N['1mu-1p']
    scaling['CC 1p 0pi']  = N['total']*f['MC']*f['CC 1p 0pi in MC']/N['CC 1p 0pi']
    
    # option 3:
    # increase the number of cosmic pairs in the overlay,
    # without changing the beam part
    if f_OverlayCosmic is not None: scaling['cosmic'] = f['MC']*f_OverlayCosmic
    
    
    if debug:
        print 'N';pp.pprint(N)
        print 'f';pp.pprint(f)
        print 'scaling:';pp.pprint(scaling)
    return scaling
# ------------------------------------------------


## ------------------------------------------------
## April-4 (last edit April-26)
#def gen_Noverlay(reducedSamples=None,cut_name=''
#                 ,f_POT = Nevents['f(POT)']
#                 ,N_On=1 # number of pairs in BeamOn before event-selection cuts
#                 ,debug=0
#                ):
#    # @return the number of events in each subsample of the overlay, POT-normalized
#    N = dict()
#    N['Cosmic'],N['mup'],N['others'] = float(len(reducedSamples[cut_name]['cosmic'])),float(len(reducedSamples[cut_name]['1mu-1p'])),float(len(reducedSamples[cut_name]['other pairs']))
#    N['MC'] = N['mup'] + N['others']
#    N['Overlay'] = N['Cosmic'] + N['MC']                     
#    N['eff Overlay'] = N['Overlay']/(len(reducedSamples['no cut']['cosmic'])
#                                     +len(reducedSamples['no cut']['1mu-1p'])                                     
#                                     +len(reducedSamples['no cut']['other pairs']))
#    N['Overlay POT Scaled'] = f_POT*(N['Cosmic'] + N['MC'])
#    # scale the cosmic in the MC
#    N['Cosmic original'] = float(len(reducedSamples['no cut']['cosmic']))
#    f_Cosmic = (1./N['Cosmic original'])*(N_On/f_POT - N['MC'])
#    N['Cosmic Scaled'] = f_Cosmic*N['Cosmic']
#    N['Overlay Cosmic Scaled'] = N['MC'] + N['Cosmic Scaled']
#    N['eff Overlay Cosmic Scaled'] = (N['Overlay Cosmic Scaled']
#                                      /(f_Cosmic*len(reducedSamples['no cut']['cosmic'])                                                    
#                                        +len(reducedSamples['no cut']['1mu-1p'])
#                                        +len(reducedSamples['no cut']['other pairs'])))
#    N['Overlay Cosmic & POT Scaled'] = f_POT*N['Overlay Cosmic Scaled']
#    if debug:   print N
#    return N,f_Cosmic
## ------------------------------------------------


# ------------------------------------------------
# March-6, 2018 (last edit May-15)
def apply_cuts_to_data(OnBeamFV=None,OffBeamFV=None
                       ,PIDa_p_min=13,do_PIDaCali=True
                       ,minPEcut = 100
                       ,maxdYZcut = 200
                       ,delta_theta_12=55  # deg.
                       ,opt_box=(50,100) # [Nwires x Nticks]
                       ,r_max_RdQ_CC1p0pi = 0.35 # sphere in U,V,Y space, apply a cut only to CC1p0pi
                       ,delta_Delta_phi=35 # deg.
                       ,Pt_max=0.35        # GeV/c
                       ,cuts_order=['no cut']
                       ,debug=0
                       # replace the cut on PIDa to a cut on chi2_proton
                       ,Chi2Proton_muCandidate_min=80
                       ,Chi2Proton_muCandidate_max=np.inf
                       ,Chi2Proton_pCandidate_min=0
                       ,Chi2Proton_pCandidate_max=30
                       ):
    reducedOnBeam = dict()
    numbers = pd.DataFrame()
    
    reducedOnBeam['no cut'] = OnBeamFV
    if OffBeamFV is not None:
        reducedOffBeam = dict()
        reducedOffBeam['no cut'] = OffBeamFV
    numbers = numbers.append(pd.DataFrame({'$N_{On}$':len(OnBeamFV),r'${\epsilon}_{On}$ [%]':100.
                                           ,'$N_{Off}$':len(OffBeamFV) if OffBeamFV is not None else -1
                                          ,r'${\epsilon}_{Off}$ [%]':100.
                                           ,'$N_{On-Off}$':(len(OnBeamFV)-OffBeam_scaling*len(OffBeamFV)) if OffBeamFV is not None else -1
                                          ,r'${\epsilon}_{On-Off}$ [%]':100
                                           ,'$N_{Off}^{scaled}$':OffBeam_scaling*len(OffBeamFV) if OffBeamFV is not None else -1
                                         },index=['preselection']))

    for i_cut,cut in zip(range(1,len(cuts_order)),cuts_order[1:]):#{
        if debug: print 'grabbing reduced data samples after (',cuts_order[i_cut-1],') and applying cut on (',cuts_order[i_cut],')'
        OnBeam_previous_cut = reducedOnBeam[cuts_order[i_cut-1]]
        if reducedOffBeam is not None: OffBeam_previous_cut = reducedOffBeam[cuts_order[i_cut-1]]
        else: OffBeam_previous_cut = None
        
        for sam,sam_name in zip([OnBeam_previous_cut,OffBeam_previous_cut] if OffBeamFV is not None else [OnBeam_previous_cut]
                                ,['OnBeam','OffBeam']):#{
            if debug: print 'len('+sam_name+'):',len(sam)
            if sam is None: continue
            
            if cut == 'PIDa':
                if do_PIDaCali:
                    sam = sam[sam['pidcali_PIDaYplane_pCandidate']>PIDa_p_min]
                else:
                    sam = sam[sam['pid_PIDaYplane_pCandidate']>PIDa_p_min]
            
            # replace the cut on PIDa to a cut on chi2_proton
            elif cut == 'Chi2Proton':
                sam = sam[ (sam['pidcali_Chi2ProtonYplane_muCandidate']>Chi2Proton_muCandidate_min)
                                         &(sam['pidcali_Chi2ProtonYplane_muCandidate']<Chi2Proton_muCandidate_max)
                                         &(sam['pidcali_Chi2ProtonYplane_pCandidate']>Chi2Proton_pCandidate_min)
                                         &(sam['pidcali_Chi2ProtonYplane_pCandidate']<Chi2Proton_pCandidate_max)]
            
            
            elif cut == 'ClosestFlash':
                sam = sam[(sam['Nflashes']>0)
                                         &(sam['ClosestFlash_TotalPE'] > minPEcut)
                                         &(sam['ClosestFlash_YZdistance'] < maxdYZcut)]
            
            # replace the cut on ClosestFlash to a cut on MatchedFlash
            elif cut == 'MatchedFlash':
                sam = sam[(sam['Nflashes']>0)
                                         &(sam['MatchedFlash_TotalPE'] > minPEcut)
                                         &(sam['MatchedFlash_YZdistance'] < maxdYZcut)]
                    
            elif cut == 'length':
                sam = sam[sam['l_muCandidate'] > sam['l_pCandidate']]
    
            elif cut == 'non-collinearity':
                sam = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
                    
            elif cut == 'vertex activity':
                R_str,box_str = 'RdQaroundVertex','[%d wires x %d ticks]'%(opt_box[0],opt_box[1])
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
        numbers = numbers.append(pd.DataFrame({'$N_{On}$':len(reducedOnBeam[cut]),r'${\epsilon}_{On}$ [%]':(100.*float(len(reducedOnBeam[cut]))/len(OnBeamFV))
                                               ,'$N_{Off}$':len(reducedOffBeam[cut]) if OffBeamFV is not None else -1
                                              ,r'${\epsilon}_{Off}$ [%]':(100.*float(len(reducedOffBeam[cut]))/len(OffBeamFV)) if OffBeamFV is not None else -1
                                              ,'$N_{On-Off}$':(len(reducedOnBeam[cut])-OffBeam_scaling*len(reducedOffBeam[cut])) if OffBeamFV is not None else -1
                                              ,r'${\epsilon}_{On-Off}$ [%]':(100*(len(reducedOnBeam[cut])-OffBeam_scaling*len(reducedOffBeam[cut]))                                                           
                                                                      /(len(OnBeamFV)-OffBeam_scaling*len(OffBeamFV))) if OffBeamFV is not None else -1
                                               ,'$N_{Off}^{scaled}$':OffBeam_scaling*len(reducedOffBeam[cut]) if OffBeamFV is not None else -1
                                              },index=[cut]))
    #}
    if OffBeamFV is not None: reducedOnBeam,reducedOffBeam,numbers
    return reducedOnBeam,numbers
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -

#
## ------------------------------------------------
## March-6, 2018 (last edit May 14)
#def apply_cuts_to_overlay(OverlaySamples=None
#                          ,N_On=1 # number of pairs in BeamOn before event-selection cuts
#                          ,PIDa_p_min=13
#                          ,do_PIDaCali=True
#                          ,minPEcut = 100
#                          ,maxdYZcut = 200
#                          ,delta_theta_12=55  # deg.
#                          ,opt_box=(50,100) # [Nwires x Nticks]
#                          ,r_max_RdQ_CC1p0pi = 0.35 # sphere in U,V,Y space, apply a cut only to CC1p0pi
#                          ,delta_Delta_phi=35 # deg.
#                          ,Pt_max=0.35        # GeV/c
#                          ,cuts_order=['no cut']
#                          ,debug=0
#                          ,f_POT=Nevents['f(POT)']
#                          # replace the cut on PIDa to a cut on chi2_proton
#                          ,Chi2Proton_muCandidate_min=80
#                          ,Chi2Proton_muCandidate_max=np.inf
#                          ,Chi2Proton_pCandidate_min=0
#                          ,Chi2Proton_pCandidate_max=30
#                          ):
#    reducedSamples = dict()
#    pureffOverlay,numbers = pd.DataFrame(),pd.DataFrame()
#
#
#    cut_name = 'no cut'
#    reducedSamples[cut_name] = dict()
#    for pair_type in pair_types: reducedSamples['no cut'][pair_type] = OverlaySamples[pair_type]
#    pureffOverlay = get_pureff_cut(OverlaySamples=OverlaySamples
#                                   ,pureff=pureffOverlay,cut_name='no cut',reduced=reducedSamples['no cut'])
#    Noverlay,_ = gen_Noverlay(reducedSamples=reducedSamples,cut_name=cut_name,N_On=N_On,debug=debug,f_POT=f_POT )
#    numbers = numbers.append(pd.DataFrame({ r'$N_{Overlay}$':Noverlay['Overlay']
#                                          ,r'${\epsilon}_{Overlay}$ [%]':100
#                                          ,r'$N_{Overlay, cosmic-scaled}$':Noverlay['Overlay Cosmic Scaled']
#                                          ,r'$N_{Overlay, cosmic-scaled}^{POT-scaled}$':Noverlay['Overlay Cosmic & POT Scaled']
#                                          ,r'${\epsilon}_{Overlay, cosmic-scaled}$ [%]':100
#                                          ,r'$N_{cosmic}$':Noverlay['Cosmic']
#                                          ,r'$N_{cosmic, cosmic-scaled}$':Noverlay['Cosmic Scaled']
#                                          ,r'$N_{cosmic, cosmic-scaled}^{POT-scaled}$':f_POT*Noverlay['Cosmic Scaled']
#                                          },index=['preselection']))
#
#
#    for i_cut,cut in zip(range(1,len(cuts_order)),cuts_order[1:]):#{
#        reduced = dict()
#        
#        if debug: print 'grabbing reduced samples after (',cuts_order[i_cut-1],') and applying cut on (',cuts_order[i_cut],')'
#        samples_previous_cut = reducedSamples[cuts_order[i_cut-1]]
#        if debug: print 'len(samples_previous_cut):',len(samples_previous_cut)
#        
#        for pair_type in pair_types:#{
#            sam = samples_previous_cut[pair_type]
#            
#            if cut == 'PIDa':
#                if do_PIDaCali:
#                    reduced[pair_type] = sam[sam['pidcali_PIDaYplane_pCandidate']>PIDa_p_min]
#                else:
#                    reduced[pair_type] = sam[sam['pid_PIDaYplane_pCandidate']>PIDa_p_min]
#        
#            # replace the cut on PIDa to a cut on chi2_proton
#            elif cut == 'Chi2Proton':
#                reduced[pair_type] = sam[ (sam['pidcali_Chi2ProtonYplane_muCandidate']>Chi2Proton_muCandidate_min)
#                                         &(sam['pidcali_Chi2ProtonYplane_muCandidate']<Chi2Proton_muCandidate_max)
#                                         &(sam['pidcali_Chi2ProtonYplane_pCandidate']>Chi2Proton_pCandidate_min)
#                                         &(sam['pidcali_Chi2ProtonYplane_pCandidate']<Chi2Proton_pCandidate_max)]
#            
#            
#            elif cut == 'ClosestFlash':
#                reduced[pair_type] = sam[(sam['Nflashes']>0)
#                                         &(sam['ClosestFlash_TotalPE'] > minPEcut)
#                                         &(sam['ClosestFlash_YZdistance'] < maxdYZcut)]
#            
#            # replace the cut on ClosestFlash to a cut on MatchedFlash
#            elif cut == 'MatchedFlash':
#                reduced[pair_type] = sam[(sam['Nflashes']>0)
#                                         &(sam['MatchedFlash_TotalPE'] > minPEcut)
#                                         &(sam['MatchedFlash_YZdistance'] < maxdYZcut)]
#         
#         
#
#
#            elif cut == 'length':
#                reduced[pair_type] = sam[sam['l_muCandidate'] > sam['l_pCandidate']]
#    
#            elif cut == 'non-collinearity':
#                reduced[pair_type] = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
#            
#            elif cut == 'vertex activity':
#                R_str,box_str = 'RdQaroundVertex','[%d wires x %d ticks]'%(opt_box[0],opt_box[1])
#                Ru,Rv,Ry = R_str+'[plane 0]'+box_str,R_str+'[plane 1]'+box_str,R_str+'[plane 2]'+box_str
#                reduced[pair_type] = sam[(sam[Ru]==1) | (sam[Rv]==1) | (sam[Ry]==1)
#                                         |
#                                         (sqrt( np.square(sam[Ru]-1) + square(sam[Rv]-1) + square(sam[Ry]-1) )
#                                          <= r_max_RdQ_CC1p0pi) ]
#            elif cut == 'delta phi':
#                reduced[pair_type] = sam[np.abs(sam['delta_phi']-180.)<delta_Delta_phi]
#
#            elif cut == 'soft Pt':
#                reduced[pair_type] = sam[sam['reco_Pt']<Pt_max]
#        #}
#        reducedSamples[cut] = reduced
#        pureffOverlay = get_pureff_cut(OverlaySamples=OverlaySamples
#                                       ,pureff=pureffOverlay,cut_name=cut,reduced=reduced)
#        Noverlay,_ = gen_Noverlay(reducedSamples=reducedSamples,cut_name=cut,N_On=N_On,debug=debug,f_POT=f_POT )
#        numbers = numbers.append(pd.DataFrame({ r'$N_{Overlay}$':Noverlay['Overlay']
#                                          ,r'${\epsilon}_{Overlay}$ [%]':100.*Noverlay['eff Overlay']
#                                          ,r'$N_{Overlay, cosmic-scaled}$':Noverlay['Overlay Cosmic Scaled']
#                                          ,r'$N_{Overlay, cosmic-scaled}^{POT-scaled}$':Noverlay['Overlay Cosmic & POT Scaled']
#                                          ,r'${\epsilon}_{Overlay, cosmic-scaled}$ [%]':100.*Noverlay['eff Overlay Cosmic Scaled']
#                                          ,r'$N_{cosmic}$':Noverlay['Cosmic']
#                                          ,r'$N_{cosmic, cosmic-scaled}$':Noverlay['Cosmic Scaled']
#                                          ,r'$N_{cosmic, cosmic-scaled}^{POT-scaled}$':f_POT*Noverlay['Cosmic Scaled']
#                                          },index=[cut]))
#                                          
#    #}
#    return reducedSamples,pureffOverlay,numbers
## -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
#
#


#
#
## ------------------------------------------------
## March-6, 2018
#def get_Nreduced(OverlaySamples=None,reduced = dict()):
#    Noriginal , Nreduced , freduced = dict() , dict() , dict()
#    for pair_type in pair_types:
#        sam = OverlaySamples[pair_type]
#        Noriginal[pair_type] = len(OverlaySamples[pair_type])
#        Nreduced[pair_type] = float(len(reduced[pair_type]))
#        freduced[pair_type] = 100.0 * Nreduced[pair_type]/Noriginal[pair_type]
#    return Nreduced , freduced
## ------------------------------------------------
#
## ------------------------------------------------
## March 5, 2018
#def get_pureff_cut(OverlaySamples=None,reduced=None,pureff=None, cut_name = 'PIDa'):
#    eff,pur = dict(),dict()
#    Nreduced , freduced = get_Nreduced(OverlaySamples=OverlaySamples,reduced=reduced)
#    Ntot = (Nreduced['1mu-1p']+Nreduced['cosmic']+Nreduced['other pairs'])
#    eff['1mu-1p'] = freduced['1mu-1p']
#    pur['1mu-1p'] = 100.*Nreduced['1mu-1p']/Ntot if Ntot>0 else 0
#    eff['CC 1p 0pi'] = freduced['CC 1p 0pi']
#    pur['CC 1p 0pi'] = 100.*Nreduced['CC 1p 0pi']/Ntot if Ntot>0 else 0
#    pureff_cut = pd.DataFrame({'label':cut_name
#                              ,'$\mu p$ eff.':'%.1f'%eff['1mu-1p']+'%'
#                              ,'$\mu p$ pur.':'%.1f'%pur['1mu-1p']+'%'
#                              ,'CC$0\pi 1 p$ eff.':'%.1f'%freduced['CC 1p 0pi']+'%'
#                              ,'CC$0\pi 1 p$ pur.':'%.1f'%(100.*Nreduced['CC 1p 0pi']/Ntot if Ntot>0 else 0)+'%'}
#                              , index=[cut_name]
#                              )
#
#    for pair_type in pair_types: pureff_cut[pair_type] = '%.1f'%freduced[pair_type]+'%' +' (%.0f)'%Nreduced[pair_type]
#    pureff = pureff.append(pureff_cut)
#    return pureff
## ------------------------------------------------


# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Dec-6,2017 (last edit April-7)
def plot_OnBeam(OnBeamSample=None,OnBeamFV=None
                , var='PIDa_assigned_proton' , x_label='$PID_a^p$'                 
                , bins=np.linspace(0,30,31)                 
                , ax=None, figsize=(14,6),fontsize=25                
                , color=OnBeamColor
                , do_add_legend=False , legend_loc='best',y_label='counts'
                , remove_ticks_x=False, remove_ticks_y=False):
    bin_width = bins[1]-bins[0]
    mid = 0.5*(bins[:-1]+bins[1:])

    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    x = OnBeamSample[var]
    h_OnBeam,edges = np.histogram( x , bins=bins )
    h_OnBeam_err = np.sqrt(h_OnBeam)
    
    plt.errorbar( x = mid, xerr=bin_width/2., markersize=12
                 , y=h_OnBeam , yerr=h_OnBeam_err
                 , fmt='o', color=color , ecolor=color
                 , label='BNB (%d=%.1f'%(len(OnBeamSample),100*float(len(OnBeamSample))/len(OnBeamFV))+'%)'
                )
    plt.plot([0,0],[0,0],'--',color='black',linewidth=2)
    
    set_axes(ax,x_label=x_label,y_label='counts',do_add_grid=True,fontsize=fontsize
             ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
             ,remove_ticks_x=remove_ticks_x, remove_ticks_y=remove_ticks_y
            )
    if do_add_legend: 
        if legend_loc=='bbox':
            leg=plt.legend(bbox_to_anchor=(1.,1.05),fontsize=fontsize,loc=2)
        else:
            leg=plt.legend(fontsize=fontsize,loc=legend_loc)
    plt.tight_layout()
    if do_add_legend is False: return ax,h_OnBeam    
    return ax,leg,h_OnBeam
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -



# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Dec-6,2017 (last edit April-7, 2018)
def plot_OffBeam(OffBeamSample=None,OffBeamFV=None
                , var='PIDa_assigned_proton' , x_label='$PID_a^p$'                 
                , bins=np.linspace(0,30,31)                 
                , ax=None, figsize=(14,6),fontsize=25                
                , color=OffBeamColor
                , do_add_legend=False , legend_loc='best',y_label='counts'                 
                , do_OffBeam_scaling=True, remove_ticks_x=False, remove_ticks_y=False):
    bin_width = bins[1]-bins[0]
    mid = 0.5*(bins[:-1]+bins[1:])

    if ax is None: fig,ax=plt.subplots(figsize=figsize)
    x = OffBeamSample[var]
    
    h_OffBeam,edges = np.histogram( x , bins=bins )
    h_OffBeam_err = np.sqrt(h_OffBeam)
    
    h_OffBeam = OffBeam_scaling*h_OffBeam if do_OffBeam_scaling else h_OffBeam
    h_OffBeam_err = OffBeam_scaling*h_OffBeam_err if do_OffBeam_scaling else h_OffBeam_err
    
    plt.errorbar( x = mid, xerr=bin_width/2., markersize=12
                 , y=h_OffBeam
                 , yerr=h_OffBeam_err
                 , fmt='s', color=color , ecolor=color
                 , label='extBNB (%.1f=%.1f'%(OffBeam_scaling*len(OffBeamSample),100*float(len(OffBeamSample))/len(OffBeamFV))+'%)'
                )
    plt.plot([0,0],[0,0],'--',color='black',linewidth=2)
    
    set_axes(ax,x_label=x_label,y_label=y_label,do_add_grid=True,fontsize=fontsize
             ,xlim=(np.min(bins)-bin_width,np.max(bins)+bin_width)
             ,remove_ticks_x=remove_ticks_x, remove_ticks_y=remove_ticks_y
            )
    if do_add_legend: 
        if legend_loc=='bbox':
            leg=plt.legend(bbox_to_anchor=(1.,1.05),fontsize=fontsize,loc=2)
        else:
            leg=plt.legend(fontsize=fontsize,loc=legend_loc)
    plt.tight_layout()
    if do_add_legend is False: return ax,h_OffBeam
    return ax,leg,h_OffBeam
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -







# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Dec-4,2017
def draw_var(cut_name=None,var=None,x_label=None,bins=None,debug=0
             ,OnBeamFV=None,OffBeamFV=None
             ,do_cosmic_only=True,do_bnb_only=True
             ,chi2_xrange=None,chi2_xy=None
             ,reduced_OnBeam=None,reduced_OffBeam=None
             ,OverlaySamples=None, reduced_OverlaySamples=None
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
                                     ,OffBeamFV=OffBeamFV,OverlaySamples=OverlaySamples
                                     ,cosmic_overlay_sample = reduced_OverlaySamples[cut_name]['cosmic']
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
                                   ,MCsamples=reduced_OverlaySamples[cut_name]
                                   ,OverlaySamples=OverlaySamples
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
                            , do_add_cosmic_overlay=True , cosmic_overlay_sample=None, OverlaySamples=None
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
                            , alpha=overlay_alpha,OverlaySamples=OverlaySamples)
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
                        , alpha=0.5, Int_OffBeam=1,OverlaySamples=None):
    x = cosmic_overlay_sample[var]
    bin_width = bins[1]-bins[0]
        
    h_overlay,edges = np.histogram( x , bins=bins )
    Int_overlay = np.sum([h_overlay[i]*bin_width for i in range(len(h_overlay))])
    SUMoverlay = np.sum(h_overlay)    

    cosmic_overlay_scaling = Int_OffBeam/Int_overlay    
    ax.hist( x , weights=cosmic_overlay_scaling*np.ones(len(x))
                     , bins=bins-0.5*bin_width, width=bin_width
                     , color=MCcolors[2]
                     , label='cosmic overlay (%d=%.1f'%(len(cosmic_overlay_sample),100*float(len(cosmic_overlay_sample))/len(OverlaySamples['cosmic']))+'%)'
                     , alpha=alpha)       
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -        





# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# written Nov-9,2017  (last edit April-4, 2018)
def OnBeam_minus_OffBeam_1d( OnBeamSample=None , OffBeamSample=None , debug=0
                            , var='PIDa_assigned_proton' , x_label=r'$PID_a^p$' , y_label='counts'
                            , bins=np.linspace(0,30,31) 
                            , ax=None, figsize=(14,6),fontsize=25
                            , color=OnBeamColor
                            , do_add_MCoverlay=True , MCsamples=None, OverlaySamples=None
                            , MC_scaling=Nevents['f(POT)']
                            , do_add_legend=True , legend_loc='best', MCalpha=0.7
                            , do_add_chi2_MC_data=False , chi2_xrange=None, chi2_xy=(0,0), doOffBeam_scaling=True
                            , OriginalOnBeamSample=None , OriginalOffBeamSample=None
                            , remove_ticks_x=False, remove_ticks_y=False):
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


    plt.errorbar( x=edges[:-1], xerr=bin_width/2., markersize=12
                 , y=h_OnBeam_minus_OffBeam , yerr=h_OnBeam_minus_OffBeam_err
                 , fmt='o', color=color , ecolor='black'
                 , label=r'(On-Off) Beam ($\int=$%.1f=%.1f'%(Integral,100*Integral/Integral_Original)+'%)'
                )
    if debug>1: print "OnBeam-OffBeam (bins[:-1]):\n",bins[:-1]
    plt.plot([np.min(ax.get_xlim()),np.min(ax.get_xlim())],[0,0],'--',color='black',linewidth=2)
    
    if do_add_MCoverlay:        
        h_MC, bins_MC = plot_stacked_MCsamples( ax=ax, debug=debug
                               , MCsamples = MCsamples , OverlaySamples=OverlaySamples
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
             ,remove_ticks_x=remove_ticks_x, remove_ticks_y=remove_ticks_y
            #              ,ylim=(np.min([0,np.min(h_OnBeam_minus_OffBeam-1.1*h_OnBeam_minus_OffBeam_err)])#                     ,np.max(h_OnBeam_minus_OffBeam+1.1*h_OnBeam_minus_OffBeam_err))
            )
    plt.tight_layout()
    if do_add_legend==False: return ax
    return ax,leg
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -










# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# Nov-20,2017 (last editted April-26, 2018)
def plot_stacked_MCsamples( reducedOverlay=None
                           , ax=None, debug=0,overlay_scaling=None,cut_name='no cut'
                           , MC_scaling=Nevents['f(POT)']
                           , var=None, x_label='',y_label='', bins=None , alpha=0.8, fontsize=25
                           , remove_ticks_x=False, remove_ticks_y=False 
                           , xlim=None
                           , do_add_legend=False
                           , do_individual_histograms=True
                           , stackColor='black',stackLabel='overlay'):
    '''
    return: h, bins
            stacked histogram values and bins,
            of the samples from the overlay: 
            (cosmic, other-pairs) + 1mu 1p pairs
    '''
    bin_width = bins[1]-bins[0]
    mid = 0.5*(bins[:-1]+bins[1:])
    h,labels,colors,N = dict(),dict(),dict(),dict()
    
    for i_pair_type,pair_type in enumerate(pair_types):
        sample = reducedOverlay[cut_name][pair_type]
        N[pair_type] = float(len(sample))
        Noriginal = len(reducedOverlay['no cut'][pair_type])
        labels[pair_type] = MClabels[i_pair_type]+' (%.1f=%.1f'%(N[pair_type],100.*N[pair_type]/Noriginal)+'%)'
        colors[pair_type] = MCcolors[i_pair_type];
        x = sample[var]; x = x[x<1e5];
        h[pair_type],edges = np.histogram(x,bins=bins)
        h[pair_type+' scaled'] = overlay_scaling[pair_type]*h[pair_type]
    # -- - - - --------- - - -- ---- -  - --- -- -- -- --
    if do_individual_histograms:#{

        # mu-p
        plt.bar(mid,h['cosmic scaled']+h['other pairs scaled']+h['1mu-1p scaled'] , width=bin_width
            ,color=colors['1mu-1p'],alpha=alpha, label=labels['1mu-1p'])
        # CC 1p 0pi
        plt.bar(mid,h['cosmic scaled']+h['other pairs scaled']+h['CC 1p 0pi scaled'] , width=bin_width
                ,color=colors['CC 1p 0pi'],alpha=alpha, label=labels['CC 1p 0pi'])
        # other pairs
        plt.bar(mid,h['cosmic scaled']+h['other pairs scaled'] , width=bin_width
                        ,color=colors['other pairs'],alpha=alpha , label=labels['other pairs'])
                        # cosmic
        plt.bar(mid, h['cosmic scaled'] , width=bin_width
                        ,color=colors['cosmic'],alpha=alpha, label=labels['cosmic'])
    #}
    # all
    h_stack = h['cosmic scaled']+h['other pairs scaled']+h['1mu-1p scaled']
    plt.step(mid+0.5*bin_width,h_stack ,color=stackColor,alpha=alpha, label=stackLabel)

    if np.max(h_stack)>np.max(ax.get_ylim()): ax.set_ylim(np.min(ax.get_ylim()),1.05*np.max(h_stack))
    set_axes(ax,x_label=x_label,y_label=y_label,do_add_grid=True,fontsize=fontsize          
             ,xlim=(np.min(bins)-0.5*bin_width,np.max(bins)+0.5*bin_width) if xlim is None else xlim
             ,remove_ticks_x=remove_ticks_x             
             ,remove_ticks_y=remove_ticks_y                 
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


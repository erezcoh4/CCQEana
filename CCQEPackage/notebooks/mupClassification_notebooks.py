from ccqe_notebook_tools import * 


vertices_files_path = '/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/ccqe_candidates/'

MCbnbDATAcosmicSamples=dict()
reduced_MCbnbDATAcosmicSamples=dict(dict())
pureff_MCbnbDATAcosmic = pd.DataFrame()
pureff_MCbnbDATAcosmic_numbers = pd.DataFrame()

MCbnbMCcosmicSamples=dict()
reduced_MCbnbMCcosmicSamples=dict(dict())
pureff_MCbnbMCcosmic = pd.DataFrame()
pureff_MCbnbMCcosmic_numbers = pd.DataFrame()


# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# May-29,2018
def Gaussian_resolution(df=None,vgen='',vrec=''
                        ,bins=None
                        ,truncation_limit=np.inf
                        ,p0=[0,1] # initial guess
                        ,do_plot_bestfit=False,precision=3
                        ):
    delta = df[vgen] - df[vrec]
    histo,_ = np.histogram( delta , bins=bins , normed=1);
    trunc_df = df[np.abs(df[vgen] - df[vrec])<truncation_limit]
    trunc_delta = trunc_df[vgen] - trunc_df[vrec]
    h_trun,bins = np.histogram( trunc_delta , bins=bins , normed=1);
    mid = 0.5*(bins[1:]+bins[:-1])
    xdata,ydata = mid,h_trun
    pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), xdata, ydata, p0=p0)
    results = dict({'mid':mid,'histo':histo,'bin_width':(bins[1]-bins[0])
                   ,'mu':pars[0],'mu_err':np.sqrt(cov[0,0])
                   ,'sigma':pars[1],'sigma_err':np.sqrt(cov[1,1])
                   })
    if do_plot_bestfit:
        x = linspace(np.min(bins),np.max(bins),100)
        plt.plot(x, norm.pdf(x,*pars), 'k--',linewidth = 2, label=None)

    return results
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -

# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# May-29,2018
def Gaussian_relative_resolution(df=None,vgen='',vrec=''
                                 ,bins=linspace(-100,100,100)
                                 ,truncation_limit=np.inf
                                 ,p0=[0,1] # initial guess
                                 ,do_plot_bestfit=False,precision=3
                                 ):
    R = 100.*(df[vgen] - df[vrec])/df[vgen]
    histo,_ = np.histogram( R , bins=bins , normed=1);
    trunc_df = df[np.abs(df[vgen] - df[vrec])<truncation_limit]
    trunc_R = 100.*(trunc_df[vgen] - trunc_df[vrec])/trunc_df[vgen]
    h_trun,bins = np.histogram( trunc_R , bins=bins , normed=1);
    mid = 0.5*(bins[1:]+bins[:-1])
    xdata,ydata = mid,h_trun
    pars, cov = curve_fit(lambda x, mu, sig : norm.pdf(x, loc=mu, scale=sig), xdata, ydata, p0=p0)
    results = dict({'mid':mid,'histo':histo,'bin_width':(bins[1]-bins[0])
                   ,'mu':pars[0],'mu_err':np.sqrt(cov[0,0])
                   ,'sigma':pars[1],'sigma_err':np.sqrt(cov[1,1])
                   })

    if do_plot_bestfit:
        x = linspace(np.min(bins),np.max(bins),100)
        plt.plot(x, norm.pdf(x,*pars), 'k--',linewidth = 2, label=None)

    return results
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -



# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -
# May-29, 2018
def plot_purity( OverlaySamples=None
                , ax=None, debug=0,overlay_scaling=None
                , purity_pair_type='CC1p0pi'
                , var=None, bins=None
                , x_label=''
                , y_label='purity [%]'
                , fontsize=25,markersize=15
                , xlim=None
                , color='blue',ecolor='black',label=r'$CC1p0\pi$ purity'
                ):
    '''
        return: x , pur
        purity of the sample 'purity_subsample' in the overlay for the given purity_pair_type
        '''
    bin_width = bins[1]-bins[0]
    mid = 0.5*(bins[:-1]+bins[1:])
    h = dict()
    for i_pair_type,pair_type in enumerate(pair_types):
        h[pair_type],edges = np.histogram(OverlaySamples[pair_type][var],bins=bins)
        h[pair_type+' scaled'] = overlay_scaling[pair_type]*h[pair_type] if overlay_scaling else h[pair_type]
    # -- - - - --------- - - -- ---- -  - --- -- -- -- --
    h_ovrelay = h['cosmic scaled']+h['other pairs scaled']+h['1mu-1p scaled']
    # for the statistical uncertainty of the purity to really decrease with increasing statistics
    # in the overlay,
    # we need to bring back the "real" statistics in the overlay,
    # by multiplying back with overlay_scaling['N(Ovelay)/N(On)']
    numerator = np.array(h[purity_pair_type+' scaled'])*overlay_scaling['N(Ovelay)/N(On)']
    denominator = np.array([h_ovrelay[i] if h_ovrelay[i] else 1000 for i in range(len(h_ovrelay))])*overlay_scaling['N(Ovelay)/N(On)']
    purity = numerator/denominator
    purity_err = np.array([purity[i]*np.sqrt((1./numerator[i] if numerator[i]>0.5 else 0)
                                             + (1./denominator[i] if denominator[i]>0.5 else 0))
                           for i in range(len(purity))])
    if debug:
        print 'h_ovrelay:',h_ovrelay
        print 'h['+purity_pair_type+' scaled]:',h[purity_pair_type+' scaled']
        print 'purity:',purity
        print 'purity_err:',purity_err


    plt.errorbar(x=mid,xerr=0.5*bin_width
             , y=100.*purity
             , yerr=100.*purity_err
             ,color=color, label=label
             ,fmt='o',markersize=markersize,capthick=3,capsize=3,ecolor=ecolor)
    
    set_axes(ax,x_label=x_label,y_label=y_label,do_add_grid=True,alpha_grid=1,fontsize=fontsize
             ,ylim=(0,100)
             ,xlim=(np.min(bins)-0.5*bin_width,np.max(bins)+0.5*bin_width) if xlim is None else xlim
             )
    return mid,purity,purity_err
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -


# ------------------------------------------------
# March-6, 2018 (last edit Aug-12, 2018)
def get_Nreduced(OverlaySamples=None,reduced = dict(),OriginalOverlaySamples=None):
    Noriginal , Nreduced , freduced = dict() , dict() , dict()
    for pair_type in pair_types:
        sam = OverlaySamples[pair_type]
        Noriginal[pair_type] = len(OriginalOverlaySamples[pair_type])
        Nreduced[pair_type] = float(len(reduced[pair_type]))
        freduced[pair_type] = 100.0 * Nreduced[pair_type]/Noriginal[pair_type]
    return Nreduced , freduced
# ------------------------------------------------


# ------------------------------------------------
# March 5, 2018 (last edit Aug-15, 2018)
def get_pureff_cut(OverlaySamples=None,reduced=None,pureff=None, cut_name = 'PIDa',OriginalOverlaySamples=None):
    eff,pur = dict(),dict()
    Nreduced , freduced = get_Nreduced(OverlaySamples=OverlaySamples,reduced=reduced,OriginalOverlaySamples=OriginalOverlaySamples)
    Ntot = (Nreduced['1mu-1p']+Nreduced['cosmic']+Nreduced['other pairs'])
    eff['1mu-1p'] = freduced['1mu-1p']
    pur['1mu-1p'] = 100.*Nreduced['1mu-1p']/Ntot if Ntot>0 else 0
    eff['CC1p0pi'] = freduced['CC1p0pi']
    pur['CC1p0pi'] = 100.*Nreduced['CC1p0pi']/Ntot if Ntot>0 else 0
    eff['CC1p'] = freduced['CC1p']
    pur['CC1p'] = 100.*Nreduced['CC1p']/Ntot if Ntot>0 else 0
    pureff_cut = pd.DataFrame({'label':cut_name
                              ,'$\mu p$ eff.':'%.1f'%eff['1mu-1p']+'%'
                              ,'$\mu p$ pur.':'%.1f'%pur['1mu-1p']+'%'
                              ,'CC$0\pi 1 p$ eff.':'%.1f'%freduced['CC1p0pi']+'%'
                              ,'CC$0\pi 1 p$ pur.':'%.1f'%(100.*Nreduced['CC1p0pi']/Ntot if Ntot>0 else 0)+'%'
                              ,'CC$1 p$ eff.':'%.1f'%freduced['CC1p']+'%'
                              ,'CC$1 p$ pur.':'%.1f'%(100.*Nreduced['CC1p']/Ntot if Ntot>0 else 0)+'%'}
                              , index=[cut_name]
                              )
    for pair_type in pair_types: pureff_cut[pair_type] = '%.1f'%freduced[pair_type]+'%' +' (%.0f)'%Nreduced[pair_type]
    pureff = pureff.append(pureff_cut)
    return pureff
# ------------------------------------------------





# ------------------------------------------------
# March-6, 2018 (last edit Aug-15,2018)
def apply_cuts_to_overlay(OverlaySamples=None
                          ,PIDa_p_min=13,do_PIDaCali=True
                          ,minPEcut = 150
                          ,maxdYZcut = 200
                          ,delta_theta_12=55  # deg.
                          ,opt_box=(50,100) # [Nwires x Nticks]
                          ,r_max_RdQ_CC1p = 0.43 # sphere in U,V,Y space, apply a cut only to CC1p0pi
                          ,delta_Delta_phi=35 # deg.
                          ,Pt_max=0.35        # GeV/c
                          ,theta_pq_max=25    # deg.
                          # replace the cut on PIDa to a cut on chi2_proton
                          ,Chi2Proton_muCandidate_min=80
                          ,Chi2Proton_muCandidate_max=np.inf
                          ,Chi2Proton_pCandidate_min=0
                          ,Chi2Proton_pCandidate_max=30
                          # hyper-parameters
                          ,overlay_scaling=None
                          ,cuts_order=['no cut']
                          ,debug=0
                          ):
    reducedSamples = dict()
    pureffOverlay,numbers = pd.DataFrame(),pd.DataFrame()
    
    cut_name = 'no cut'
    reducedSamples['no cut'] = dict()
    reducedSamples[cut_name] = dict()
    for pair_type in pair_types: reducedSamples['no cut'][pair_type] = OverlaySamples[pair_type]
    pureffOverlay = get_pureff_cut(OverlaySamples=OverlaySamples,OriginalOverlaySamples=reducedSamples['no cut']
                                   ,pureff=pureffOverlay,cut_name=cut_name,reduced=reducedSamples[cut_name])
    Noverlay = get_Noverlay(reducedSamples=reducedSamples,cut_name=cut_name,debug=debug,overlay_scaling=overlay_scaling )
    numbers = numbers.append(pd.DataFrame({ r'$N_{Overlay}$':Noverlay['Overlay']
                                          ,r'${\epsilon}_{Overlay}$ [%]':100
                                          
                                          ,r'$N_{cosmic}$':Noverlay['cosmic']
                                          ,r'$N_{other pairs}$':Noverlay['other pairs']
                                          ,r'$N_{1mu-1p}$':Noverlay['1mu-1p']
                                          ,r'$N_{CC1p0pi}$':Noverlay['CC1p0pi']
                                          ,r'$N_{CC1p}$':Noverlay['CC1p']

                                          ,r'${\epsilon}_{cosmic}$ [%]':100
                                          ,r'${\epsilon}_{CC1p0pi}$ [%]':100
                                          ,r'${\epsilon}_{CC1p}$ [%]':100
                                          ,r'${\epsilon}_{1mu-1p}$ [%]':100
                                          ,r'${\epsilon}_{other pairs}$ [%]':100
                                           
                                          ,r'${\mathcal{p}}_{1mu-1p}$ [%]':100.*Noverlay['pur 1mu-1p']
                                          ,r'${\mathcal{p}}_{CC1p0pi}$ [%]':100.*Noverlay['pur CC1p0pi']
                                          ,r'${\mathcal{p}}_{CC1p}$ [%]':100.*Noverlay['pur CC1p']
                                          
                                          # scaling and re-weighting
                                          ,r'$N_{Overlay scaled}$':Noverlay['Overlay scaled']
                                          ,r'${\epsilon}_{Overlay scaled}$ [%]':100
                                          
                                          ,r'$N_{cosmic scaled}$':Noverlay['cosmic scaled']
                                          ,r'$N_{other pairs scaled}$':Noverlay['other pairs scaled']
                                          ,r'$N_{1mu-1p scaled}$':Noverlay['1mu-1p scaled']
                                          ,r'$N_{CC1p0pi scaled}$':Noverlay['CC1p0pi scaled']
                                          ,r'$N_{CC1p scaled}$':Noverlay['CC1p scaled']
                                           },index=['preselection']))
                                          
    for i_cut,cut in zip(range(1,len(cuts_order)),cuts_order[1:]):#{
        reduced = dict()
        if debug: print 'grab reduced samples after (',cuts_order[i_cut-1],') and apply cut on (',cuts_order[i_cut],')'
        samples_previous_cut = reducedSamples[cuts_order[i_cut-1]]
        
        for pair_type in pair_types:#{
            sam = samples_previous_cut[pair_type]
            if debug: print 'sam('+pair_type+'):',len(sam)
                                    
                                    
            if cut == 'PIDa':#{
                if do_PIDaCali:
                    reduced[pair_type] = sam[sam['pidcali_PIDaYplane_pCandidate']>PIDa_p_min]
                else:
                    reduced[pair_type] = sam[sam['pid_PIDaYplane_pCandidate']>PIDa_p_min]
            #}

            # replace the cut on PIDa to a cut on chi2_proton
            elif cut == 'Chi2Proton': #{
                reduced[pair_type] = sam[ (sam['pidcali_Chi2ProtonYplane_muCandidate']>Chi2Proton_muCandidate_min)
                             &(sam['pidcali_Chi2ProtonYplane_muCandidate']<Chi2Proton_muCandidate_max)
                             &(sam['pidcali_Chi2ProtonYplane_pCandidate']>Chi2Proton_pCandidate_min)
                             &(sam['pidcali_Chi2ProtonYplane_pCandidate']<Chi2Proton_pCandidate_max)]
        
            #}
            elif cut == 'Nflashes':#{
                reduced[pair_type] = sam[(sam['Nflashes']>0)]
            #}
            elif cut == 'ClosestFlash':#{
                reduced[pair_type] = sam[(sam['Nflashes']>0)
                                         &(sam['ClosestFlash_TotalPE'] > minPEcut)
                                         &(sam['ClosestFlash_YZdistance'] < maxdYZcut)]
            #}
            # replace the cut on ClosestFlash to a cut on MatchedFlash
            elif cut == 'MatchedFlash':#{
                reduced[pair_type] = sam[(sam['Nflashes']>0)
                             &(sam['MatchedFlash_TotalPE'] > minPEcut)
                             &(sam['MatchedFlash_YZdistance'] < maxdYZcut)]
            #}

            elif cut == 'length':#{
                reduced[pair_type] = sam[sam['l_muCandidate'] > sam['l_pCandidate']]
            #}

            elif cut == 'non-collinearity':#{
                reduced[pair_type] = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
            #}
            elif cut == 'vertex activity':#{
                R_str,box_str = 'RdQaroundVertex','[%d wires x %d ticks]'%(opt_box[0],opt_box[1])
                Ru,Rv,Ry = R_str+'[plane 0]'+box_str,R_str+'[plane 1]'+box_str,R_str+'[plane 2]'+box_str
                reduced[pair_type] = sam[(sam[Ru]==1) | (sam[Rv]==1) | (sam[Ry]==1)
                                         |
                                         (np.sqrt( np.square(sam[Ru]-1) + np.square(sam[Rv]-1) + np.square(sam[Ry]-1) )
                                          <= r_max_RdQ_CC1p) ]
            #}
            elif cut == 'delta phi':#{
                reduced[pair_type] = sam[np.abs(sam['delta_phi']-180.)<delta_Delta_phi]
            #}
            elif cut == 'Pt no delta phi':#{
                sam = reducedSamples['vertex activity'][pair_type]
                reduced[pair_type] = sam[sam['reco_Pt']<Pt_max]
            #}
            elif cut == 'Pt & delta phi':#{
                reduced[pair_type] = sam[(sam['reco_Pt']<Pt_max)
                                         &(np.abs(sam['delta_phi']-180.)<delta_Delta_phi)]
            #}
            elif cut == 'theta_pq & delta phi':#{
                sam = reducedSamples['delta phi'][pair_type]
                reduced[pair_type] = sam[sam['reco_theta_pq']<theta_pq_max]
            #}
            elif cut == 'tight Pt':#{
                sam = reducedSamples['Pt & delta phi'][pair_type]
                reduced[pair_type] = sam[(sam['reco_Pt']<0.15)
                                         &(np.abs(sam['delta_phi']-180.)<delta_Delta_phi)]
            #}
        #}
        reducedSamples[cut] = reduced
        pureffOverlay = get_pureff_cut(OverlaySamples=OverlaySamples,pureff=pureffOverlay,cut_name=cut
                                       ,reduced=reduced,OriginalOverlaySamples=reducedSamples['no cut'])
        Noverlay = get_Noverlay(reducedSamples=reducedSamples,cut_name=cut,debug=debug,overlay_scaling=overlay_scaling )
        numbers = numbers.append(pd.DataFrame({ r'$N_{Overlay}$':Noverlay['Overlay']
                                              ,r'${\epsilon}_{Overlay}$ [%]':100.*Noverlay['eff Overlay']
                                              
                                              ,r'$N_{cosmic}$':Noverlay['cosmic']
                                              ,r'$N_{other pairs}$':Noverlay['other pairs']
                                              ,r'$N_{1mu-1p}$':Noverlay['1mu-1p']
                                              ,r'$N_{CC1p0pi}$':Noverlay['CC1p0pi']
                                              ,r'$N_{CC1p}$':Noverlay['CC1p']

                                              ,r'${\epsilon}_{cosmic}$ [%]':100*Noverlay['eff cosmic']
                                              ,r'${\epsilon}_{CC1p0pi}$ [%]':100*Noverlay['eff CC1p0pi']
                                              ,r'${\epsilon}_{CC1p}$ [%]':100*Noverlay['eff CC1p']
                                              ,r'${\epsilon}_{1mu-1p}$ [%]':100*Noverlay['eff 1mu-1p']
                                              ,r'${\epsilon}_{other pairs}$ [%]':100*Noverlay['eff other pairs']
                                           
                                              ,r'${\mathcal{p}}_{1mu-1p}$ [%]':100.*Noverlay['pur 1mu-1p']
                                              ,r'${\mathcal{p}}_{CC1p0pi}$ [%]':100.*Noverlay['pur CC1p0pi']
                                              ,r'${\mathcal{p}}_{CC1p}$ [%]':100.*Noverlay['pur CC1p']

                                              # scaling and re-weighting
                                              ,r'$N_{Overlay scaled}$':Noverlay['Overlay scaled']
                                              ,r'${\epsilon}_{Overlay scaled}$ [%]':100*Noverlay['eff Overlay scaled']
                                              
                                              ,r'$N_{cosmic scaled}$':Noverlay['cosmic scaled']
                                              ,r'$N_{other pairs scaled}$':Noverlay['other pairs scaled']
                                              ,r'$N_{1mu-1p scaled}$':Noverlay['1mu-1p scaled']
                                              ,r'$N_{CC1p0pi scaled}$':Noverlay['CC1p0pi scaled']
                                              ,r'$N_{CC1p scaled}$':Noverlay['CC1p scaled']

                                          },index=[cut]))
    #}
    if debug>1: os.system('say "I have completed applying % cuts to the overlay samples, from %s to %s".'%(len(cuts_order),cuts_order[0],cuts_order[-1]))
    return reducedSamples,pureffOverlay,numbers
#}
# -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- - -- - - -- -- - -- -





# ------------------------------------------------
# April-4 (last edit Aug-12,2018)
def get_Noverlay(reducedSamples=None,cut_name=''
                 ,overlay_scaling=None
                 ,debug=0
                 ):
    # @return the number of events in each subsample of the overlay, POT-normalized
    N = dict()
    for pair_type in pair_types: #{
        N[pair_type+' original'] = float(len(reducedSamples['no cut'][pair_type]))
        N[pair_type] = float(len(reducedSamples[cut_name][pair_type]))
        N['eff '+pair_type] = N[pair_type]/N[pair_type+' original']
    #}
    N['MC'] = N['1mu-1p'] + N['other pairs']
    N['Overlay'] = N['cosmic'] + N['MC']
    N['MC original'] = N['1mu-1p original'] + N['other pairs original']
    N['Overlay original'] = N['cosmic original'] + N['MC original']
    for pair_type in pair_types: N['pur '+pair_type] = N[pair_type]/N['Overlay']
    N['eff Overlay'] = N['Overlay']/N['Overlay original']
                                     

    # scaling and re-weighting
    for pair_type in pair_types:#{
        N[pair_type+' original scaled'] = overlay_scaling[pair_type]*N[pair_type+' original'] if overlay_scaling else N[pair_type+' original']
        N[pair_type+' scaled'] = overlay_scaling[pair_type]*N[pair_type] if overlay_scaling else N[pair_type]
        N['eff '+pair_type+' scaled'] = N[pair_type+' scaled']/N[pair_type+' original scaled']
    #}
    N['MC scaled'] = N['1mu-1p scaled'] + N['other pairs scaled']
    N['Overlay scaled'] = N['cosmic scaled'] + N['MC scaled']
    N['MC original scaled'] = N['1mu-1p original scaled'] + N['other pairs original scaled']
    N['Overlay original scaled'] = N['cosmic original scaled'] + N['MC original scaled']
    for pair_type in pair_types: N['pur '+pair_type+' scaled'] = N[pair_type+' scaled']/N['Overlay scaled']
    N['eff Overlay scaled'] = N['Overlay scaled']/N['Overlay original scaled']

    if debug: print 'Noverlay in',cut_name,'\n',N
    return N
# ------------------------------------------------


# ------------------------------------------------
# July-24, 2018 (last edit Aug-16,2018)
def load_samples(date='2018_04_28'
                 ,filename='ecohen_physical_files_adi_prodgenie_bnb_nu_uboone_overlay_cosmic_data_100K_reco2_2018_02_17_vertices'
                 ,only_in_FV=False):
    '''
        return:
                samples: 
                dict() of pairs, pd.DataFrame() of different into pair_types
        '''
    pairs = pd.read_csv(vertices_files_path+'/'+date+'/'+filename+'.csv')
    pairsFV = sample_in_FV(pairs)
    print len(pairs),'ccqe candidate pairs,',len(pairsFV),'in FV'
    samples = dict()
    for pair_type in pair_types:#{
        if only_in_FV: samples[pair_type] = pairsFV[pairsFV[pair_type]==True]
        else: samples[pair_type] = pairs[pairs[pair_type]==True]
        if 'CC1p' in pair_type: print_line()
        print len(samples[pair_type]),'are '+pair_type+', %.1f'%(100.*float(len(samples[pair_type]))/len(pairsFV))+'%'
    #}
    print "I finished loading overlay samples. We have in total %d pairs"%len(pairs)
    return samples
# ------------------------------------------------


# ------------------------------------------------
# sep-12
def get_pair_hpars(index):#{
    pair_type = pair_types[index];
    label = MClabels[index]; 
    cmap = MCcmaps[index]; 
    color = MCcolors[index]
    return pair_type,label,cmap,color
#}
# ------------------------------------------------



# ------------------------------------------------
# Aug-30, 2017
def get_Nreduced_MCbnbDATAcosmic(reduced_MCbnbDATAcosmic = dict()):
    Noriginal , Nreduced , freduced = dict() , dict() , dict()
    for pair_type in pair_types:
        sam = MCbnbDATAcosmicSamples[pair_type]
        Noriginal[pair_type] = len(MCbnbDATAcosmicSamples[pair_type])
        Nreduced[pair_type] = float(len(reduced_MCbnbDATAcosmic[pair_type]))
        freduced[pair_type] = 100.0 * Nreduced[pair_type]/Noriginal[pair_type]
    return Nreduced , freduced
# ------------------------------------------------


      
# ------------------------------------------------
# Oct. 13, 2017 (last edit Aug-15,2018)
def get_pureff_MCbnbDATAcosmic_cut(cut_name = 'PIDa', cut_label=None , reduced_MCbnbDATAcosmic = dict()):
    ''' 
        return
        eff (mu-p) , pur (mu-p), eff (CC1p0pi) , pur (CC1p0pi)
    '''
    global pureff_MCbnbDATAcosmic
    eff = dict()
    pur = dict()
    Nreduced_MCbnbDATAcosmic , freduced_MCbnbDATAcosmic = get_Nreduced_MCbnbDATAcosmic(reduced_MCbnbDATAcosmic=reduced_MCbnbDATAcosmic)
    Ntot = (Nreduced_MCbnbDATAcosmic['1mu-1p']+Nreduced_MCbnbDATAcosmic['cosmic']+Nreduced_MCbnbDATAcosmic['other pairs'])
    
    eff['1mu-1p'] = freduced_MCbnbDATAcosmic['1mu-1p']
    pur['1mu-1p'] = 100.*Nreduced_MCbnbDATAcosmic['1mu-1p']/Ntot if Ntot>0 else 0
    
    eff['CC1p0pi'] = freduced_MCbnbDATAcosmic['CC1p0pi']
    pur['CC1p0pi'] = 100.*Nreduced_MCbnbDATAcosmic['CC1p0pi']/Ntot if Ntot>0 else 0

    eff['CC1p'] = freduced_MCbnbDATAcosmic['CC1p']
    pur['CC1p'] = 100.*Nreduced_MCbnbDATAcosmic['CC1p']/Ntot if Ntot>0 else 0

    pureff_MCbnbDATAcosmic_cut = pd.DataFrame({'label':cut_label
                               ,'$\mu p$ eff.':'%.1f'%eff['1mu-1p']+'%'
                               ,'$\mu p$ pur.':'%.1f'%pur['1mu-1p']+'%'
                               ,'CC$0\pi 1 p$ eff.':'%.1f'%freduced_MCbnbDATAcosmic['CC1p0pi']+'%'
                               ,'CC$0\pi 1 p$ pur.':'%.1f'%(100.*Nreduced_MCbnbDATAcosmic['CC1p0pi']/Ntot if Ntot>0 else 0)+'%'
                               ,'CC$1 p$ eff.':'%.1f'%freduced_MCbnbDATAcosmic['CC1p']+'%'
                               ,'CC$1 p$ pur.':'%.1f'%(100.*Nreduced_MCbnbDATAcosmic['CC1p']/Ntot if Ntot>0 else 0)+'%'
                                              }
                               , index=[cut_name]
                              )
    for pair_type in pair_types: pureff_MCbnbDATAcosmic_cut[pair_type] = '%.1f'%freduced_MCbnbDATAcosmic[pair_type]+'%' +' (%.0f)'%Nreduced_MCbnbDATAcosmic[pair_type]
    pureff_MCbnbDATAcosmic = pureff_MCbnbDATAcosmic.append(pureff_MCbnbDATAcosmic_cut)
    reduced_MCbnbDATAcosmicSamples[cut_name] = reduced_MCbnbDATAcosmic  
    Ntot = Nreduced_MCbnbDATAcosmic['1mu-1p']+Nreduced_MCbnbDATAcosmic['cosmic']+Nreduced_MCbnbDATAcosmic['other pairs']
    return freduced_MCbnbDATAcosmic['1mu-1p'],(100.*Nreduced_MCbnbDATAcosmic['1mu-1p']/Ntot if Ntot>0 else 0),freduced_MCbnbDATAcosmic['CC1p'],(100.*Nreduced_MCbnbDATAcosmic['CC1p']/Ntot if Ntot>0 else 0)
# ------------------------------------------------

# ------------------------------------------------
# Oct. 13, 2017
def get_pureff_MCbnbDATAcosmic_numbers(cut_name = 'PIDa', cut_label=None , reduced_MCbnbDATAcosmic = dict()):
    global pureff_MCbnbDATAcosmic_numbers
    Nreduced_MCbnbDATAcosmic , freduced_MCbnbDATAcosmic = get_Nreduced_MCbnbDATAcosmic(reduced_MCbnbDATAcosmic=reduced_MCbnbDATAcosmic)
    Ntot = Nreduced_MCbnbDATAcosmic['cosmic']+Nreduced_MCbnbDATAcosmic['other pairs']+Nreduced_MCbnbDATAcosmic['1mu-1p']
    pureff_MCbnbDATAcosmic_numbers_cut = pd.DataFrame({'cut name':cut_name,
                                                      'cut label':cut_label,
                                                      'cosmic':Nreduced_MCbnbDATAcosmic['cosmic'],
                                                      'other pairs':Nreduced_MCbnbDATAcosmic['other pairs'],
                                                      '\mup':Nreduced_MCbnbDATAcosmic['1mu-1p'],
                                                      '\CCIpOpi':Nreduced_MCbnbDATAcosmic['CC1p0pi'],
                                                      'eff \mup':freduced_MCbnbDATAcosmic['1mu-1p'],
                                                      'eff \CCIpOpi':freduced_MCbnbDATAcosmic['CC1p0pi'],
                                                      'pur \mup':float(100*Nreduced_MCbnbDATAcosmic['1mu-1p'])/Ntot if Ntot>0 else 0,
                                                      'pur \CCIpOpi':float(100*Nreduced_MCbnbDATAcosmic['CC1p0pi'])/Ntot if Ntot>0 else 0}
                                                      , index=[cut_name]
                                                      )
    pureff_MCbnbDATAcosmic_numbers = pureff_MCbnbDATAcosmic_numbers.append(pureff_MCbnbDATAcosmic_numbers_cut)
# ------------------------------------------------


#---------------------------------------------------------------------------------------------
# July-11, 2017 (last edited May-15 2018)
def plot_feature_pairs(reduced_samples=None,cut_name='PIDa',
                       var='l_long',x_label='$l_{long}$ [cm]',mul=1,
                       bins=np.linspace(0,300,50),
                       figsize=(12,8),legend_fontsize=25,fontsize=25,
                       do_add_legend=False,legend_loc='upper center',
                       ticks_color='black'):
    fig,ax = plt.subplots(figsize=figsize)
    max_h=0
    text_colors=[]
    for i,(pair_type,label,cmap,color) in enumerate(zip(pair_types,MClabels,MCcmaps,MCcolors)):
        sample = reduced_samples[pair_type]
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
# July-11, 2017 (edited May-16, 2018)
def plot_cut_samples (reduced_samples=None,markers_size=5,
                      cut_name='maximal distance between tracks',mul=1,
                      cut_var ='distance',
                      cut_type= 'max',
                      x_label = 'maximal tracks distance [cm]', y_label='% of sample',
                      xcenter=0,figsize=(12,8),fontsize=25,
                      xmin=0.1, xmax=10 , Nbins=10, do_add_legend=True, legend_loc='bbox',legend_fontsize=25,
                      ticks_color='black'):
    fig,ax=plt.subplots(figsize=figsize)
    for i,(pair_type,label,cmap,color) in enumerate(zip(pair_types,MClabels,MCcmaps,MCcolors)):
        sample = reduced_samples[pair_type]
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
        for color,text in zip(MCcolors,leg.get_texts()):
            text.set_color(color)
    ax.set_ylim(0,101)
    ax.set_xlim(xmin,xmax)
    set_axes(ax,x_label=x_label,y_label=y_label,fontsize=fontsize,ticks_color=ticks_color,xticks=np.linspace(xmin,xmax,7),yticks=[25,50,75,100])
    ax.grid(linestyle='--',alpha=0.75)
    plt.tight_layout()
    return ax,leg
#---------------------------------------------------------------------------------------------




# ------------------------------------------------
# Sep 9
def load_MCbnbMCcosmicSamples(filename='prodgenie_bnb_nu_cosmic_uboone_mcc8.2_reco2_vertices'):
    '''
        return:
        MCbnbMCcosmicPairsFV:  pandas.DataFrame() of prodgenie_bnb_nu_cosmic_uboone_mcc8.2_reco2_vertices in FV
        samples: dict() of pairs broken into pair_types
        '''
    global MCbnbMCcosmicSamples
    pairs = pd.read_csv(vertices_files_path+filename+'.csv')
    MCbnbMCcosmicPairsFV = sample_in_FV(pairs)
    print len(pairs),'pairs from MC-BNB + cosmic MC overlay'
    print len(MCbnbMCcosmicPairsFV),'pairs in FV'
    for pair_type in pair_types:#{
        MCbnbMCcosmicSamples[pair_type] = MCbnbMCcosmicPairsFV[MCbnbMCcosmicPairsFV[pair_type]==True]
        Ntype = len(MCbnbMCcosmicSamples[pair_type])
        if pair_type=='CC1p0pi': print_line()
        print Ntype,'are '+pair_type+', %.1f'%(100.*float(Ntype)/len(MCbnbMCcosmicPairsFV))+'%'
    #}
    return MCbnbMCcosmicPairsFV, MCbnbMCcosmicSamples
# ------------------------------------------------




# ------------------------------------------------
# Sep-9, 2017
def get_Nreduced_MCbnbMCcosmic(reduced_MCbnbMCcosmic = dict()):
    Noriginal , Nreduced , freduced = dict() , dict() , dict()
    for pair_type in pair_types:
        sam = MCbnbMCcosmicSamples[pair_type]
        Noriginal[pair_type] = len(MCbnbMCcosmicSamples[pair_type])
        Nreduced[pair_type] = float(len(reduced_MCbnbMCcosmic[pair_type]))
        freduced[pair_type] = 100.0 * Nreduced[pair_type]/Noriginal[pair_type]
    return Nreduced , freduced
# ------------------------------------------------



# ------------------------------------------------
# Oct. 13, 2017
def get_pureff_MCbnbMCcosmic_cut(cut_name = 'PIDa', cut_label=None , reduced_MCbnbMCcosmic = dict()):
    '''
        return
        eff (mu-p) , pur (mu-p), eff (CC1p0pi) , pur (CC1p0pi)
        '''
    global pureff_MCbnbMCcosmic
    eff = dict()
    pur = dict()
    Nreduced_MCbnbMCcosmic , freduced_MCbnbMCcosmic = get_Nreduced_MCbnbMCcosmic(reduced_MCbnbMCcosmic=reduced_MCbnbMCcosmic)
    Ntot = (Nreduced_MCbnbMCcosmic['1mu-1p']+Nreduced_MCbnbMCcosmic['cosmic']+Nreduced_MCbnbMCcosmic['other pairs'])
    
    eff['1mu-1p'] = freduced_MCbnbMCcosmic['1mu-1p']
    pur['1mu-1p'] = 100.*Nreduced_MCbnbMCcosmic['1mu-1p']/Ntot if Ntot>0 else 0
    
    eff['CC1p0pi'] = freduced_MCbnbMCcosmic['CC1p0pi']
    pur['CC1p0pi'] = 100.*Nreduced_MCbnbMCcosmic['CC1p0pi']/Ntot if Ntot>0 else 0
    
    pureff_MCbnbMCcosmic_cut = pd.DataFrame({'label':cut_label
                                            ,'$\mu p$ eff.':'%.1f'%eff['1mu-1p']+'%'
                                            ,'$\mu p$ pur.':'%.1f'%pur['1mu-1p']+'%'
                                            ,'CC$0\pi 1 p$ eff.':'%.1f'%freduced_MCbnbMCcosmic['CC1p0pi']+'%'
                                            ,'CC$0\pi 1 p$ pur.':'%.1f'%(100.*Nreduced_MCbnbMCcosmic['CC1p0pi']/Ntot if Ntot>0 else 0)+'%'}
                                            , index=[cut_name]
                                            )

    for pair_type in pair_types: pureff_MCbnbMCcosmic_cut[pair_type] = '%.1f'%freduced_MCbnbMCcosmic[pair_type]+'%' +' (%.0f)'%Nreduced_MCbnbMCcosmic[pair_type]
    pureff_MCbnbMCcosmic = pureff_MCbnbMCcosmic.append(pureff_MCbnbMCcosmic_cut)
    reduced_MCbnbMCcosmicSamples[cut_name] = reduced_MCbnbMCcosmic
    Ntot = Nreduced_MCbnbMCcosmic['1mu-1p']+Nreduced_MCbnbMCcosmic['cosmic']+Nreduced_MCbnbMCcosmic['other pairs']
    
    return freduced_MCbnbMCcosmic['1mu-1p'],(100.*Nreduced_MCbnbMCcosmic['1mu-1p']/Ntot if Ntot>0 else 0),freduced_MCbnbMCcosmic['CC1p0pi'],(100.*Nreduced_MCbnbMCcosmic['CC1p0pi']/Ntot if Ntot>0 else 0)
# ------------------------------------------------


# ------------------------------------------------
# Oct. 13, 2017
def get_pureff_MCbnbMCcosmic_numbers(cut_name = 'PIDa', cut_label=None , reduced_MCbnbMCcosmic = dict()):
    global pureff_MCbnbMCcosmic_numbers
    Nreduced_MCbnbMCcosmic , freduced_MCbnbMCcosmic = get_Nreduced_MCbnbMCcosmic(reduced_MCbnbMCcosmic=reduced_MCbnbMCcosmic)
    Ntot = Nreduced_MCbnbMCcosmic['cosmic']+Nreduced_MCbnbMCcosmic['other pairs']+Nreduced_MCbnbMCcosmic['1mu-1p']
    pureff_MCbnbMCcosmic_numbers_cut = pd.DataFrame({'cut name':cut_name,
                                                      'cut label':cut_label,
                                                      'cosmic':Nreduced_MCbnbMCcosmic['cosmic'],
                                                      'other pairs':Nreduced_MCbnbMCcosmic['other pairs'],
                                                      '\mup':Nreduced_MCbnbMCcosmic['1mu-1p'],
                                                      '\CCIpOpi':Nreduced_MCbnbMCcosmic['CC1p0pi'],
                                                      'eff \mup':freduced_MCbnbMCcosmic['1mu-1p'],
                                                      'eff \CCIpOpi':freduced_MCbnbMCcosmic['CC1p0pi'],
                                                      'pur \mup':float(100*Nreduced_MCbnbMCcosmic['1mu-1p'])/Ntot if Ntot>0 else 0,
                                                      'pur \CCIpOpi':float(100*Nreduced_MCbnbMCcosmic['CC1p0pi'])/Ntot if Ntot>0 else 0}
                                                      , index=[cut_name]
                                                      )
    pureff_MCbnbMCcosmic_numbers = pureff_MCbnbMCcosmic_numbers.append(pureff_MCbnbMCcosmic_numbers_cut)
# ------------------------------------------------




from ccqe_notebook_tools import * 


vertices_files_path = '/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/ccqe_candidates/'
MCsamples=dict()
reduced_MCsamples=dict(dict())
pur_eff = pd.DataFrame()
pur_eff_numbers = pd.DataFrame()


# ------------------------------------------------
# Aug 30
def apply_cuts( PIDa_p_min=8
               , delta_theta_12=60  # deg.
               , delta_Delta_phi=40 # deg.
               , theta_pq_max=25    # deg.
               , Pt_max=0.35        # GeV/c
               , i_optimal_box_size=9 , Rmin = (0.25,0.25,0.6) # U,V,Y
               ):
#    # initiate the pur_eff dataframe
#    pur_eff = pd.DataFrame()
#    pur_eff_numbers = pd.DataFrame()
    # --- -- --- - -- -- --- --
    reduced = dict()
    for pair_type in pair_types: reduced[pair_type] = MCsamples[pair_type]
    reduced_MCsamples['no cut'] = reduced

    # before cuts
    reduced = dict()
    for pair_type in pair_types:
        sam = reduced_MCsamples['no cut'][pair_type]
        reduced[pair_type] = sam
    get_pur_eff_cut(cut_name = 'no cut', reduced = reduced)
    get_pur_eff_numbers(cut_name = 'no cut', reduced = reduced)

    # cut 1: PIDa
    reduced = dict()
    for pair_type in pair_types:
        sam = reduced_MCsamples['no cut'][pair_type]
        reduced[pair_type] = sam[sam['PIDa_assigned_proton']>PIDa_p_min]
    get_pur_eff_cut(cut_name = 'PIDa', cut_label='${PID}_a>%.0f$'%PIDa_p_min, reduced = reduced)
    get_pur_eff_numbers(cut_name = '\CutPIDa', cut_label='${PID}_a>%.0f$'%PIDa_p_min, reduced = reduced)

    # cut 2: require that the longer track is the one with larger PIDa
    reduced = dict()
    for pair_type in pair_types:
        sam = reduced_MCsamples['PIDa'][pair_type]
        reduced[pair_type] = sam[sam['PIDa_long'] < sam['PIDa_short']]
    get_pur_eff_cut(cut_name = 'length', cut_label='$l_{\\mu}>l_{p}$', reduced = reduced)
    get_pur_eff_numbers(cut_name = '\Cutlmup', cut_label='$l_{\\mu}>l_{p}$', reduced = reduced)

    # cut 3: |\theta_{1,2}-90^0|<60^0$
    reduced = dict()
    for pair_type in pair_types:
        sam = reduced_MCsamples['length'][pair_type]
        reduced[pair_type] = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
    get_pur_eff_cut(cut_name='non-collinearity' , cut_label='$|\theta_{1,2}-90^0|<%.0f^0$'%delta_theta_12, reduced = reduced)
    get_pur_eff_numbers(cut_name = 'non-collinearity', cut_label='$|\theta_{1,2}-90^0|<%.0f^0$'%delta_theta_12, reduced = reduced)

    # cut 4: vertex $\Delta Q$
    N_box_sizes = 30
    MinNwiresBox,dNwiresBox = 5,5
    MinNticksBox,dNticksBox = 10,10
    NwiresBox,NticksBox=[],[]
    for i_box_size in range(N_box_sizes):#{
        NwiresBox.append(MinNwiresBox + i_box_size * dNwiresBox)
        NticksBox.append(MinNticksBox + i_box_size * dNticksBox)
    #}
    Ru = 'RdQaroundVertex[plane 0][%d wires x %d ticks]'%(NwiresBox[i_optimal_box_size],NticksBox[i_optimal_box_size])
    Rv = 'RdQaroundVertex[plane 1][%d wires x %d ticks]'%(NwiresBox[i_optimal_box_size],NticksBox[i_optimal_box_size])
    Ry = 'RdQaroundVertex[plane 2][%d wires x %d ticks]'%(NwiresBox[i_optimal_box_size],NticksBox[i_optimal_box_size])
    RuMin = Rmin[0]-0.01
    RvMin = Rmin[1]-0.01
    RyMin = Rmin[2]-0.01
    cut_label='$(%.2f,%.2f,%.2f)<R_{\Delta Q}^{(U,V,Y)}$'%(RuMin,RvMin,RyMin)
    reduced = dict()
    for pair_type in pair_types:
        sam = reduced_MCsamples['non-collinearity'][pair_type]
        reduced[pair_type] = sam[(sam[Ru]>RuMin)&(sam[Ru]<1.01)
                                 &(sam[Rv]>RvMin)&(sam[Rv]<1.01)
                                 &(sam[Ry]>RyMin)&(sam[Ry]<1.01)]
    get_pur_eff_cut(cut_name = 'vertex activity' ,cut_label=cut_label , reduced = reduced)
    get_pur_eff_numbers(cut_name = 'vertex activity', cut_label=cut_label, reduced = reduced)


    # cut 5: $\Delta phi$
    reduced = dict()
    for pair_type in pair_types:
        sam = reduced_MCsamples['vertex activity'][pair_type]
        reduced[pair_type] = sam[np.abs(sam['delta_phi']-180.)<delta_Delta_phi]
    get_pur_eff_cut(cut_name = 'delta phi', cut_label='$|\Delta \phi - \pi|<%.0f^0$'%delta_Delta_phi, reduced = reduced)
    get_pur_eff_numbers(cut_name = '\CutDeltaPhi', cut_label='$|\Delta \phi - \pi|<%.0f^0$'%delta_Delta_phi, reduced = reduced)


    # cut 6: $\theta_{pq}<25$
    reduced = dict()
    for pair_type in pair_types:
        sam = reduced_MCsamples['delta phi'][pair_type]
        reduced[pair_type] = sam[sam['reco_theta_pq']<theta_pq_max]
    get_pur_eff_cut(cut_name ='theta_pq' , cut_label= '$\theta_{pq}<%.0f^0$'%theta_pq_max, reduced = reduced)
    get_pur_eff_numbers(cut_name ='theta_pq' , cut_label= '$\\theta_{pq}<%.0f^0$'%theta_pq_max, reduced = reduced)


    # modified cut 6: $p_{t}<0.35$
    reduced = dict()
    for pair_type in pair_types:
        sam = reduced_MCsamples['delta phi'][pair_type]
        reduced[pair_type] = sam[sam['reco_Pt']<Pt_max]
    get_pur_eff_cut(cut_name='soft Pt', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced = reduced)
    get_pur_eff_numbers(cut_name='soft Pt', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced = reduced)


    # tight Pt cut for good Ev reconstruction
    reduced = dict()
    for pair_type in pair_types:
        sam = reduced_MCsamples['theta_pq'][pair_type]
        reduced[pair_type] = sam[sam['reco_Pt']<0.15]
    get_pur_eff_cut(cut_name ='tight Pt', cut_label= '$p_{t}<0.15$ GeV/c', reduced = reduced)
    get_pur_eff_numbers(cut_name ='tight Pt', cut_label= '$p_{t}<0.15$ GeV/c', reduced = reduced)


    
    return pur_eff,pur_eff_numbers
# ------------------------------------------------



# ------------------------------------------------
# Aug 30
def load_pairs_as_samples():
    '''
    return: 
           pairsFV:  pandas.DataFrame() of prodgenie_bnb_nu_uboone_overlay_mcc8_reco2_vertices in FV
           samples: dict() of pairs broken into pair_types
    '''
    global MCsamples
    pairs = pd.read_csv(vertices_files_path+'prodgenie_bnb_nu_uboone_overlay_mcc8_reco2_vertices.csv')    
    pairsFV = sample_in_FV(pairs)
    print len(pairs),'pairs from MC-BNB + cosmic DATA overlay'
    print len(pairsFV),'pairs in FV'
    for pair_type in pair_types:
        MCsamples[pair_type] = pairsFV[pairsFV[pair_type]==True]
        Ntype = len(MCsamples[pair_type])
        if pair_type=='CC 1p 0pi': print_line()
        print Ntype,'are '+pair_type+', %.1f'%(100.*float(Ntype)/len(pairsFV))+'%'

    return pairsFV, MCsamples
# ------------------------------------------------


        
# ------------------------------------------------
# Aug-30, 2017
def get_Nreduced(reduced = dict()):
    Noriginal , Nreduced , freduced = dict() , dict() , dict()
    for pair_type in pair_types:
        sam = MCsamples[pair_type]
        Noriginal[pair_type] = len(MCsamples[pair_type])
        Nreduced[pair_type] = float(len(reduced[pair_type]))
        freduced[pair_type] = 100.0 * Nreduced[pair_type]/Noriginal[pair_type]
    return Nreduced , freduced
# ------------------------------------------------
           
    
    
# ------------------------------------------------
# Aug-30, 2017
def get_pur_eff_cut(cut_name = 'PIDa', cut_label=None , reduced = dict()):
    ''' 
        return
        eff (mu-p) , pur (mu-p), eff (CC 1p 0pi) , pur (CC 1p 0pi)
    '''
    
    global pur_eff
    eff = dict()
    pur = dict()
    Nreduced , freduced = get_Nreduced(reduced=reduced)
    Ntot = (Nreduced['1mu-1p']+Nreduced['cosmic']+Nreduced['other pairs'])
    
    eff['1mu-1p'] = freduced['1mu-1p']
    pur['1mu-1p'] = 100.*Nreduced['1mu-1p']/Ntot if Ntot>0 else 0
    
    eff['CC 1p 0pi'] = freduced['CC 1p 0pi']
    pur['CC 1p 0pi'] = 100.*Nreduced['CC 1p 0pi']/Ntot if Ntot>0 else 0
    
    pur_eff_cut = pd.DataFrame({'label':cut_label
                               ,'$\mu p$ eff.':'%.1f'%eff['1mu-1p']+'%'
                               ,'$\mu p$ pur.':'%.1f'%pur['1mu-1p']+'%'
                               ,'CC$0\pi 1 p$ eff.':'%.1f'%freduced['CC 1p 0pi']+'%'
                               ,'CC$0\pi 1 p$ pur.':'%.1f'%(100.*Nreduced['CC 1p 0pi']/Ntot if Ntot>0 else 0)+'%'}
                               , index=[cut_name]
                              )
    for pair_type in pair_types: pur_eff_cut[pair_type] = '%.1f'%freduced[pair_type]+'%' +' (%.0f)'%Nreduced[pair_type]
    pur_eff = pur_eff.append(pur_eff_cut)
    reduced_MCsamples[cut_name] = reduced  
    Ntot = Nreduced['1mu-1p']+Nreduced['cosmic']+Nreduced['other pairs']
    return freduced['1mu-1p'],(100.*Nreduced['1mu-1p']/Ntot if Ntot>0 else 0),freduced['CC 1p 0pi'],(100.*Nreduced['CC 1p 0pi']/Ntot if Ntot>0 else 0)
# ------------------------------------------------

# ------------------------------------------------
# Aug-30, 2017
def get_pur_eff_numbers(cut_name = 'PIDa', cut_label=None , reduced = dict()):
    global pur_eff_numbers
    Nreduced , freduced = get_Nreduced(reduced=reduced)
    Ntot = Nreduced['cosmic']+Nreduced['other pairs']+Nreduced['1mu-1p']
    pur_eff_numbers_cut = pd.DataFrame({'cut name':cut_name,
                                       'cut label':cut_label,
                                        'cosmic':Nreduced['cosmic'],                                
                                        'other pairs':Nreduced['other pairs'],                                
                                        '\mup':Nreduced['1mu-1p'],                               
                                        '\CCIpOpi':Nreduced['CC 1p 0pi'],
                                        'eff \mup':freduced['1mu-1p'],                               
                                        'eff \CCIpOpi':freduced['CC 1p 0pi'],
                                        'pur \mup':float(100*Nreduced['1mu-1p'])/Ntot if Ntot>0 else 0,
                                        'pur \CCIpOpi':float(100*Nreduced['CC 1p 0pi'])/Ntot if Ntot>0 else 0}
                                       , index=[cut_name]
                                      )
    pur_eff_numbers = pur_eff_numbers.append(pur_eff_numbers_cut)
# ------------------------------------------------




# ------------------------------------------------
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
    for i,(pair_type,label,cmap,color) in enumerate(zip(pair_types,labels,cmaps,colors)):
        sample = reduced_MCsamples[cut_name][pair_type]
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
def plot_cut_samples (reduced_cut_name='${PID}_A$',markers_size=5,
                      cut_name='maximal distance between tracks',mul=1,
                      cut_var ='distance',
                      cut_type= 'max',
                      x_label = 'maximal tracks distance [cm]', y_label='% of sample',
                      xcenter=0,figsize=(12,8),fontsize=25,
                      xmin=0.1, xmax=10 , Nbins=10, do_add_legend=True, legend_loc='bbox',legend_fontsize=25,
                      ticks_color='black'):
    fig,ax=plt.subplots(figsize=figsize)
    for i,(pair_type,label,cmap,color) in enumerate(zip(pair_types,labels,cmaps,colors)):
        sample = reduced_MCsamples[reduced_cut_name][pair_type]
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
        for color,text in zip(colors,leg.get_texts()):
            text.set_color(color)
    ax.set_ylim(0,101)
    ax.set_xlim(xmin,xmax)    
    set_axes(ax,x_label=x_label,y_label=y_label,fontsize=fontsize,ticks_color=ticks_color,xticks=np.linspace(xmin,xmax,7),yticks=[25,50,75,100])
    ax.grid(linestyle='--',alpha=0.75)
    plt.tight_layout()    
    return ax,leg
#---------------------------------------------------------------------------------------------




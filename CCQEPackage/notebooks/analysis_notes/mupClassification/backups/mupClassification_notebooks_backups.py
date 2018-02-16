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

# ------------------------------------------------
# sep-12
def get_pair_hpars(index):
    pair_type = pair_types[index];
    label = MClabels[index]; 
    cmap = MCcmaps[index]; 
    color = MCcolors[index]
    return pair_type,label,cmap,color
# ------------------------------------------------



# ------------------------------------------------
# last edit:
# Feb. 12, 2018
def apply_cuts_MCbnbDATAcosmic( PIDa_p_min=12
                               , delta_theta_12=60  # deg.
                               , delta_Delta_phi=40 # deg.
                               , theta_pq_max=25    # deg.
                               , Pt_max=0.35        # GeV/c
                               , Pmiss_max=0.3      # GeV/c
                               , opt_box=(50,100) # [Nwires x Nticks]
                               , r_max_RdQ_CC1p0pi = 0.35 # sphere in U,V,Y space, apply a cut only to CC1p0pi
                               , minPEcut = 100
                               , maxdYZcut = 200                              

                               
                               # cuts sensitivity
                               , PIDa_p_min_minus = 11
                               , PIDa_p_min_plus = 13
                               , delta_theta_12_minus = 55 
                               , delta_theta_12_plus = 65 
                               , r_max_RdQ_CC1p0pi_minus = 0.25
                               , r_max_RdQ_CC1p0pi_plus = 0.4                               
                               , delta_Delta_phi_minus = 35
                               , delta_Delta_phi_plus = 45
                               , Pt_max_minus=0.3
                               , Pt_max_plus=0.4
                               , theta_pq_max_minus=15
                               , theta_pq_max_plus=35 
               ):
    '''
        return:      
            pureff_MCbnbDATAcosmic
            pureff_MCbnbDATAcosmic_numbers
    '''
    # --- -- --- - -- -- --- --
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        reduced_MCbnbDATAcosmic[pair_type] = MCbnbDATAcosmicSamples[pair_type]
    #}
    reduced_MCbnbDATAcosmicSamples['no cut'] = reduced_MCbnbDATAcosmic

    # before cuts
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['no cut'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name = 'no cut', cut_label='no cut', reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name = 'no cut', cut_label='no cut', reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    # cut 3 before everything: |\theta_{1,2}-90^0|<60^0$
    cut_name, cut_label='non-collinearity first' ,'$|\theta_{1,2}-90^0|<%.0f^0$'%delta_theta_12
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['no cut'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    # cut 1: PIDa
    cut_name , cut_label = 'PIDa','${PID}_a>%.0f$'%PIDa_p_min
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['no cut'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['PIDa_assigned_proton']>PIDa_p_min]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    
    # -- -- -- -- -- ---- - ---- 
    # Optical filtering
    # -- -- -- -- -- ---- - ---- 
    cut_name , cut_label = 'flash', 'optical filter'
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['PIDa'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[(sam['ClosestFlash_TotalPE'] > minPEcut)
                                                 &(sam['ClosestFlash_YZdistance'] < maxdYZcut)]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name,cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)


    
    # cut 2: require that the longer track is the one with larger PIDa
    cut_name , cut_label = 'length', r'$l_{\mu}>l_{p}$'
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        #sam = reduced_MCbnbDATAcosmicSamples['PIDa'][pair_type]
        sam = reduced_MCbnbDATAcosmicSamples['flash'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['PIDa_long'] < sam['PIDa_short']]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    
    # cut 3: |\theta_{1,2}-90^0|<60^0$
    cut_name, cut_label='non-collinearity' ,'$|\theta_{1,2}-90^0|<%.0f^0$'%delta_theta_12
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['length'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    # -- - -- -- - -- -- - - -- --- 
    # VERTEX ACTIVITY 
    # -- - -- -- - -- -- - - -- --- 
    # cut 6: vertex activity (RdQ)
    # for the 1mu-1p sample we do not apply a cut, which means that we stop here for 1mu-1p
    # for the CC 1p 0pi we apply a cut:
    # a sphere around RdQ=1 of radius r_min_RdQ_CC1p0pi
    # where r_min_RdQ_CC1p0pi is taken from analysis_notes/RdQ/RdQaroundVertex_cut_selection
    cut_name , cut_label = 'vertex activity' , '$\sqrt{\sum_{p=0,1,2}(R_{\Delta Q}^{p}-1)^2}<%.2f$'%r_max_RdQ_CC1p0pi
    box_str='[%d wires x %d ticks]'%(opt_box[0],opt_box[1])
    Ru = 'RdQaroundVertex[plane 0]'+box_str
    Rv = 'RdQaroundVertex[plane 1]'+box_str
    Ry = 'RdQaroundVertex[plane 2]'+box_str
    
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{         
        sam = reduced_MCbnbDATAcosmicSamples['non-collinearity'][pair_type]
        sam = sam[(sam[Ru]==1) | (sam[Rv]==1) | (sam[Ry]==1) 
                  | 
                  (np.sqrt( np.square(sam[Ru]-1) + np.square(sam[Rv]-1) + np.square(sam[Ry]-1) ) <= r_max_RdQ_CC1p0pi) ]
        reduced_MCbnbDATAcosmic[pair_type] = sam        
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name ,cut_label=cut_label , reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    # cut 5: $\Delta phi$
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['vertex activity'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[np.abs(sam['delta_phi']-180.)<delta_Delta_phi]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name = 'delta phi', cut_label='$|\Delta \phi - \pi|<%.0f^0$'%delta_Delta_phi, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name = '\CutDeltaPhi', cut_label='$|\Delta \phi - \pi|<%.0f^0$'%delta_Delta_phi, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    # modified cut 6: $\theta_{pq}<25$
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['delta phi'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_theta_pq']<theta_pq_max]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name ='theta_pq' , cut_label= '$\theta_{pq}<%.0f^0$'%theta_pq_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name ='theta_pq' , cut_label= '$\\theta_{pq}<%.0f^0$'%theta_pq_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    # cut 5.5: $p_{t}<0.35$ without application of the $\Delta phi$ cut
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['vertex activity'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_Pt']<Pt_max]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name='Pt no Delta phi', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name='Pt no Delta phi', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    
    
    # modified cut 6: $p_{t}<0.35$
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['delta phi'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_Pt']<Pt_max]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name='soft Pt', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name='soft Pt', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    # another modified cut 6: $p_{miss}<0.3$
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['delta phi'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_Pmiss']<Pmiss_max]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name='soft Pmiss', cut_label='$p_{miss}<%.2f$ GeV/c'%Pmiss_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name='soft Pmiss', cut_label='$p_{miss}<%.2f$ GeV/c'%Pmiss_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)

    
    
    
    # tight Pt cut for good Ev reconstruction
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['theta_pq'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_Pt']<0.15]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name ='tight Pt', cut_label= '$p_{t}<0.15$ GeV/c', reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name ='tight Pt', cut_label= '$p_{t}<0.15$ GeV/c', reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)


    # cut on l_p>10 cm
    cut_name,cut_label='l_p_min_10',r'$l_p>10$ cm'
    reduced_MCbnbDATAcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['soft Pt'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['l_assigned_proton']>10]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)


    # -- -- - -- - -- - -- - -- - -- - -- - -- - -- 
    # cut-sensitivity
    # -- -- - -- - -- - -- - -- - -- - -- - -- - -- 
    # cut 1: PIDa    
    for cut_name,PIDa_p_min_cut in zip(['PIDa-','PIDa+'],[PIDa_p_min_minus,PIDa_p_min_plus]):#{
        cut_label = '${PID}_a>%.0f$'%PIDa_p_min_cut
        reduced_MCbnbDATAcosmic = dict()
        for pair_type in pair_types:#{
            sam = reduced_MCbnbDATAcosmicSamples['no cut'][pair_type]
            reduced_MCbnbDATAcosmic[pair_type] = sam[sam['PIDa_assigned_proton']>PIDa_p_min_cut]
        #}
        get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic);
        get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic);
    #}
    
    # cut 3: collinearity
    for cut_name,delta_theta_12_cut in zip(['theta12-','theta12+'],[delta_theta_12_minus,delta_theta_12_plus]):#{
        cut_label='$|\theta_{1,2}-90^0|<%.0f^0$'%delta_theta_12_cut
        reduced_MCbnbDATAcosmic = dict()
        for pair_type in pair_types:#{
            sam = reduced_MCbnbDATAcosmicSamples['length'][pair_type]
            reduced_MCbnbDATAcosmic[pair_type] = sam[np.abs(sam['theta_12']-90)<delta_theta_12_cut]
        #}
        get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
        get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic);
    #}

    # cut 6: vertex activity (RdQ)
    for cut_name,r_max_RdQ_CC1p0pi_cut in zip(['RdQ-','RdQ+'],[r_max_RdQ_CC1p0pi_minus,r_max_RdQ_CC1p0pi_plus]):#{
        cut_label = '$\sqrt{\sum_{p=0,1,2}(R_{\Delta Q}^{p}-1)^2}<%.2f$'%r_max_RdQ_CC1p0pi
        reduced_MCbnbDATAcosmic = dict()
        for pair_type in pair_types:#{         
            sam = reduced_MCbnbDATAcosmicSamples['non-collinearity'][pair_type]
            sam = sam[(sam[Ru]==1) | (sam[Rv]==1) | (sam[Ry]==1) 
                      | 
                      (np.sqrt( np.square(sam[Ru]-1) + np.square(sam[Rv]-1) + np.square(sam[Ry]-1) ) <= r_max_RdQ_CC1p0pi_cut) ]
            reduced_MCbnbDATAcosmic[pair_type] = sam        
        #}
        get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name ,cut_label=cut_label , reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
        get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic);
    #}

    # cut 5: $\Delta phi$
    for cut_name,delta_Delta_phi_cut in zip(['delta_phi-','delta_phi+'],[delta_Delta_phi_minus,delta_Delta_phi_plus]):#{
        cut_label ='$|\Delta \phi - \pi|<%.0f^0$'%delta_Delta_phi_cut
        reduced_MCbnbDATAcosmic = dict()
        for pair_type in pair_types:#{
            sam = reduced_MCbnbDATAcosmicSamples['vertex activity'][pair_type]
            reduced_MCbnbDATAcosmic[pair_type] = sam[np.abs(sam['delta_phi']-180.)<delta_Delta_phi_cut]
        #}
        get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
        get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic);
    #}

    # cut 6: Pt
    for cut_name,Pt_max_cut in zip(['Pt-','Pt+'],[Pt_max_minus,Pt_max_plus]):#{
        cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max_cut
        reduced_MCbnbDATAcosmic = dict()
        for pair_type in pair_types:#{
            sam = reduced_MCbnbDATAcosmicSamples['delta phi'][pair_type]
            reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_Pt']<Pt_max_cut]
        #}
        get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
        get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic);
    #}
   
    # cut 6: \theta_pq
    for cut_name,theta_pq_max_cut in zip(['theta_pq-','theta_pq+'],[theta_pq_max_minus,theta_pq_max_plus]):#{
        cut_label=r'$\theta_{pq}<%.0f^0$ deg.'%theta_pq_max_cut
        reduced_MCbnbDATAcosmic = dict()
        for pair_type in pair_types:#{
            sam = reduced_MCbnbDATAcosmicSamples['delta phi'][pair_type]
            reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_theta_pq']<theta_pq_max_cut]
        #}
        get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
        get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic);
    #}

    return pureff_MCbnbDATAcosmic,pureff_MCbnbDATAcosmic_numbers
# ------------------------------------------------


# ------------------------------------------------
# last edit:
# Oct. 13, 2017
def apply_cuts( PIDa_p_min=12
               , delta_theta_12=60  # deg.
               , delta_Delta_phi=40 # deg.
               , theta_pq_max=25    # deg.
               , Pt_max=0.35        # GeV/c
               , opt_box=(50,100) # [Nwires x Nticks]
               , r_max_RdQ_CC1p0pi = 0.32 # sphere in U,V,Y space, apply a cut only to CC1p0pi
               , min_PE_ClosestFlash = 50
               , max_dYZ_ClosestFlash = 150
               ):
    '''
        return:
        
            pureff_MCbnbDATAcosmic
            pureff_MCbnbDATAcosmic_numbers
            pureff_MCbnbMCcosmic
            pureff_MCbnbMCcosmic_numbers
    '''
    # --- -- --- - -- -- --- --
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        reduced_MCbnbDATAcosmic[pair_type] = MCbnbDATAcosmicSamples[pair_type]
        reduced_MCbnbMCcosmic[pair_type] = MCbnbMCcosmicSamples[pair_type]
    #}
    reduced_MCbnbDATAcosmicSamples['no cut'] = reduced_MCbnbDATAcosmic
    reduced_MCbnbMCcosmicSamples['no cut'] = reduced_MCbnbMCcosmic

    # before cuts
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['no cut'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam
        sam = reduced_MCbnbMCcosmicSamples['no cut'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name = 'no cut', cut_label='no cut', reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name = 'no cut', cut_label='no cut', reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name = 'no cut', cut_label='no cut', reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name = 'no cut', cut_label='no cut', reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)


    # cut 1: PIDa
    cut_name , cut_label = 'PIDa','${PID}_a>%.0f$'%PIDa_p_min
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['no cut'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['PIDa_assigned_proton']>PIDa_p_min]
        sam = reduced_MCbnbMCcosmicSamples['no cut'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam[sam['PIDa_assigned_proton']>PIDa_p_min]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)


    # cut 2: require that the longer track is the one with larger PIDa
    cut_name , cut_label = 'length', '$l_{\\mu}>l_{p}$'
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['PIDa'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['PIDa_long'] < sam['PIDa_short']]
        sam = reduced_MCbnbMCcosmicSamples['PIDa'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam[sam['PIDa_long'] < sam['PIDa_short']]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    
    # cut 3: |\theta_{1,2}-90^0|<60^0$
    cut_name, cut_label='non-collinearity' ,'$|\theta_{1,2}-90^0|<%.0f^0$'%delta_theta_12
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['length'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
        sam = reduced_MCbnbMCcosmicSamples['length'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam[np.abs(sam['theta_12']-90)<delta_theta_12]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)


    # -- - -- -- - -- -- - - -- --- 
    # FLASH A: look at the closest flash to the vertex
    # -- - -- -- - -- -- - - -- --- 
    # cut 4: number of photoelectrons in the closest flash to the vertex (ClosestFlash_TotalPE)
    cut_name, cut_label='closest-flash PE' ,'$N_{PE}^{closest-flash}>%.0f$'%min_PE_ClosestFlash
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['non-collinearity'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['ClosestFlash_TotalPE'] > min_PE_ClosestFlash]
        sam = reduced_MCbnbMCcosmicSamples['non-collinearity'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam[sam['ClosestFlash_TotalPE'] > min_PE_ClosestFlash]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)


    # cut 5: distance of the closest flash to the vertex in YZ plane (ClosestFlash_YZdistance)
    cut_name, cut_label='closest-flash dYZ' ,'$d_{YZ-plane}^{closest-flash}<%.0f$ cm'%max_dYZ_ClosestFlash
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['closest-flash PE'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['ClosestFlash_YZdistance'] < max_dYZ_ClosestFlash]
        sam = reduced_MCbnbMCcosmicSamples['closest-flash PE'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam[sam['ClosestFlash_YZdistance'] < max_dYZ_ClosestFlash]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)

    
    
    # -- - -- -- - -- -- - - -- --- 
    # VERTEX ACTIVITY 
    # -- - -- -- - -- -- - - -- --- 
    # cut 6: vertex activity (RdQ)
    # for the 1mu-1p sample we do not apply a cut, which means that we stop here for 1mu-1p
    # for the CC 1p 0pi we apply a cut:
    # a sphere around RdQ=1 of radius r_min_RdQ_CC1p0pi
    # where r_min_RdQ_CC1p0pi is taken from analysis_notes/RdQ/RdQaroundVertex_cut_selection
    cut_name , cut_label = 'vertex activity' , '$\sqrt{\sum_{p=0,1,2}(R_{\Delta Q}^{p}-1)^2}<%.2f$'%r_max_RdQ_CC1p0pi
    box_str='[%d wires x %d ticks]'%(opt_box[0],opt_box[1])
    Ru = 'RdQaroundVertex[plane 0]'+box_str
    Rv = 'RdQaroundVertex[plane 1]'+box_str
    Ry = 'RdQaroundVertex[plane 2]'+box_str
    
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()    
    for pair_type in pair_types:#{ 
        
        sam = reduced_MCbnbDATAcosmicSamples['non-collinearity'][pair_type]#closest-flash dYZ
        sam = sam[(sam[Ru]==1) | (sam[Rv]==1) | (sam[Ry]==1) 
                  | 
                  (np.sqrt( np.square(sam[Ru]-1) + np.square(sam[Rv]-1) + np.square(sam[Ry]-1) ) <= r_max_RdQ_CC1p0pi) ]
        reduced_MCbnbDATAcosmic[pair_type] = sam
        
        sam = reduced_MCbnbMCcosmicSamples['closest-flash dYZ'][pair_type]
        sam = sam[(sam[Ru]==1) | (sam[Rv]==1) | (sam[Ry]==1) 
                  | 
                  (np.sqrt( np.square(sam[Ru]-1) + np.square(sam[Rv]-1) + np.square(sam[Ry]-1) ) <= r_max_RdQ_CC1p0pi) ]
        reduced_MCbnbMCcosmic[pair_type] = sam
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name=cut_name ,cut_label=cut_label , reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name=cut_name,cut_label=cut_label , reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name=cut_name, cut_label=cut_label, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)


    # cut 5: $\Delta phi$
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['vertex activity'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[np.abs(sam['delta_phi']-180.)<delta_Delta_phi]
        sam = reduced_MCbnbMCcosmicSamples['vertex activity'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam[np.abs(sam['delta_phi']-180.)<delta_Delta_phi]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name = 'delta phi', cut_label='$|\Delta \phi - \pi|<%.0f^0$'%delta_Delta_phi, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name = '\CutDeltaPhi', cut_label='$|\Delta \phi - \pi|<%.0f^0$'%delta_Delta_phi, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name = 'delta phi', cut_label='$|\Delta \phi - \pi|<%.0f^0$'%delta_Delta_phi, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name = '\CutDeltaPhi', cut_label='$|\Delta \phi - \pi|<%.0f^0$'%delta_Delta_phi, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)


    # cut 6: $\theta_{pq}<25$
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['delta phi'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_theta_pq']<theta_pq_max]
        sam = reduced_MCbnbMCcosmicSamples['delta phi'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam[sam['reco_theta_pq']<theta_pq_max]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name ='theta_pq' , cut_label= '$\theta_{pq}<%.0f^0$'%theta_pq_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name ='theta_pq' , cut_label= '$\\theta_{pq}<%.0f^0$'%theta_pq_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name ='theta_pq' , cut_label= '$\theta_{pq}<%.0f^0$'%theta_pq_max, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name ='theta_pq' , cut_label= '$\\theta_{pq}<%.0f^0$'%theta_pq_max, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)


    # modified cut 6: $p_{t}<0.35$
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['delta phi'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_Pt']<Pt_max]
        sam = reduced_MCbnbMCcosmicSamples['delta phi'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam[sam['reco_Pt']<Pt_max]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name='soft Pt', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name='soft Pt', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name='soft Pt', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name='soft Pt', cut_label='$p_{t}<%.2f$ GeV/c'%Pt_max, reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)


    # tight Pt cut for good Ev reconstruction
    reduced_MCbnbDATAcosmic = dict()
    reduced_MCbnbMCcosmic = dict()
    for pair_type in pair_types:#{
        sam = reduced_MCbnbDATAcosmicSamples['theta_pq'][pair_type]
        reduced_MCbnbDATAcosmic[pair_type] = sam[sam['reco_Pt']<0.15]
        sam = reduced_MCbnbMCcosmicSamples['theta_pq'][pair_type]
        reduced_MCbnbMCcosmic[pair_type] = sam[sam['reco_Pt']<0.15]
    #}
    get_pureff_MCbnbDATAcosmic_cut(cut_name ='tight Pt', cut_label= '$p_{t}<0.15$ GeV/c', reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbDATAcosmic_numbers(cut_name ='tight Pt', cut_label= '$p_{t}<0.15$ GeV/c', reduced_MCbnbDATAcosmic = reduced_MCbnbDATAcosmic)
    get_pureff_MCbnbMCcosmic_cut(cut_name ='tight Pt', cut_label= '$p_{t}<0.15$ GeV/c', reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)
    get_pureff_MCbnbMCcosmic_numbers(cut_name ='tight Pt', cut_label= '$p_{t}<0.15$ GeV/c', reduced_MCbnbMCcosmic = reduced_MCbnbMCcosmic)



    return pureff_MCbnbDATAcosmic,pureff_MCbnbDATAcosmic_numbers,pureff_MCbnbMCcosmic,pureff_MCbnbMCcosmic_numbers
# ------------------------------------------------

# ------------------------------------------------
# Sep 9
def load_MCbnbMCcosmicSamples():
    '''
    return: 
           MCbnbMCcosmicPairsFV:  pandas.DataFrame() of prodgenie_bnb_nu_cosmic_uboone_mcc8.2_reco2_vertices in FV
           samples: dict() of pairs broken into pair_types
    '''
    global MCbnbMCcosmicSamples
    pairs = pd.read_csv(vertices_files_path+'prodgenie_bnb_nu_cosmic_uboone_mcc8.2_reco2_vertices.csv')
    MCbnbMCcosmicPairsFV = sample_in_FV(pairs)
    print len(pairs),'pairs from MC-BNB + cosmic MC overlay'
    print len(MCbnbMCcosmicPairsFV),'pairs in FV'
    for pair_type in pair_types:#{
        MCbnbMCcosmicSamples[pair_type] = MCbnbMCcosmicPairsFV[MCbnbMCcosmicPairsFV[pair_type]==True]
        Ntype = len(MCbnbMCcosmicSamples[pair_type])
        if pair_type=='CC 1p 0pi': print_line()
        print Ntype,'are '+pair_type+', %.1f'%(100.*float(Ntype)/len(MCbnbMCcosmicPairsFV))+'%'
    #}
    return MCbnbMCcosmicPairsFV, MCbnbMCcosmicSamples
# ------------------------------------------------


# ------------------------------------------------
# Aug 30
def load_MCbnbDATAcosmicSamples():
    '''
        return:
        MCbnbDATAcosmicSamples:  pandas.DataFrame() of prodgenie_bnb_nu_uboone_overlay_mcc8_reco2_vertices in FV
        samples: dict() of pairs broken into pair_types
        '''
    global MCbnbDATAcosmicSamples
    # old overlay: prodgenie_bnb_nu_uboone_overlay_mcc8_reco2
    #pairs = pd.read_csv(vertices_files_path+'prodgenie_bnb_nu_uboone_overlay_mcc8_reco2_vertices.csv')
    # new overlay: ccqe_ana_MCBNBCosmicDATA_2018_01_30
    pairs = pd.read_csv(vertices_files_path+'ecohen_physical_files_adi_prodgenie_bnb_nu_uboone_overlay_cosmic_data_100K_reco2_2018_02_10_vertices.csv')
    MCbnbDATAcosmicPairsFV = sample_in_FV(pairs)
    print len(pairs),'pairs from MC-BNB + cosmic DATA overlay'
    print len(MCbnbDATAcosmicPairsFV),'pairs in FV'
    for pair_type in pair_types:#{
        MCbnbDATAcosmicSamples[pair_type] = MCbnbDATAcosmicPairsFV[MCbnbDATAcosmicPairsFV[pair_type]==True]
        Ntype = len(MCbnbDATAcosmicSamples[pair_type])
        if pair_type=='CC 1p 0pi': print_line()
        print Ntype,'are '+pair_type+', %.1f'%(100.*float(Ntype)/len(MCbnbDATAcosmicPairsFV))+'%'
    #}
    return MCbnbDATAcosmicPairsFV, MCbnbDATAcosmicSamples
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
# Oct. 13, 2017
def get_pureff_MCbnbMCcosmic_cut(cut_name = 'PIDa', cut_label=None , reduced_MCbnbMCcosmic = dict()):
    '''
        return
        eff (mu-p) , pur (mu-p), eff (CC 1p 0pi) , pur (CC 1p 0pi)
        '''
    
    global pureff_MCbnbMCcosmic
    eff = dict()
    pur = dict()
    Nreduced_MCbnbMCcosmic , freduced_MCbnbMCcosmic = get_Nreduced_MCbnbMCcosmic(reduced_MCbnbMCcosmic=reduced_MCbnbMCcosmic)
    Ntot = (Nreduced_MCbnbMCcosmic['1mu-1p']+Nreduced_MCbnbMCcosmic['cosmic']+Nreduced_MCbnbMCcosmic['other pairs'])
    
    eff['1mu-1p'] = freduced_MCbnbMCcosmic['1mu-1p']
    pur['1mu-1p'] = 100.*Nreduced_MCbnbMCcosmic['1mu-1p']/Ntot if Ntot>0 else 0
    
    eff['CC 1p 0pi'] = freduced_MCbnbMCcosmic['CC 1p 0pi']
    pur['CC 1p 0pi'] = 100.*Nreduced_MCbnbMCcosmic['CC 1p 0pi']/Ntot if Ntot>0 else 0
    
    pureff_MCbnbMCcosmic_cut = pd.DataFrame({'label':cut_label
                                              ,'$\mu p$ eff.':'%.1f'%eff['1mu-1p']+'%'
                                              ,'$\mu p$ pur.':'%.1f'%pur['1mu-1p']+'%'
                                              ,'CC$0\pi 1 p$ eff.':'%.1f'%freduced_MCbnbMCcosmic['CC 1p 0pi']+'%'
                                              ,'CC$0\pi 1 p$ pur.':'%.1f'%(100.*Nreduced_MCbnbMCcosmic['CC 1p 0pi']/Ntot if Ntot>0 else 0)+'%'}
                                              , index=[cut_name]
                                              )
    for pair_type in pair_types: pureff_MCbnbMCcosmic_cut[pair_type] = '%.1f'%freduced_MCbnbMCcosmic[pair_type]+'%' +' (%.0f)'%Nreduced_MCbnbMCcosmic[pair_type]
    pureff_MCbnbMCcosmic = pureff_MCbnbMCcosmic.append(pureff_MCbnbMCcosmic_cut)
    reduced_MCbnbMCcosmicSamples[cut_name] = reduced_MCbnbMCcosmic
    Ntot = Nreduced_MCbnbMCcosmic['1mu-1p']+Nreduced_MCbnbMCcosmic['cosmic']+Nreduced_MCbnbMCcosmic['other pairs']
    return freduced_MCbnbMCcosmic['1mu-1p'],(100.*Nreduced_MCbnbMCcosmic['1mu-1p']/Ntot if Ntot>0 else 0),freduced_MCbnbMCcosmic['CC 1p 0pi'],(100.*Nreduced_MCbnbMCcosmic['CC 1p 0pi']/Ntot if Ntot>0 else 0)
# ------------------------------------------------

  
      
# ------------------------------------------------
# Oct. 13, 2017
def get_pureff_MCbnbDATAcosmic_cut(cut_name = 'PIDa', cut_label=None , reduced_MCbnbDATAcosmic = dict()):
    ''' 
        return
        eff (mu-p) , pur (mu-p), eff (CC 1p 0pi) , pur (CC 1p 0pi)
    '''
    
    global pureff_MCbnbDATAcosmic
    eff = dict()
    pur = dict()
    Nreduced_MCbnbDATAcosmic , freduced_MCbnbDATAcosmic = get_Nreduced_MCbnbDATAcosmic(reduced_MCbnbDATAcosmic=reduced_MCbnbDATAcosmic)
    Ntot = (Nreduced_MCbnbDATAcosmic['1mu-1p']+Nreduced_MCbnbDATAcosmic['cosmic']+Nreduced_MCbnbDATAcosmic['other pairs'])
    
    eff['1mu-1p'] = freduced_MCbnbDATAcosmic['1mu-1p']
    pur['1mu-1p'] = 100.*Nreduced_MCbnbDATAcosmic['1mu-1p']/Ntot if Ntot>0 else 0
    
    eff['CC 1p 0pi'] = freduced_MCbnbDATAcosmic['CC 1p 0pi']
    pur['CC 1p 0pi'] = 100.*Nreduced_MCbnbDATAcosmic['CC 1p 0pi']/Ntot if Ntot>0 else 0
    
    pureff_MCbnbDATAcosmic_cut = pd.DataFrame({'label':cut_label
                               ,'$\mu p$ eff.':'%.1f'%eff['1mu-1p']+'%'
                               ,'$\mu p$ pur.':'%.1f'%pur['1mu-1p']+'%'
                               ,'CC$0\pi 1 p$ eff.':'%.1f'%freduced_MCbnbDATAcosmic['CC 1p 0pi']+'%'
                               ,'CC$0\pi 1 p$ pur.':'%.1f'%(100.*Nreduced_MCbnbDATAcosmic['CC 1p 0pi']/Ntot if Ntot>0 else 0)+'%'}
                               , index=[cut_name]
                              )
    for pair_type in pair_types: pureff_MCbnbDATAcosmic_cut[pair_type] = '%.1f'%freduced_MCbnbDATAcosmic[pair_type]+'%' +' (%.0f)'%Nreduced_MCbnbDATAcosmic[pair_type]
    pureff_MCbnbDATAcosmic = pureff_MCbnbDATAcosmic.append(pureff_MCbnbDATAcosmic_cut)
    reduced_MCbnbDATAcosmicSamples[cut_name] = reduced_MCbnbDATAcosmic  
    Ntot = Nreduced_MCbnbDATAcosmic['1mu-1p']+Nreduced_MCbnbDATAcosmic['cosmic']+Nreduced_MCbnbDATAcosmic['other pairs']
    return freduced_MCbnbDATAcosmic['1mu-1p'],(100.*Nreduced_MCbnbDATAcosmic['1mu-1p']/Ntot if Ntot>0 else 0),freduced_MCbnbDATAcosmic['CC 1p 0pi'],(100.*Nreduced_MCbnbDATAcosmic['CC 1p 0pi']/Ntot if Ntot>0 else 0)
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
                                                      '\CCIpOpi':Nreduced_MCbnbMCcosmic['CC 1p 0pi'],
                                                      'eff \mup':freduced_MCbnbMCcosmic['1mu-1p'],
                                                      'eff \CCIpOpi':freduced_MCbnbMCcosmic['CC 1p 0pi'],
                                                      'pur \mup':float(100*Nreduced_MCbnbMCcosmic['1mu-1p'])/Ntot if Ntot>0 else 0,
                                                      'pur \CCIpOpi':float(100*Nreduced_MCbnbMCcosmic['CC 1p 0pi'])/Ntot if Ntot>0 else 0}
                                                      , index=[cut_name]
                                                      )
    pureff_MCbnbMCcosmic_numbers = pureff_MCbnbMCcosmic_numbers.append(pureff_MCbnbMCcosmic_numbers_cut)
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
                                        '\CCIpOpi':Nreduced_MCbnbDATAcosmic['CC 1p 0pi'],
                                        'eff \mup':freduced_MCbnbDATAcosmic['1mu-1p'],                               
                                        'eff \CCIpOpi':freduced_MCbnbDATAcosmic['CC 1p 0pi'],
                                        'pur \mup':float(100*Nreduced_MCbnbDATAcosmic['1mu-1p'])/Ntot if Ntot>0 else 0,
                                        'pur \CCIpOpi':float(100*Nreduced_MCbnbDATAcosmic['CC 1p 0pi'])/Ntot if Ntot>0 else 0}
                                       , index=[cut_name]
                                      )
    pureff_MCbnbDATAcosmic_numbers = pureff_MCbnbDATAcosmic_numbers.append(pureff_MCbnbDATAcosmic_numbers_cut)
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
        sample = reduced_MCbnbDATAcosmicSamples[cut_name][pair_type]
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
def plot_cut_samples (reduced_samples=None,
                      reduced_cut_name='PIDa',markers_size=5,
                      cut_name='maximal distance between tracks',mul=1,
                      cut_var ='distance',
                      cut_type= 'max',
                      x_label = 'maximal tracks distance [cm]', y_label='% of sample',
                      xcenter=0,figsize=(12,8),fontsize=25,
                      xmin=0.1, xmax=10 , Nbins=10, do_add_legend=True, legend_loc='bbox',legend_fontsize=25,
                      ticks_color='black'):
    fig,ax=plt.subplots(figsize=figsize)
    for i,(pair_type,label,cmap,color) in enumerate(zip(pair_types,MClabels,MCcmaps,MCcolors)):
        sample = reduced_MCbnbDATAcosmicSamples[reduced_cut_name][pair_type]
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




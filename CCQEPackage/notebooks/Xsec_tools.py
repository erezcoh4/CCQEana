import sys, pandas as pd, numpy as np, matplotlib.pyplot as plt
sys.path.insert(0, '/Users/erezcohen/larlite/UserDev/mySoftware/MySoftwarePackage/mac')
from plot_tools import *

import sys; sys.path.insert(0, '../..');
from ccqe_notebook_tools import *
from mupClassification_notebooks import *
from onbeam_offbeam_notebooks import *

from mpl_toolkits.mplot3d import Axes3D


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
    other-pairs
    cosmic in the overlay (which are correlated with the beam)
    cosmic not in the overlay (off-beam data)
    '''
flux = 3.601e10 # cm^-2
flux_err = 0
Ntargets = 1.25e31
Ntargets_err = 0

Limits=dict({
            'Pmu':(0.2,1.405)
            ,'cos(theta(mu))':(-0.5,0.95)
            ,'phi(mu)':(-180,180)
            ,'Pp':(0.3,1.0)
            ,'cos(theta(p))':(0.25,0.9)
            ,'phi(p)':(-180,180)
            })
NBins=5
Bins = dict()
for key in Limits.keys(): Bins[key] = np.linspace(Limits[key][0],Limits[key][1],NBins+1)
Centers = dict()
for key in Limits.keys(): Centers[key] = 0.5*(Bins[key][1:] + Bins[key][:-1])


vlabels = dict({
               'Pmu':r"p_{\mu}"
               ,'theta(mu)':r"\theta_{\mu}"
               ,'phi(mu)':r"\phi_{\mu}"
               ,'cos(theta(mu))':r"\cos(\theta_{\mu})"
               ,'Pp':r"p_{p}"
               ,'theta(p)':r"\theta_{p}"
               ,'cos(theta(p))':r"\cos(\theta_{p})"
               ,'phi(p)':r"\phi_{p}"
               })

Vlabels = dict()
for key in vlabels.keys(): Vlabels[key] = r"$"+vlabels[key]+"$"

Units = dict({
             'Pmu':r"GeV/c"
             ,'theta(mu)':r"deg."
             ,'cos(theta(mu))':None
             ,'phi(mu)':r"deg."
             ,'Pp':r"GeV/c"
             ,'theta(p)':r"deg."
             ,'cos(theta(p))':None
             ,'phi(p)':r"deg."
             })

Colors = dict({
              'overlay':'forestgreen'
              ,'CC1p':'blue'
              ,'beam off':'black'
              ,'beam on':'tomato'
              })

Xsec_path = '/Users/erezcohen/Desktop/uBoone/CCQEanalysis/Xsec/'
Paths = dict({'selected events':Xsec_path+'selected_events/'
             ,'migration maps':Xsec_path+'migration_maps/'
             ,'background maps':Xsec_path+'background_maps/'
             ,'efficiency maps':Xsec_path+'efficeincy_maps/'
             ,'1d Xsec':Xsec_path+'1d_Xsec/'})


#
#
## ----------------------------------------------------------
## Sep-08, 2018
#def compute_Xsec_in_3d(beam_on=None,beam_off=None
#                       ,generated_CC1p=None,selected_CC1p=None,overlay=None
#                       ,NBins=5):
#    # the return is a dictionary of results
#    for key in Limits.keys(): Bins[key] = np.linspace(Limits[key][0],Limits[key][1],NBins+1)
#    global bins1,bins2,bins3,N1,N2,N3
#    bins1,bins2,bins3 = Bins['Pmu'], Bins['cos(theta(mu))'] , Bins['phi(mu)']
#    N1,N2,N3 = len(bins1)-1,len(bins2)-1,len(bins3)-1
#    
#    h = dict()
#    keys = ['on','off scaled','generated','CC1p','CC1p scaled','B','eff','eff err'
#            ,'Xsec','Xsec err','mc-Xsec','mc-Xsec err','generated-Xsec','generated-Xsec err']
#    for key in keys: h[key] = np.zeros((N1,N2,N3))
#    N = dict()
#    for i_P in range(N1):#{
#        Pmin,Pmax = bins1[i_P],bins1[i_P+1]
#        P_bin_width = Pmax - Pmin
#                        
#        for i_cos_theta in range(N2):#{
#            cos_theta_min,cos_theta_max = bins2[i_cos_theta],bins2[i_cos_theta+1]
#            cos_theta_bin_width = cos_theta_max - cos_theta_min
#                    
#            for i_phi in range(N3):#{
#                phi_min,phi_max = bins3[i_phi],bins3[i_phi+1]
#                phi_bin_width = phi_max - phi_min
#                
#                bin_width = P_bin_width * cos_theta_bin_width * phi_bin_width
#                                                    
#
#    beam_on_in_bin = sam_in_3d_bin(beam_on,
#                               'reco_Pmu_mcs',Pmin,Pmax,
#                               'reco_Pmu_cos_theta',cos_theta_min,cos_theta_max,
#                               'reco_Pmu_mcs_phi',phi_min,phi_max)
#    N['on'] = len(beam_on_in_bin)
#        h['on'][i_P][i_cos_theta][i_phi] = N['on']
#            
#            Xsec_in_bin , Xsec_err_in_bin = 0 , 0
#                mc_Xsec_in_bin , mc_Xsec_err_in_bin = 0 , 0
#                    gen_Xsec_in_bin , gen_Xsec_err_in_bin = 0 , 0
#                        
#                        N['off'] = len(sam_in_3d_bin(beam_off,
#                                                     'reco_Pmu_mcs',Pmin,Pmax,
#                                                     'reco_Pmu_cos_theta',cos_theta_min,cos_theta_max,
#                                                     'reco_Pmu_mcs_phi',phi_min,phi_max))
#                            N['off scaled'] = N['off']*OffBeam_scaling
#                                h['off scaled'][i_P][i_cos_theta][i_phi] = N['off scaled']
#                                    
#                                    N['generated'] = len(sam_in_3d_bin(generated_CC1p,
#                                                                       'truth_Pmu',Pmin,Pmax,
#                                                                       'truth_Pmu_cos_theta',cos_theta_min,cos_theta_max,
#                                                                       'truth_Pmu_phi',phi_min,phi_max) )
#                                        N['generated scaled'] = N['generated']*Nevents['f(POT)']
#                                            h['generated'][i_P][i_cos_theta][i_phi] = N['generated']
#                                                
#                                                N['CC1p'] = len(sam_in_3d_bin(selected_CC1p,
#                                                                              'reco_Pmu_mcs',Pmin,Pmax,
#                                                                              'reco_Pmu_cos_theta',cos_theta_min,cos_theta_max,
#                                                                              'reco_Pmu_mcs_phi',phi_min,phi_max) )
#                                                    h['CC1p'][i_P][i_cos_theta][i_phi] = N['CC1p']
#                                                        
#                                                        N['CC1p scaled'] = N['CC1p']*Nevents['f(POT)']
#                                                            h['CC1p scaled'][i_P][i_cos_theta][i_phi] = N['CC1p scaled']
#                                                                
#                                                                
#                                                                N['ovrelay'] = len(sam_in_3d_bin(overlay,
#                                                                                                 'reco_Pmu_mcs',Pmin,Pmax,
#                                                                                                 'reco_Pmu_cos_theta',cos_theta_min,cos_theta_max,
#                                                                                                 'reco_Pmu_mcs_phi',phi_min,phi_max) )
#                                                                    N['ovrelay scaled'] = N['ovrelay']*Nevents['f(POT)']
#                                                                        
#                                                                        B, B_err = N['ovrelay scaled'] - N['CC1p scaled'] , np.sqrt(N['ovrelay scaled'] - N['CC1p scaled'])
#                                                                            h['B'][i_P][i_cos_theta][i_phi] = B
#                                                                                
#                                                                                eff, eff_err = get_eff(Ngen=N['generated'] , Nsel=N['CC1p'])
#                                                                                    h['eff'][i_P][i_cos_theta][i_phi] = eff
#                                                                                        h['eff err'][i_P][i_cos_theta][i_phi] = eff_err
#                                                                                            
#                                                                                            Xsec_in_bin,Xsec_err_in_bin = compute_Xsec(Non=N['on'], Noff=N['off'], B=B, eff=eff,
#                                                                                                                                       bin_width = bin_width,
#                                                                                                                                       Non_err = np.sqrt(N['on']),
#                                                                                                                                       Noff_err= np.sqrt(N['off']),
#                                                                                                                                       B_err   = B_err, eff_err = eff_err)
#                                                                                                
#                                                                                                mc_Xsec_in_bin,mc_Xsec_err_in_bin = compute_Xsec(Non=N['CC1p scaled'], eff=eff,
#                                                                                                                                                 bin_width = bin_width,
#                                                                                                                                                 Non_err = np.sqrt(N['CC1p'])*Nevents['f(POT)'],
#                                                                                                                                                 eff_err = eff_err)
#                                                                                                    
#                                                                                                    gen_Xsec_in_bin,gen_Xsec_err_in_bin = compute_Xsec(Non=N['generated scaled'], eff=1,
#                                                                                                                                                       bin_width = bin_width,
#                                                                                                                                                       Non_err = np.sqrt(N['generated'])*Nevents['f(POT)'])
#                                                                                                        h['Xsec'][i_P][i_cos_theta][i_phi] = Xsec_in_bin
#                                                                                                            h['Xsec err'][i_P][i_cos_theta][i_phi] = Xsec_err_in_bin
#                                                                                                                h['mc-Xsec'][i_P][i_cos_theta][i_phi] = mc_Xsec_in_bin
#                                                                                                                    h['mc-Xsec err'][i_P][i_cos_theta][i_phi] = mc_Xsec_err_in_bin
#                                                                                                                        h['generated-Xsec'][i_P][i_cos_theta][i_phi] = gen_Xsec_in_bin
#                                                                                                                            h['generated-Xsec err'][i_P][i_cos_theta][i_phi] = gen_Xsec_err_in_bin
## if N['on']>0:#{#} if N['on'] > 0
##} i_Pmu_phi
##} i_Pmu_cos_theta
##} i_Pmu
#print 'done.'
#    return h
## ----------------------------------------------------------
## Sep-08, 2018
#



# ----------------------------------------------------------
# Sep-04, 2018
def get_labels(observable=''):
    bins=Bins[observable]; vlabel=vlabels[observable]; Vlabel=Vlabels[observable]; units=Units[observable]
    xlabel=Vlabel+' ['+units+']' if units is not None else Vlabel
    mid = 0.5*(bins[1:]+bins[:-1]); bin_width=(mid[1]-mid[0])
    return bins,mid,bin_width,vlabel,xlabel,units
# ----------------------------------------------------------

# ----------------------------------------------------------
# Sep-02, 2018
def sam_in_3d_bin(sam,
                  Pvar,Pmin,Pmax,
                  cos_theta_var,cos_theta_min,cos_theta_max,
                  phi_var,phi_min,phi_max):
    
    return sam[(Pmin<=sam[Pvar])&(sam[Pvar] < Pmax)
               &
               (cos_theta_min<=sam[cos_theta_var])&(sam[cos_theta_var] < cos_theta_max)
               &
               (phi_min <= 180./np.pi*sam[phi_var])&(180./np.pi*sam[phi_var] < phi_max)
               ]
# ----------------------------------------------------------


# ----------------------------------------------------------
# Aug-29, 2018
def compute_Xsec(Non=1, Noff=0, B=0, eff=1, bin_width=1,
                 Non_err=1, Noff_err=0, B_err=0, eff_err=0, eff_cutoff=0.05):
    '''
        input:
        ------
        Non         number of beam on events
        Noff        number of beam off events
        B           background estimation from overlay
        eff         efficiency (built-in cut-off <eff_cutoff=5%>: if eff<eff_cutoff return 0 +/- 0 )
        [err]       uncertainties
        
        return:
        ------
        Xsec        cross-section in units of (10^-39) cm2 / bin_units
        Xsec_err    cross-section uncertainty in units of (10^-39) cm2 / bin_units
        '''
    
    if eff<eff_cutoff:  return 0,0
    

    num = Non - Noff - B
    den = eff * Ntargets * flux * bin_width
    Xsec = np.max( [0, num/den] )

    num_err = np.sqrt( np.square(Non_err) + np.square(Noff_err) + np.square(B_err) )
    den_err = den * np.sqrt( np.square(eff_err/eff) + np.square(Ntargets_err/Ntargets) + np.square(flux_err/flux) )
    Xsec_err = Xsec * np.max( [0 , np.sqrt(  (np.square(num_err/num) if num>0 else 0) + (np.square(den_err/den) if den>0 else 0) ) ] )
    
    return Xsec*1e39, Xsec_err*1e39
# ----------------------------------------------------------



# ----------------------------------------------------------
# Aug-27, 2018
def full_chain_Xsec_diff_1d(selected_beam_off=None,selected_beam_on=None,selected_overlay_concat=None,selected_CC1p=None
                            ,observable='Pmu',truth_var='truth_Pmu',recovar='reco_Pmu_mcs',smearedvar='',mul=1,ax=None
                            ,mc_scale_factor=1,do_draw_all=False,debug=0,extra_name=''):
    bins=Bins[observable]; vlabel=vlabels[observable]; Vlabel=Vlabels[observable]; units=Units[observable]
    xlabel=Vlabel+' ['+units+']' if units is not None else Vlabel
    
    # (1) background subtraction
    subtrsact_bkg_1d(selected_beam_off=selected_beam_off,selected_beam_on=selected_beam_on
                     ,selected_overlay_concat=selected_overlay_concat,selected_CC1p=selected_CC1p
                     ,bins=bins,xlabel=xlabel,xvar=recovar,debug=debug,do_draw=do_draw_all,extra_name=extra_name)
                
    # (2) efficiency:
    compute_effiency(genie_CC1p=genie_CC1p,selected_CC1p=selected_CC1p
                                      ,bins=bins,xvar=smearedvar,xlabel=xlabel,ylabel=r'$\bar{\epsilon}$',do_draw=do_draw_all,mul=mul,debug=debug,extra_name=extra_name)
                     
    # (3) cross-section
    Xsec_diff_1d(observable=observable,recovar=recovar,smearedvar=smearedvar,ax=ax,mc_scale_factor=mc_scale_factor,debug=0,extra_name=extra_name)
# ----------------------------------------------------------


# ----------------------------------------------------------
# Aug-27, 2018
def Xsec_ratio_1d(observable='Pmu',recovar='reco_Pmu_mcs',smearedvar='',debug=0,ax=None,extra_name='mc_mc',extra_label='mc/mc overlay'):#{
    bins=Bins[observable]; vlabel=vlabels[observable]; Vlabel=Vlabels[observable]; units=Units[observable]
    xlabel=Vlabel+' ['+units+']' if units is not None else Vlabel
    mid = 0.5*(bins[1:]+bins[:-1]); bin_width=0.5*(mid[1]-mid[0])
    
    h,herr = dict(),dict()
    
    h['nominal Xsec'] = np.loadtxt(Paths['1d Xsec'] + "Xsec_%s_%d_bins.csv"%(smearedvar,len(bins)-1), delimiter=",")
    herr['nominal Xsec'] = np.loadtxt(Paths['1d Xsec'] + "Xec_err_%s_%d_bins.csv"%(smearedvar,len(bins)-1), delimiter=",")
    if debug: print 'read',"Xsec_%s_%d_bins.csv"%(smearedvar,len(bins)-1)
    h['modified Xsec'] = np.loadtxt(Paths['1d Xsec'] + extra_name + "Xsec_%s_%d_bins"%(smearedvar,len(bins)-1)+".csv", delimiter=",")
    herr['modified Xsec'] = np.loadtxt(Paths['1d Xsec'] + extra_name+ "Xec_err_%s_%d_bins"%(smearedvar,len(bins)-1)+".csv", delimiter=",")
    if debug: print 'read', extra_name+"Xsec_%s_%d_bins"%(smearedvar,len(bins)-1)+".csv"
    
    h['nominal mc'] = np.loadtxt(Paths['1d Xsec'] + "mc_cc1p_Xsec_%s_%d_bins.csv"%(smearedvar,len(bins)-1), delimiter=",")
    herr['nominal mc'] = np.loadtxt(Paths['1d Xsec'] + "mc_cc1p_Xsec_err_%s_%d_bins.csv"%(smearedvar,len(bins)-1), delimiter=",")
    if debug: print 'read',"mc_cc1p_Xsec_%s_%d_bins.csv"%(smearedvar,len(bins)-1)
    h['modified mc'] = np.loadtxt(Paths['1d Xsec'] +  extra_name+"mc_cc1p_Xsec_%s_%d_bins"%(smearedvar,len(bins)-1)+".csv", delimiter=",")
    herr['modified mc'] = np.loadtxt(Paths['1d Xsec'] +  extra_name+"mc_cc1p_Xsec_err_%s_%d_bins"%(smearedvar,len(bins)-1)+".csv", delimiter=",")
    if debug: print 'read', extra_name+"mc_cc1p_Xsec_%s_%d_bins"%(smearedvar,len(bins)-1)+".csv"
    
    print 'done reading Xsec.'
    if ax is None:
        fig,ax = plt.subplots(figsize=(9.708,6))
    for slabel,color in zip(['mc','Xsec'],[Colors['CC1p'],Colors['beam on']]):
        h['ratio '+slabel] = h['modified '+slabel]/h['nominal '+slabel]
        herr['ratio '+slabel] = h['ratio '+slabel]*np.sqrt( np.square(herr['nominal '+slabel]/h['nominal '+slabel])
                                                           + np.square(herr['modified '+slabel]/h['modified '+slabel]) )
    ax.errorbar( x=mid , xerr=bin_width, y=h['ratio '+slabel], yerr=herr['ratio '+slabel]
                , fmt='o', markersize=10 , color=color, capsize=1, capthick=3, label='data')
    set_axes(ax,xlabel,y_label=r'Xsec. ratio',do_add_legend=False,title=extra_label)
    if debug:#{
        print "h[nominal "+slabel+"]:",h['nominal '+slabel]
        print "herr[nominal "+slabel+"]:",herr['nominal '+slabel]
        print "h[modified "+slabel+"]:",h['modified '+slabel]
        print "herr[modified "+slabel+"]:",herr['modified '+slabel]
        print "h[ratio "+slabel+"]:",h['ratio '+slabel]
        print "herr[ratio "+slabel+"]:",herr['ratio '+slabel]
    #}
    print 'done ploting Xsec. ratio'
#}
# ----------------------------------------------------------



# ----------------------------------------------------------
# Aug-27, 2018
def build_migration_maps(selected_CC1p=None):
    migration_maps=dict()
    migration_maps['Pmu']=build_migration_map(sam=selected_CC1p,bins=Bins['Pmu'],xvar='truth_Pmu',yvar='reco_Pmu_mcs')
    migration_maps['Pp']=build_migration_map(sam=selected_CC1p,bins=Bins['Pp'],xvar='truth_Pp',yvar='reco_Pp')
    migration_maps['Pmu_cos_theta']=build_migration_map(sam=selected_CC1p,bins=Bins['cos(theta(mu))'],xvar='truth_Pmu_cos_theta',yvar='reco_Pmu_cos_theta')
    migration_maps['Pp_cos_theta']=build_migration_map(sam=selected_CC1p,bins=Bins['cos(theta(p))'],xvar='truth_Pp_cos_theta',yvar='reco_Pp_cos_theta')
    return migration_maps
# ----------------------------------------------------------



# ----------------------------------------------------------
# Aug-27, 2018
def smear_MC_gen_sel(migration_maps=None,selected_CC1p=None,genie_CC1p=None,extra_name=''):
    smear_MC(in_sample=selected_CC1p,migration_maps=migration_maps,name='selected_CC1p',extra_name=extra_name)
    smear_MC(in_sample=genie_CC1p,migration_maps=migration_maps,name='selected_genie_CC1p',extra_name=extra_name)
# ----------------------------------------------------------



# ----------------------------------------------------------
# Aug-27, 2018
def Xsec_diff_1d(observable='Pmu',recovar='reco_Pmu_mcs',smearedvar='',ax=None
                 ,mc_scale_factor=1,debug=0,extra_name=''):#{
    
    bins=Bins[observable]; vlabel=vlabels[observable]; Vlabel=Vlabels[observable]; units=Units[observable]
    xlabel=Vlabel+' ['+units+']' if units is not None else Vlabel
    
    if smearedvar is '': smearedvar='smeared_'+observable
    mul=180./np.pi if 'theta' in observable else 1
    mid = 0.5*(bins[1:]+bins[:-1]); bin_width=0.5*(mid[1]-mid[0])
    
    h,herr = dict(),dict()
    h['N-B'] = np.loadtxt(Paths['background maps'] + extra_name + "beam_on_bkg_sbtrctd_%s_%d_bins.csv"%(recovar,len(bins)-1), delimiter=",")
    herr['N-B'] = np.loadtxt(Paths['background maps'] + extra_name + "beam_on_bkg_sbtrctd_err_%s_%d_bins.csv"%(recovar,len(bins)-1), delimiter=",")
    if debug: print 'read ', extra_name + "beam_on_bkg_sbtrctd_%s_%d_bins.csv"%(recovar,len(bins)-1)
    
    h['mc'] = np.loadtxt(Paths['background maps'] + extra_name + "mc_cc1p_%s_%d_bins.csv"%(recovar,len(bins)-1), delimiter=",")
    herr['mc'] = np.loadtxt(Paths['background maps'] + extra_name + "mc_cc1p_err_%s_%d_bins.csv"%(recovar,len(bins)-1), delimiter=",")
    if debug: print 'read ',"mc_cc1p_%s_%d_bins.csv"%(recovar,len(bins)-1)
    
    h['eff'] = np.loadtxt(Paths['efficiency maps'] + extra_name + "eff_%s_%d_bins.csv"%(smearedvar,len(bins)-1), delimiter=",")
    herr['eff'] = np.loadtxt(Paths['efficiency maps'] + extra_name + "eff_err_%s_%d_bins.csv"%(smearedvar,len(bins)-1), delimiter=",")
    if debug: print 'read efficiency',extra_name + "eff_%s_%d_bins.csv"%(smearedvar,len(bins)-1)
    
    h['Xsec'] = h['N-B']/(h['eff']*Ntargets*flux*bin_width)
    herr['Xsec'] = np.sqrt( np.square(herr['N-B']/(h['eff']*Ntargets*flux*bin_width))
                           +np.square(h['N-B']*herr['eff']/(h['eff']*h['eff']*Ntargets*flux*bin_width))
                           +np.square(h['N-B']*Ntargets_err/(h['eff']*Ntargets*Ntargets*flux*bin_width))
                           +np.square(h['N-B']*flux_err/(h['eff']*Ntargets*flux*flux*bin_width)))
    h['Xsec'],herr['Xsec'] = h['Xsec']*1e39,herr['Xsec']*1e39
                           

    h['mc'] = h['mc']/(h['eff']*Ntargets*flux*bin_width)
    
    herr['mc'] = np.sqrt( np.square(herr['mc']/(h['eff']*Ntargets*flux*bin_width))
                                                +np.square(h['mc']*herr['eff']/(h['eff']*h['eff']*Ntargets*flux*bin_width))
                                                +np.square(h['mc']*Ntargets_err/(h['eff']*Ntargets*Ntargets*flux*bin_width))
                                                +np.square(h['mc']*flux_err/(h['eff']*Ntargets*flux*flux*bin_width)))
    h['mc'],herr['mc'] = h['mc']*1e39*mc_scale_factor,herr['mc']*1e39*mc_scale_factor
    
    # draw the cross-section
    if ax is None:#{
        fig,ax = plt.subplots(figsize=(9.708,6))
    #}
    ax.bar( x=mid , height=2*herr['mc'], bottom=h['mc']-herr['mc'], width=2*bin_width, color=Colors['CC1p'], label=r'mc ($\Delta$stat.)')
    ax.errorbar( x=mid , xerr=bin_width, y=h['Xsec'], yerr=herr['Xsec'] , fmt='o', markersize=10
                    , color=Colors['beam on'], capsize=1, capthick=3, label='data')
    set_axes(ax,xlabel
             ,y_label=(r'$\frac{d\sigma}{d'+vlabel+'}$'
                       +r' $\left[\times 10^{-39}\right]$ '
                       +(r'$\left[\frac{cm^{2}}{(%s)}\right]$'%units
                         if units is not None
                         else r'[cm$^{2}$]'))
             ,do_add_legend=True
             ,title='' if mc_scale_factor == 1 else 'mc scaled by %.3f'%mc_scale_factor)
    if mc_scale_factor == 1:#{
        np.savetxt(Paths['1d Xsec'] + extra_name + "Xsec_%s_%d_bins.csv"%(smearedvar,len(bins)-1), h['Xsec'], delimiter=",")
        np.savetxt(Paths['1d Xsec'] + extra_name + "Xec_err_%s_%d_bins.csv"%(smearedvar,len(bins)-1), herr['Xsec'], delimiter=",")
        if debug: print 'saved',extra_name + "Xsec_%s_%d_bins.csv"%(smearedvar,len(bins)-1)

        np.savetxt(Paths['1d Xsec'] + extra_name + "mc_cc1p_Xsec_%s_%d_bins.csv"%(smearedvar,len(bins)-1), h['mc'], delimiter=",")
        np.savetxt(Paths['1d Xsec'] + extra_name + "mc_cc1p_Xsec_err_%s_%d_bins.csv"%(smearedvar,len(bins)-1), herr['mc'], delimiter=",")
        if debug: print 'saved',extra_name + "mc_cc1p_Xsec_%s_%d_bins.csv"%(smearedvar,len(bins)-1)
    #}
    print 'done computing Xsec.'
#}
# ----------------------------------------------------------


# ----------------------------------------------------------
# Aug-27, 2018
def find_bin( x , bins ):#{
    for i in range(len(bins)-1):#{
        if bins[i]<x and x<bins[i+1]:#{
            return i
        #}
    #}
    return 0
#}
# ----------------------------------------------------------



# ----------------------------------------------------------
# Aug-27, 2018
def subtrsact_bkg_1d(selected_beam_off=None,selected_beam_on=None,selected_overlay_concat=None,selected_CC1p=None
                     ,bins=Bins['Pmu'],xlabel='',xvar='reco_Pmu_mcs',mul=1,do_draw=True,debug=0
                     ,extra_name=''):#{
    mid = 0.5*(bins[1:]+bins[:-1]); bin_width=0.5*(mid[1]-mid[0])
    h,herr = dict(),dict()
    
    # DATA
    for sam,label in zip([selected_beam_off,selected_beam_on],['beam off','beam on']):#{
        h[label],_ = np.histogram( mul*sam[xvar] , bins=bins )
        herr[label] = np.sqrt(h[label])
    #}
    # scale beam-off to beam-on exposure time
    h['beam off'] = h['beam off']*OffBeam_scaling
    herr['beam off'] = h['beam off']*OffBeam_scaling

    # MC
    for sam,label in zip([selected_overlay_concat,selected_CC1p],['overlay','CC1p']):#{
        h[label],_ = np.histogram( mul*sam[xvar] , bins=bins )
        h[label] = h[label]*Nevents['f(POT)']
        herr[label] = np.sqrt(h[label])*Nevents['f(POT)']
    #}
    h['ovrelay + beam off'] = h['overlay'] + h['beam off']
    herr['overlay + beam off'] = np.sqrt(np.square(herr['overlay']) + np.square(herr['beam off']))
    h['background'] = h['ovrelay + beam off'] - h['CC1p']
    herr['background'] = np.sqrt(np.square(herr['overlay + beam off']) + np.square(herr['CC1p']))
    h['beam on bkg subtracted'] = h['beam on'] - h['background']
    herr['beam on bkg subtracted'] = np.sqrt(np.square(herr['beam on']) + np.square(herr['background']))
    if do_draw:#{
        fig = plt.figure(figsize=(20,6))
        ax=fig.add_subplot(1,2,1)
        ax.bar( mid , h['overlay']+h['beam off'] , width=2*bin_width , color=Colors['beam off'], label='beam-off',alpha=0.7)
        ax.bar( mid , h['overlay'] , width=2*bin_width, color=Colors['overlay'], label='overlay',alpha=0.7)
        ax.bar( mid , h['CC1p'], width=2*bin_width, color=Colors['CC1p'], label=r'CC1p signal')
        ax.errorbar( x=mid , xerr=bin_width, y=h['beam on'], yerr=herr['beam on'] , fmt='o', markersize=10
                    , color=Colors['beam on'], capsize=1, capthick=3, label='beam-on')
        set_axes(ax,xlabel,do_add_legend=True)
        ylim = ax.get_ylim()
        ax=fig.add_subplot(1,2,2)
        ax.bar( mid , h['CC1p'], width=2*bin_width, color=Colors['CC1p'])
        ax.errorbar( x=mid , xerr=bin_width, y=h['beam on bkg subtracted'], yerr=herr['beam on bkg subtracted'] , fmt='o', markersize=10
                                , color=Colors['beam on'], capsize=1, capthick=3, label='beam-on')
        set_axes(ax,xlabel,do_add_legend=True,ylim=ylim)
    #}
    np.savetxt(Paths['background maps'] + extra_name + "beam_on_bkg_sbtrctd_%s_%d_bins.csv"%(xvar,len(bins)-1), h['beam on bkg subtracted'], delimiter=",")
    np.savetxt(Paths['background maps'] + extra_name + "beam_on_bkg_sbtrctd_err_%s_%d_bins.csv"%(xvar,len(bins)-1), herr['beam on bkg subtracted'], delimiter=",")
    if debug: print 'saved',extra_name + "beam_on_bkg_sbtrctd_%s_%d_bins.csv"%(xvar,len(bins)-1)
    
    np.savetxt(Paths['background maps'] + extra_name + "mc_cc1p_%s_%d_bins.csv"%(xvar,len(bins)-1), h['CC1p'], delimiter=",")
    np.savetxt(Paths['background maps'] + extra_name + "mc_cc1p_err_%s_%d_bins.csv"%(xvar,len(bins)-1), herr['CC1p'], delimiter=",")
    if debug: print 'saved',extra_name + "mc_cc1p_%s_%d_bins.csv"%(xvar,len(bins)-1)
    if debug: print 'done background subtraction from data'
# ----------------------------------------------------------


# ----------------------------------------------------------
# Aug-27, 2018
def build_migration_map(sam=None,bins=None,xvar='',yvar='',mul=1,debug=0,do_draw=True
                        ,extra_name=''):#{
    xbins = ybins = bins*(1./mul)
    nbins = len(bins)-1
    migration_map = np.zeros((nbins,nbins))
    for ix in range(nbins):#{
        xmin,xmax = xbins[ix],xbins[ix+1]
        sum_in_column = 0
        for iy in range(nbins):#{
            ymin,ymax = ybins[iy],ybins[iy+1]
            sam_bin = sam[(sam[xvar]>xmin)&(sam[xvar]<xmax)&(sam[yvar]>ymin)&(sam[yvar]<ymax)]
            migration_map[ix][iy] = float(len(sam_bin))
            sum_in_column += migration_map[ix][iy]
        #}
        # now normalize the entire column to 1
        for iy in range(nbins):#{
            migration_map[ix][iy] /= sum_in_column
        #}
    #}
    migration_map=migration_map.T
    if debug: print 'done computing migration matrix'
    if do_draw:#{
        fig,ax = plt.subplots(figsize=(9.708,6))
        sns.heatmap(migration_map,annot=True,fmt=".2f",cbar=False,cmap='jet')
        set_axes(ax,'True bin $j$','Reconstructed bin $i$')
        ax.invert_yaxis()
        plt.tight_layout()
    #}
    filename = Paths['migration maps'] + extra_name + "%s_vs_%s_%d_bins.csv"%(xvar,yvar,nbins)
    np.savetxt(filename, migration_map, delimiter=",")
    if debug: print 'done. saved migration map into',filename
    return migration_map
#}
# ----------------------------------------------------------


# ----------------------------------------------------------
# Aug-27, 2018
def smear_MC(in_sample=None,migration_maps=None,debug=0,name='smeared_mc_tmp'
             ,extra_name=''):#{
    sam = in_sample
    debug=0
    smeared_Pmu_array,smeared_Pp_array,smeared_Pmu_cos_theta_array,smeared_Pp_cos_theta_array = [],[],[],[]
    for i,row in sam.iterrows():#{
        
        # smearing p(muon)
        true_bin_j = find_bin(row['truth_Pmu'], Bins['Pmu'])
        smeared_Pmu_array.append(choice(a=0.5*(Bins['Pmu'][1:]+Bins['Pmu'][:-1]) , p=migration_maps['Pmu'][:,true_bin_j]))
        # smearing cos(theta(muon))
        true_bin_j = find_bin(row['truth_Pmu_cos_theta'], Bins['cos(theta(mu))'])
        smeared_Pmu_cos_theta_array.append(choice(a=0.5*(Bins['cos(theta(mu))'][1:]+Bins['cos(theta(mu))'][:-1]) , p=migration_maps['Pmu_cos_theta'][:,true_bin_j]))
        # smearing p(proton)
        true_bin_j = find_bin(row['truth_Pp'], Bins['Pp'])
        smeared_Pp_array.append(choice(a=0.5*(Bins['Pp'][1:]+Bins['Pp'][:-1]) , p=migration_maps['Pp'][:,true_bin_j]))
        # smearing cos(theta(p))
        true_bin_j = find_bin(row['truth_Pp_cos_theta'], Bins['cos(theta(p))'])
        smeared_Pp_cos_theta_array.append(choice(a=0.5*(Bins['cos(theta(p))'][1:]+Bins['cos(theta(p))'][:-1]) , p=migration_maps['Pp_cos_theta'][:,true_bin_j]))
    #}
    print 'done looping.'
    sam['smeared_Pmu']=smeared_Pmu_array
    sam['smeared_Pmu_cos_theta']=smeared_Pmu_cos_theta_array
    sam['smeared_Pp']=smeared_Pp_array
    sam['smeared_Pp_cos_theta']=smeared_Pp_cos_theta_array
    print 'done smearing',name
    prefix = Paths['selected events'] + versions['Overlay'] + '_' + versions['overlay date'] + '_' + extra_name
    sam.to_csv(prefix + name + '.csv')
    print 'saved ',name,' to\n',(prefix + name + '.csv')
    return sam
# ----------------------------------------------------------

# ----------------------------------------------------------
# Aug-27, 2018 (last edit Aug-31,2018)
def get_eff(Ngen=1,Nsel=1,debug=0):#{
    '''
        return: eff, eff_err
        '''
    eff = np.min( [1 , float(Nsel)/Ngen if Ngen>0 else 0] )
    eff_err = eff * np.min( [1 , np.sqrt((1./Nsel if Nsel>0 else 0) + (1./Ngen if Ngen>0 else 0))] )
    if debug: print 'eff = %.4f +/ %.4f'%(eff,eff_err)
    return eff,eff_err
#}
def get_eff_samples(generated=None,selected=None,debug=0):#{
    '''
        return: eff, eff_err
        '''
    Ngen = float(len(generated))
    Nsel = float(len(selected))
    return get_eff(Ngen=Ngen,Nsel=Nsel,debug=debug)
#}
#---------------------------------------------------------------------------------------------

# ----------------------------------------------------------
# Aug-27, 2018 (last edit Sep-3, 2018)
def compute_effiency(genie_CC1p=None
                     ,selected_CC1p=None
                     ,bins=Bins['Pmu']
                     ,xvar='smeared_Pmu'
                     ,xlabel=(Vlabels['Pmu']+Units['Pmu'])
                     ,ylabel=r'$\bar{\epsilon}$'
                     ,do_draw=True,debug=0
                     ,mul=1
                     ,extra_name=''
                     ):#{
    mid = 0.5*(bins[1:]+bins[:-1]); bin_width=0.5*(mid[1]-mid[0])
    h = dict()
    h['mid'] = mid
    h['bin width'] = bin_width
    h['generated'],_ = np.histogram(mul*genie_CC1p[xvar],bins=bins)
    h['selected'],_ = np.histogram(mul*selected_CC1p[xvar],bins=bins)
    
    h['eff'],h['eff err'] = [],[]
    for i in range(len(bins)-1):#{
        eff,eff_err = get_eff(Ngen=h['generated'][i],Nsel=h['selected'][i])
        h['eff'].append(eff)
        h['eff err'].append(eff_err)
    #}
    if debug: print 'done computing efficiency.'

    if do_draw:#{
        fig=plt.figure(figsize=(20,6))
        ax=fig.add_subplot(1,2,1)
        for label,color in zip(['generated','selected'],['forestgreen','royalblue']):
            plt.errorbar(x=mid,xerr=bin_width,y=h[label],yerr=np.sqrt(h[label])
                         ,color=color,capsize=10,fmt='.',markersize=0,label=label)
        set_axes(ax,xlabel,'counts',do_add_grid=True,do_add_legend=True)
        ax=fig.add_subplot(1,2,2)
        plt.errorbar(x=mid,xerr=bin_width,y=h['eff'],yerr=h['eff err']
                 ,color='black',capsize=10,fmt='.',markersize=0)
        set_axes(ax,xlabel,ylabel,do_add_grid=True,ylim=(0,1.05*np.max(h['eff']+h['eff err'])))
        if debug: print 'done drawing.'
    #}

    # save to csv
    np.savetxt(Paths['efficiency maps'] + extra_name + "eff_%s_%d_bins.csv"%(xvar,len(bins)-1), h['eff'], delimiter=",")
    np.savetxt(Paths['efficiency maps'] + extra_name + "eff_err_%s_%d_bins.csv"%(xvar,len(bins)-1), h['eff err'], delimiter=",")
    if debug: print 'saved efficiency into',Paths['efficiency maps'] + extra_name + "eff_%s_%d_bins.csv"%(xvar,len(bins)-1),'\n and uncertainty'
    return h
#}
# ----------------------------------------------------------



# ----------------------------------------------------------
# Aug-27, 2018 (last edit Aug-31, 2018)
def sample_in_limits(sam=None
                     ,varPmu='reco_Pmu_mcs',varPmu_cos_theta='reco_Pmu_cos_theta',varPmu_phi='reco_Pmu_mcs_phi'
                     ,varPp='reco_Pp',varPp_cos_theta='reco_Pp_cos_theta',varPp_phi='reco_Pp_phi'
                     ):#{
    return sam[ (Limits['Pmu'][0] <= sam[varPmu]) & (sam[varPmu]<=Limits['Pmu'][1])
               & (Limits['cos(theta(mu))'][0] <= sam[varPmu_cos_theta]) & (sam[varPmu_cos_theta]<=Limits['cos(theta(mu))'][1])
               & (Limits['Pp'][0] <= sam[varPp])  & (sam[varPp] <= Limits['Pp'][1])
               & (Limits['cos(theta(p))'][0] <= sam[varPp_cos_theta]) & (sam[varPp_cos_theta] <= Limits['cos(theta(p))'][1])
               ]
#}
#---------------------------------------------------------------------------------------------



# ----------------------------------------------------------
# Aug-27, 2018 (last edit Aug-31, 2018)
def load_mc_and_data(extra_name=''
                     ,minPEcut = 150
                     ,maxdYZcut = 200
                     ,delta_theta_12 = 55
                     ,r_max_RdQ_CC1p = 0.43
                     ,delta_Delta_phi = 35
                     ,Pt_max = 0.35
                     ,Chi2Proton_muCandidate_min = 80
                     ,Chi2Proton_pCandidate_max = 30):#{
    # ----------------------------------------------------------
    ## (1) MC
    prefix = Paths['selected events'] + versions['Overlay'] + '_' + versions['overlay date'] + '_' + extra_name
    selected_cosmic_filename = 'selected_cosmic.csv'
    selected_overlay=dict()
    
    cuts_order  = ['no cut','Chi2Proton','Nflashes','MatchedFlash','length'
                   ,'non-collinearity','vertex activity'
                   ,'delta phi','Pt & delta phi']
        
    if os.path.isfile(prefix+selected_cosmic_filename):#{
        print 'found '+selected_cosmic_filename+', loading it...'
        for pair_type in pair_types:#{
            selected_overlay[pair_type]=pd.read_csv(prefix+'selected_'+pair_type+'.csv')
        #}
        selected_overlay_concat = pd.concat([selected_overlay['1mu-1p'],selected_overlay['cosmic'],selected_overlay['other-pairs']])
    #}
    else:#{
        print 'did not find '+selected_cosmic_filename+', so creating it...'
        OverlaySamples = load_samples(date=versions['overlay date'],filename=versions['Overlay']+'_'+versions['overlay date']+'_vertices')
        reducedOverlay,pureffOverlay,pureffNumbers = apply_cuts_to_overlay(OverlaySamples=OverlaySamples, cuts_order=cuts_order
                                                                           ,minPEcut = minPEcut
                                                                           ,maxdYZcut = maxdYZcut
                                                                           ,delta_theta_12 = delta_theta_12
                                                                           ,r_max_RdQ_CC1p = r_max_RdQ_CC1p
                                                                           ,delta_Delta_phi=delta_Delta_phi
                                                                           ,Pt_max=Pt_max
                                                                           ,Chi2Proton_muCandidate_min=Chi2Proton_muCandidate_min
                                                                           ,Chi2Proton_pCandidate_max=Chi2Proton_pCandidate_max)
        print 'applied cuts to overlay'
        for pair_type in pair_types:#{
            selected_overlay[pair_type] = sample_in_limits(sam=reducedOverlay['Pt & delta phi'][pair_type])
            outcsvname = prefix+'selected_'+pair_type+'.csv'
            selected_overlay[pair_type].to_csv(outcsvname)
            print 'saved selected',pair_type,'to',outcsvname
        #}
        # overlay scaling
        summary = pd.read_csv('/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/summary/'
                              +versions['overlay date']+'/'
                              +versions['Overlay']+'_'+versions['overlay date']+'_summary.csv')
        Nevents['OnBeam POT']   = 4.908e+19
        Nevents['overlay']      = np.sum(summary.Nevents)
        Nevents['overlay POT']  = np.sum(summary.POT)
        Nevents['f(POT)']       = Nevents['OnBeam POT']/Nevents['overlay POT']
        print "Nevents['f(POT)']:",Nevents['f(POT)']
        selected_overlay_concat = pd.concat([selected_overlay['1mu-1p'],selected_overlay['cosmic'],selected_overlay['other-pairs']])
        print len(selected_overlay_concat),'events in the overlay'
    #}
    # ----------------------------------------------------------
    ## (2) DATA
    data_prefix = Paths['selected events'] + versions['beam on'] + '_' + versions['data date'] + '_' + extra_name
    if os.path.isfile(data_prefix+'selected_beam_on.csv'):#{
        print 'checked',data_prefix+'selected_on_beam.csv and found the file...'
        selected_beam_on = pd.read_csv(data_prefix+'selected_beam_on.csv')
        selected_beam_off = pd.read_csv(data_prefix+'selected_beam_off.csv')
    #}
    else:#{
        print 'checked',prefix+'selected_on_beam.csv and there was no file there...'
        OnBeam = pd.read_csv(vertices_files_path+'/'+versions['data date']+'/'+versions['beam on']+'_'+versions['data date']+'_vertices.csv')
        print 'loaded beam-on'
        OffBeam = pd.read_csv(vertices_files_path+'/'+versions['data date']+'/'+versions['beam off']+'_'+versions['data date']+'_vertices.csv')
        print 'loaded beam-off'
        reducedOnBeam,reducedOffBeam,numbers = apply_cuts_to_data(OnBeam=OnBeam,OffBeam=OffBeam,cuts_order=cuts_order
                                                                  ,minPEcut = minPEcut
                                                                  ,maxdYZcut = maxdYZcut
                                                                  ,delta_theta_12 = delta_theta_12
                                                                  ,r_max_RdQ_CC1p = r_max_RdQ_CC1p
                                                                  ,delta_Delta_phi=delta_Delta_phi
                                                                  ,Pt_max=Pt_max
                                                                  ,Chi2Proton_muCandidate_min=Chi2Proton_muCandidate_min
                                                                  ,Chi2Proton_pCandidate_max=Chi2Proton_pCandidate_max)
        print 'applied cuts to data'
        selected_beam_on = sample_in_limits(sam=reducedOnBeam['Pt & delta phi'])
        outcsvname = data_prefix+'selected_beam_on.csv'
        selected_beam_on.to_csv(outcsvname)
        print 'saved %d'%len(selected_beam_on),'events in beam on to',outcsvname
        
        selected_beam_off = sample_in_limits(sam=reducedOffBeam['Pt & delta phi'])
        outcsvname = data_prefix+'selected_beam_off.csv'
        selected_beam_off.to_csv(outcsvname)
        print 'saved %d'%len(selected_beam_off),'events in beam off to',outcsvname
    #}
    # ----------------------------------------------------------
    ## (3) GENIE
    if os.path.isfile(prefix + 'selected_genie_CC1p.csv'):#{
        print 'checked',prefix+'selected_genie_CC1p.csv and found the file...'
        genie_CC1p = pd.read_csv(prefix+'selected_genie_CC1p.csv')
    #}
    else:#{
        overlay_genie = pd.read_csv('/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/genie/'
                                    +versions['overlay date']+'/'
                                    +versions['Overlay']+'_'+versions['overlay date']+'_genie.csv')

        print len(overlay_genie),'events in genie from overlay'
        overlay_genie_CC1p = overlay_genie[(overlay_genie.IsCC_1p_200MeVc==True)
                                           & ((overlay_genie.truth_x>3) & (overlay_genie.truth_x<256))
                                           & ((overlay_genie.truth_y>-115) & (overlay_genie.truth_y<115))
                                           & ((overlay_genie.truth_z>5) & (overlay_genie.truth_z<1037))
                                           ]
                                           
        overlay_genie_CC1p_in_limits = sample_in_limits(sam=overlay_genie_CC1p
                                                ,varPmu='truth_Pmu',varPmu_cos_theta='truth_Pmu_cos_theta'
                                                ,varPp='truth_Pp',varPp_cos_theta='truth_Pp_cos_theta'
                                                )
        genie_CC1p = overlay_genie_CC1p_in_limits
        outcsvname = prefix+'selected_genie_CC1p.csv'
        genie_CC1p.to_csv(outcsvname)
        print 'saved %d'%len(genie_CC1p),'CC1p events in genie_CC1p to',outcsvname
    #}
    selected_CC1p = selected_overlay['CC1p']
    print len(selected_CC1p),'selected CC1p events overlay'
    return selected_overlay,selected_overlay_concat,selected_CC1p,genie_CC1p,selected_beam_on,selected_beam_off
#}
# ----------------------------------------------------------

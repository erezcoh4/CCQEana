import sys, pandas as pd, numpy as np, matplotlib.pyplot as plt
sys.path.insert(0, '/Users/erezcohen/larlite/UserDev/mySoftware/MySoftwarePackage/mac')
from plot_tools import *

import sys; sys.path.insert(0, '../..');
from ccqe_notebook_tools import *
from mupClassification_notebooks import *
from onbeam_offbeam_notebooks import *
from mpl_toolkits.mplot3d import Axes3D
import cPickle as pickle

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
Ntargets = 1.1782e30
Ntargets_err = 0
# kinematical cuts
delta_theta_12=55  # deg.
delta_Delta_phi=35 # deg.
Pt_max=0.35        # GeV/c


Vars=dict({
          'Pmu':'reco_Pmu_mcs'
          ,'cos(theta(mu))':'reco_Pmu_cos_theta'
          ,'phi(mu)':'reco_Pmu_mcs_phi'
          ,'Pp':'reco_Pp'
          ,'cos(theta(p))':'reco_Pp_cos_theta'
          ,'phi(p)':'reco_Pp_phi'
          })

Limits=dict({
            'Pmu':(0.1,1.5)
            ,'cos(theta(mu))':(-0.65,0.95)
            ,'phi(mu)':(-180,180)
            ,'Pp':(0.3,1.0)
            ,'cos(theta(p))':(0.15,1.0)
            ,'phi(p)':(-180,180)
            })
NBins=7
Bins = dict()
for key in Limits.keys(): Bins[key] = np.linspace(Limits[key][0],Limits[key][1],NBins+1)
# Bins['cos(theta(mu))'] = np.array([-0.65,-0.4,-0.15,0.1,0.35,0.6,0.9,0.95])
# Oct-2, 2018
Bins['cos(theta(mu))'] = np.array([-0.65, -0.4083,-0.167,0.075,0.317,0.5583, 0.8, 0.95])
Bins['cos(theta(p))'] = np.array([ 0.15,0.2583,0.366,0.475 ,0.583,0.692,0.8,0.95])
bins1,bins2,bins3 = Bins['Pmu'], Bins['cos(theta(mu))'] , Bins['phi(mu)']
N1,N2,N3 = len(bins1)-1,len(bins2)-1,len(bins3)-1
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
Colors = dict({'overlay':'forestgreen','CC1p':'blue','beam off':'black','beam on':'tomato'})

Xsec_path = '/Users/erezcohen/Desktop/uBoone/CCQEanalysis/Xsec/'
Paths = dict({'selected events':Xsec_path+'selected_events/'
             ,'migration maps':Xsec_path+'migration_maps/'
             ,'background maps':Xsec_path+'background_maps/'
             ,'efficiency maps':Xsec_path+'efficeincy_maps/'
             ,'1d Xsec':Xsec_path+'1d_Xsec/'
             ,'systematics':Xsec_path+'systematics/'})

Xsec_ctu_titles = [r'excluding the last $\cos\theta_\mu$ bin',r'with the last $\cos\theta_\mu$ bin']
Xsec_fnames     = [r'without_last_ctu_bin',r'with_last_ctu_bin']
remove_ctu_bools= [True,False]


# ----------------------------------------------------------
# Oct-09, 2018
def extract_Xsec_full_chain(extra_name='',debug=0
                            ,minPEcut = 150,maxdYZcut = 200
                            ,delta_theta_12 = 55,r_max_RdQ_CC1p = 0.43
                            ,delta_Delta_phi = 35,Pt_max =0.35
                            ,Chi2Proton_muCandidate_min =80,Chi2Proton_pCandidate_max = 30,force_recalculated_weights=False):
    samples = load_mc_and_data(extra_name=extra_name,debug=debug
                               ,minPEcut = minPEcut,maxdYZcut = maxdYZcut
                               ,delta_theta_12 = delta_theta_12,r_max_RdQ_CC1p = r_max_RdQ_CC1p
                               ,delta_Delta_phi = delta_Delta_phi,Pt_max = Pt_max
                               ,Chi2Proton_muCandidate_min = Chi2Proton_muCandidate_min,Chi2Proton_pCandidate_max = Chi2Proton_pCandidate_max)
    selected_overlay,selected_overlay_concat,selected_CC1p,genie_CC1p,selected_beam_on,selected_beam_off = samples
    print 'done loading samples...'
    if (("Pmu weight" not in selected_beam_on.columns) or ("Pp weight" not in selected_beam_on.columns) or force_recalculated_weights):#{
        print 'no Pmu weights, computing them'
        compute_eff_weights(beam_on=selected_beam_on,beam_off=selected_beam_off,debug=debug,
                                                           generated_CC1p=genie_CC1p,selected_CC1p=selected_CC1p,overlay=selected_overlay_concat,
                                                           delta_theta_12=delta_theta_12,
                                                           delta_Delta_phi=delta_Delta_phi,
                                                           Pt_max=Pt_max)
        compute_eff_weights(beam_on=selected_beam_on,beam_off=selected_beam_off,debug=debug,
                                                               generated_CC1p=genie_CC1p,selected_CC1p=selected_CC1p,overlay=selected_overlay_concat,
                                                               delta_theta_12=delta_theta_12,
                                                               delta_Delta_phi=delta_Delta_phi,
                                                               Pt_max=Pt_max,
                                                               ob_1='Pp',ob_2='cos(theta(p))',ob_3='phi(p)',
                                                               reco_1='reco_Pp',reco_2='reco_Pp_cos_theta',reco_3='reco_Pp_phi',
                                                               true_1='truth_Pp',true_2='truth_Pp_cos_theta',true_3='truth_Pp_phi')
        print 'done assiging Pmu weights and Pp weights and saving the files.'
    #}
    else: print 'Pmu weights and Pp weights already exist.'
    # iterative process for correction around \phi~0
    if ('W(corr. phi~0)' not in selected_beam_on.columns):#{
        for sam in [selected_beam_on,selected_beam_off,selected_CC1p,selected_overlay_concat]: sam['W(corr. phi~0)'] = 1
        Xsec_dict = get_phi_Xsecs( do_corr_phi_0=False , debug=debug ,
                                  selected_beam_on=selected_beam_on,selected_beam_off=selected_beam_off,
                                  selected_overlay_concat=selected_overlay_concat,selected_CC1p=selected_CC1p)
        correction, correction_arrays = dict(),dict()
        correction_arrays['mu'] = []; correction_arrays['p'] = []; correction['mu'] = correction['p'] = 1
        current_correction = 100
        while np.abs(current_correction-1)>0.0001:#{
            for particle,var in zip(['mu','p'],['reco_Pmu_mcs_phi','reco_Pp_phi']):#{
                Xsec,Xsec_err = Xsec_dict['phi('+particle+')'],Xsec_dict['phi('+particle+') err']
                bins = Bins['phi('+particle+')']
                N = len(bins)-1; n = N/2
                Xsec_phi_0 = Xsec[n]
                Xsec_phi_not_0 = np.concatenate([Xsec[:n-1],Xsec[n+1:]])
                Xsec_err_phi_not_0 = np.concatenate([Xsec_err[:n-1],Xsec_err[n+1:]])
                mean_phi_not_0 = np.average(Xsec_phi_not_0 , weights=1./np.square(Xsec_err_phi_not_0))
                current_correction = mean_phi_not_0/Xsec_phi_0
                correction[particle] = correction[particle]*current_correction
                correction_arrays[particle].append( correction[particle] )
                for sam in [selected_beam_on,selected_beam_off,selected_CC1p,selected_overlay_concat]:#{
                    indices_phi_0 = sam.index[(bins[n] <= 180./np.pi*sam[var]) & (180./np.pi*sam[var] < bins[n+1])].tolist()
                    sam.loc[indices_phi_0,'W(corr. phi~0)'] = sam.loc[indices_phi_0,'W(corr. phi~0)'] * current_correction
                    Xsec_dict = get_phi_Xsecs( do_corr_phi_0=True , debug=debug ,
                                              selected_beam_on=selected_beam_on,selected_beam_off=selected_beam_off,
                                              selected_overlay_concat=selected_overlay_concat,selected_CC1p=selected_CC1p)
                #}
            #}
        #}
        for particle in ['mu','p']: print 'correction for phi('+particle+')~0:',correction_arrays[particle][-1]
        save_selected_samples(selected_overlay_concat , selected_CC1p , selected_beam_on , selected_beam_off, extra_name=extra_name)
        print 'done performing iterative correction for phi~0 and saved the samples...'    
    else: print 'already performed correction for phi~0 and saved the samples...'
    extract_Xsecs(debug=debug,extra_name=extra_name,particle='mu',
              selected_beam_on=selected_beam_on,selected_beam_off=selected_beam_off,                  
              selected_overlay_concat=selected_overlay_concat,selected_CC1p=selected_CC1p)
    extract_Xsecs(debug=debug,extra_name=extra_name,particle='p',
                  selected_beam_on=selected_beam_on,selected_beam_off=selected_beam_off,                  
                  selected_overlay_concat=selected_overlay_concat,selected_CC1p=selected_CC1p)
# ----------------------------------------------------------


# ----------------------------------------------------------
# Oct-08, 2018
def get_phi_Xsecs(do_corr_phi_0=False, debug=0,
                  selected_beam_on=None,selected_beam_off=None,selected_overlay_concat=None,selected_CC1p=None):#{
    Xsec_dict = dict()
    for particle in ['mu','p']:#{
        observable='phi('+particle+')'
        var,bins,mid,bin_width,vlabel,xlabel,units = get_labels(observable=observable)
        h = get_Xsec_1d(selected_beam_on,selected_beam_off,selected_overlay_concat,selected_CC1p
                        ,var=var,bins=bins,bin_width=bin_width
                        ,wname='P'+particle+' weight',mul=180./np.pi
                        ,do_corr_phi_0=do_corr_phi_0)
        Xsec_dict[observable] = h['Xsec']
        Xsec_dict[observable+' err'] = h['Xsec err']
        Xsec_dict['mc '+observable] = h['mc Xsec']
        Xsec_dict['mc '+observable+' err'] = h['mc Xsec err']
    #}
    return Xsec_dict
#}
# ----------------------------------------------------------


# ----------------------------------------------------------
# Oct-08, 2018
def extract_Xsecs(do_corr_phi_0=True, debug=0, particle='mu',
                  selected_beam_on=None,selected_beam_off=None,
                  selected_overlay_concat=None,selected_CC1p=None,
                  extra_name=""):#{
    Xsec_dicts = dict()
    for iXsec,(Xsec_title,remove_ctu_bin) in enumerate(zip(Xsec_ctu_titles,remove_ctu_bools)):#{
        Xsec_dict = get_Xsecs(do_corr_phi_0=do_corr_phi_0, debug=debug, particle=particle,
                              remove_last_cos_theta_mu_bin=remove_ctu_bin,
                              do_P=True, do_cos_theta=True, do_phi=True,
                              selected_beam_on=selected_beam_on,selected_beam_off=selected_beam_off,
                              selected_overlay_concat=selected_overlay_concat,selected_CC1p=selected_CC1p)
        Xsec_dicts[Xsec_title] = Xsec_dict
    #}
    if debug: pp.pprint(Xsec_dicts)
    outfilename = Paths['1d Xsec'] + "P"+particle + "Xsecs_1D"+extra_name+".txt"
    with open(outfilename, 'w') as outfile:#{
        outfile.write(pickle.dumps(Xsec_dicts))
    #}
    print 'saved cross-sections into',outfilename
#}
# ----------------------------------------------------------

# ----------------------------------------------------------
# Oct-08, 2018
def get_Xsec_variable(debug=0,
                      var='reco_Pt',mul=1,bins=linspace(0,1,5),
                      wname='Pmu weight',
                      do_corr_phi_0=True, remove_ctu_bin=True,
                      selected_beam_on=None,selected_beam_off=None,
                      selected_overlay_concat=None,selected_CC1p=None):#{
    
    Xsec_dict = dict()
    beam_on, beam_off, overlay, CC1p = selected_beam_on,selected_beam_off,selected_overlay_concat,selected_CC1p
    if remove_ctu_bin:#{
        beam_on = beam_on[beam_on['reco_Pmu_cos_theta']<Bins['cos(theta(mu))'][-2]]
        beam_off = beam_off[beam_off['reco_Pmu_cos_theta']<Bins['cos(theta(mu))'][-2]]
        overlay = overlay[overlay['reco_Pmu_cos_theta']<Bins['cos(theta(mu))'][-2]]
        CC1p = CC1p[CC1p['reco_Pmu_cos_theta']<Bins['cos(theta(mu))'][-2]]
    #}
    h = get_Xsec_1d(beam_on,beam_off,overlay,CC1p
                    ,var=var,bins=bins,bin_width=bins[1]-bins[0]
                    ,wname=wname,mul=mul
                    ,do_corr_phi_0=do_corr_phi_0,debug=debug)
    if debug>1:  pp.pprint(h)
    Xsec_dict[var] = h['Xsec']
    Xsec_dict[var+' err'] = h['Xsec err']
    Xsec_dict['mc '+var] = h['mc Xsec']
    Xsec_dict['mc '+var+' err'] = h['mc Xsec err']
    #}
    return Xsec_dict
# ----------------------------------------------------------


# ----------------------------------------------------------
# Oct-08, 2018
def draw_Xsec_variable(debug=0,
                       var='reco_Pt',mul=1,bins=linspace(0,1,5),vlabel='p_T',units=None,legend_loc='best',
                       wname='Pmu weight',
                       do_corr_phi_0=True,
                       selected_beam_on=None,selected_beam_off=None,
                       selected_overlay_concat=None,selected_CC1p=None,
                       filename=None,extra_name='',
                       residuals_ylim=[-1,1],residuals_yticks=[-0.5,0,0.5],residuals_ytitle=1.05,residuals_xtitle='center',figures_path='~/Desktop/'):#{
    Xsec_dicts = dict()
    for Xsec_ctu_title,remove_ctu_bin in zip(Xsec_ctu_titles,remove_ctu_bools):
        Xsec_dicts[Xsec_ctu_title] = get_Xsec_variable(debug=debug,
                                                       var=var,mul=mul,bins=bins,
                                                       wname=wname,
                                                       do_corr_phi_0=do_corr_phi_0,
                                                       remove_ctu_bin=remove_ctu_bin,
                                                       selected_beam_on=selected_beam_on,selected_beam_off=selected_beam_off,
                                                       selected_overlay_concat=selected_overlay_concat,selected_CC1p=selected_CC1p)

    fig=plt.figure(figsize=(20,8))
    mid = 0.5*(bins[1:]+bins[:-1]); bin_width=mid[1]-mid[0]
    for iXsec,(Xsec_ctu_title,iax) in enumerate(zip(Xsec_ctu_titles,[(1,3),(2,4)])):
        Xsec_dict = Xsec_dicts[Xsec_ctu_title]
        h = dict()
        ax = fig.add_subplot(3,2,iax)
        h['Xsec'],h['Xsec err'] = Xsec_dict[var],Xsec_dict[var+' err']
        plt.errorbar(x=mid,xerr=0.5*bin_width,y=h['Xsec'],yerr=h['Xsec err'],color=Colors['beam on'],fmt='o',label='data')
        h['mc Xsec'],h['mc Xsec err'] = Xsec_dict['mc '+var], Xsec_dict['mc '+var+' err']
        ax.bar( x=mid , height=2*h['mc Xsec err'], bottom=h['mc Xsec']-h['mc Xsec err'], width=bin_width, color=Colors['CC1p'],label='overlay')
        set_axes(ax,x_label='',y_label=get_Xsec_label(vlabel,units)
                 ,do_add_grid=True,remove_ticks_x=True,do_add_legend=True if iXsec==1 else False, legend_loc=legend_loc
                 ,ylim=(0,1.1*np.max(ax.get_ylim())),title=Xsec_ctu_title)
        # residuals plot
        den, den_err = h['Xsec']-h['mc Xsec'],np.sqrt(np.square(h['Xsec err'])+np.square(h['mc Xsec err']))
        num, num_err = h['Xsec'],h['mc Xsec err']
        ratio = den/num
        ax = fig.add_subplot(3,2,5+iXsec)
        ratio_err = ratio*np.sqrt(np.square(den_err/den) + np.square(num_err/num))
        plt.errorbar(x=mid,xerr=0.5*bin_width,y=ratio,yerr=ratio_err,fmt='o',markersize=0,color='black')
        chi2,ndf = chi2_two_data_curves(h1=h['Xsec'],h1err=h['Xsec err'],h2=h['mc Xsec'],h2err=h['mc Xsec err'],bins=bins,debug=debug)
        if units is None: x_label = '$'+vlabel+'$'
        elif '^' in units or '_' in units: x_label = '$'+vlabel+'$' + r' $%s$'%units
        else: x_label = '$'+vlabel+'$' + r' [%s]'%units
        set_axes(ax,x_label=x_label
                 ,y_label=r'(data-MC)/data',do_add_grid=True
                 ,ylim=residuals_ylim,yticks=residuals_yticks)
        plt.plot([np.min(bins),np.max(bins)],[0,0],'--',color='royalblue')
        plt.title(r'$\chi^2/ndf=%.2f/%d$'%(chi2,ndf), fontsize=20,y=residuals_ytitle,loc=residuals_xtitle)
    plt.tight_layout(h_pad=0.0)
    plt.subplots_adjust(hspace=0.05)
    if filename is not None: outfilename = figures_path + filename + '.pdf'
    else: outfilename = figures_path + var + '_Xsec'+extra_name+'.pdf'
    save_figure(outfilename)
    return ax
#}
# ----------------------------------------------------------

# ----------------------------------------------------------
# Oct-08, 2018
def get_Xsec_label(vlabel=None,units=None,power_factor=38):
    factor_str = r'10^{- %d}'%power_factor

    if vlabel is None:
        return r'$\times '+factor_str+' \mathrm{cm}^{2}$'

    prefix = r'$\frac{d\sigma}{d'+vlabel+'}$ '
    
    if units is not None:
        suffix = r'$\left['+factor_str+r' \frac{\mathrm{cm}^{2}}{%s} \right]$'%units
    else:
        suffix = r'$\left['+factor_str+r' \mathrm{cm}^{2}\right]$'

    return  prefix + suffix
#}
# ----------------------------------------------------------

# ----------------------------------------------------------
# Oct-04, 2018
def save_selected_samples(selected_overlay_concat , selected_CC1p , selected_beam_on , selected_beam_off,extra_name=''):
    overlay_prefix = Paths['selected events'] + versions['Overlay'] + '_' + versions['overlay date'] + '_' + extra_name
    data_prefix = Paths['selected events'] + versions['beam on'] + '_' + versions['data date'] + '_' + extra_name
    for sam,name,prefix in zip([selected_overlay_concat,selected_CC1p,selected_beam_on,selected_beam_off]
                               ,['overlay','CC1p','beam_on','beam_off']
                               ,[overlay_prefix,overlay_prefix,data_prefix,data_prefix]):#{
        outcsvname = prefix+'selected_'+name+extra_name+'.csv'
        sam.to_csv(outcsvname)
        print 'saved ',len(sam),'selected '+name+' events to',outcsvname
#}
# ----------------------------------------------------------



# ----------------------------------------------------------
# Oct-03, 2018
def get_Xsecs(do_corr_phi_0=False, debug=0, particle='mu', do_P=True, do_cos_theta=True, do_phi=True, do_print_Xsec=False,
              remove_last_cos_theta_mu_bin=False,
              selected_beam_on=None,selected_beam_off=None,selected_overlay_concat=None,selected_CC1p=None,
              extra_wname=""):#{
    Xsec_dict = dict()
    for i,(observable,true,do_var) in enumerate(zip(['P'+particle,'cos(theta('+particle+'))','phi('+particle+')']
                                                             ,['truth_P'+particle,'truth_P'+particle+'_cos_theta','truth_P'+particle+'_phi']
                                                             ,[do_P,do_cos_theta,do_phi])):#{
        if do_var==False: continue
        var,bins,mid,bin_width,vlabel,xlabel,units = get_labels(observable=observable)
        mul = 180./np.pi if 'phi' in observable else 1
        beam_on , beam_off , overlay , CC1p = selected_beam_on,selected_beam_off,selected_overlay_concat,selected_CC1p
        if remove_last_cos_theta_mu_bin:#{
            beam_on = beam_on[beam_on['reco_Pmu_cos_theta']<Bins['cos(theta(mu))'][-2]]
            beam_off = beam_off[beam_off['reco_Pmu_cos_theta']<Bins['cos(theta(mu))'][-2]]
            overlay = overlay[overlay['reco_Pmu_cos_theta']<Bins['cos(theta(mu))'][-2]]
            CC1p = CC1p[CC1p['reco_Pmu_cos_theta']<Bins['cos(theta(mu))'][-2]]
        #}
        h = get_Xsec_1d(beam_on,beam_off,overlay,CC1p
                        ,var=var,bins=bins,bin_width=bin_width,wname='P'+particle+' weight'+extra_wname,mul=mul
                        ,do_corr_phi_0=do_corr_phi_0)
        if i==0:#{
            Xsec_dict['integrated Xsec'] = np.sum(h['Xsec']*bin_width)
            Xsec_dict['integrated Xsec err'] = np.sqrt(np.sum(np.square(h['Xsec err'])*bin_width))
            Xsec_dict['mc Xsec'] = np.sum(h['mc Xsec']*bin_width)
            Xsec_dict['mc Xsec err'] = np.sqrt(np.sum(np.square(h['mc Xsec err'])*bin_width))
            if do_print_Xsec: print ('integrated Xsec: %.2f+/-%.2f'%(Xsec_dict['integrated Xsec'],Xsec_dict['integrated Xsec err']),
                                            'mc Xsec: %.2f+/-%.2f'%(Xsec_dict['mc Xsec'],Xsec_dict['mc Xsec err']))
        #}
        if debug>1:  pp.pprint(h)
        Xsec_dict[observable] = h['Xsec']
        Xsec_dict[observable+' err'] = h['Xsec err']
        Xsec_dict['mc '+observable] = h['mc Xsec']
        Xsec_dict['mc '+observable+' err'] = h['mc Xsec err']
    #}
    return Xsec_dict
#}
# ----------------------------------------------------------


# ----------------------------------------------------------
# Oct-03, 2018 (last edit Oct-04)
# computation of 1D cross-section based on weighted distirubtions
def get_Xsec_1d(beam_on=None,beam_off=None,overlay=None,CC1p=None,
                var='reco_Pmu_mcs',bins=Bins['Pmu'],bin_width=None,wname='Pmu weight',
                mul=1,
                do_corr_phi_0=False,
                debug=0):#{
    h=dict()
    for sam,slabel in zip([beam_on,beam_off,overlay,CC1p]
                          ,['beam on','beam off','overlay','CC1p']):#{
        h[slabel],h[slabel+' err']=np.zeros(len(bins)-1),np.zeros(len(bins)-1)
        for i in range(len(bins)-1):#{
            sam_in_bin = sam[(bins[i]<=mul*sam[var])& (mul*sam[var]<bins[i+1])]
            if debug>1: print "len(sam in bins[%d"%i+"]):",len(sam_in_bin)
            h[slabel][i] = np.sum(sam_in_bin[wname])
            if do_corr_phi_0:#{
                h[slabel][i] = np.sum(sam_in_bin[wname]*sam_in_bin['W(corr. phi~0)'])
            #}
            h[slabel+' err'][i] = np.sqrt(np.sum(np.square(sam_in_bin[wname+' err'])
                                                 +np.square(sam_in_bin[wname])))
        #}
    #}
    h['B'] = h['overlay'] - h['CC1p']
    h['B err'] = np.sqrt(np.square(h['overlay err']) + np.square(h['CC1p err']))
    
    h['B scaled'] = h['B']*Nevents['f(POT)']
    h['B scaled err'] = h['B err']*Nevents['f(POT)']
    
    #     print "h['beam off err']:",h['beam off err']
    h['beam off scaled'] = h['beam off']*OffBeam_scaling
    h['beam off scaled err'] = h['beam off err']*OffBeam_scaling
    
    h['N(on)-N(off)-B'] = h['beam on'] - h['beam off scaled'] - h['B scaled']
    h['N(on)-N(off)-B err'] = np.sqrt(np.square(h['beam on err'])
                                      + np.square(h['beam off scaled err'])
                                      + np.square(h['B scaled err']))

    h['Xsec'] = h['N(on)-N(off)-B']/bin_width
    h['Xsec err'] = h['N(on)-N(off)-B err']/bin_width
    
    h['Xsec beam on'] = h['beam on']/bin_width
    h['Xsec beam on err'] = h['beam on err']/bin_width
        
    # foc CC1p (mc-Xsec) we want no correction applied
    for i in range(len(bins)-1):#{
        CC1p_in_bin = CC1p[(bins[i]<=mul*CC1p[var])& (mul*CC1p[var]<bins[i+1])]
        h['CC1p'][i] = np.sum(CC1p_in_bin[wname])
        h['mc Xsec'] = h['CC1p']*Nevents['f(POT)']/bin_width
        h['mc Xsec err'] = h['CC1p err']*Nevents['f(POT)']/bin_width
    #}
    return h
#}
# ----------------------------------------------------------








# ----------------------------------------------------------
# Oct-03, 2018
# efficiency weights for cross-section
def compute_eff_weights(beam_on=None,beam_off=None
                        ,generated_CC1p=None,selected_CC1p=None,overlay=None
                        ,ob_1='Pmu',ob_2='cos(theta(mu))',ob_3='phi(mu)'
                        ,reco_1='reco_Pmu_mcs',reco_2='reco_Pmu_cos_theta',reco_3='reco_Pmu_mcs_phi'
                        ,true_1='truth_Pmu',true_2='truth_Pmu_cos_theta',true_3='truth_Pmu_phi'
                        # for fixing muon bins
                        ,do_in_kin_cuts=True
                        ,delta_theta_12=55  # deg.
                        ,delta_Delta_phi=35 # deg.
                        ,Pt_max=0.35        # GeV/c
                        ,debug=0
                        ,eff_cutoff=0.01
                        # for different binning...
                        ,do_different_binning = False
                        ,NBins=7
                        ,bins_cos_theta_mu = None
                        ,bins_cos_theta_p = None
                        # different options for systematical studies
                        ,option=""
                        ,power_factor=38
                        ):
    # for different binning...
    if do_different_binning:#{
        for key in Limits.keys(): Bins[key] = np.linspace(Limits[key][0],Limits[key][1],NBins+1)
        if bins_cos_theta_mu is not None: Bins['cos(theta(mu))'] = bins_cos_theta_mu;
        if bins_cos_theta_p is not None: Bins['cos(theta(p))'] = bins_cos_theta_p
    #}
    global bins1,bins2,bins3,N1,N2,N3
    bins1,bins2,bins3 = Bins[ob_1], Bins[ob_2] , Bins[ob_3]
    N1,N2,N3 = len(bins1)-1,len(bins2)-1,len(bins3)-1
    N,List = dict(),dict()
    for i_P in range(N1):#{
        Pmin,Pmax = bins1[i_P],bins1[i_P+1]
        for i_cos_theta in range(N2):#{
            cos_theta_min,cos_theta_max = bins2[i_cos_theta],bins2[i_cos_theta+1]
            for i_phi in range(N3):#{
                phi_min,phi_max = bins3[i_phi],bins3[i_phi+1]
                
                N['on'],List['on'] = len_sam_in_3d_bin(beam_on,
                                                       reco_1,Pmin,Pmax,
                                                       reco_2,cos_theta_min,cos_theta_max,
                                                       reco_3,phi_min,phi_max)
                                                       
                N['off'],List['off'] = len_sam_in_3d_bin(beam_off,
                                                         reco_1,Pmin,Pmax,
                                                         reco_2,cos_theta_min,cos_theta_max,
                                                         reco_3,phi_min,phi_max)

                N['generated'],List['generated'] = len_sam_in_3d_bin(generated_CC1p,
                                                                     true_1,Pmin,Pmax,
                                                                     true_2,cos_theta_min,cos_theta_max,
                                                                     true_3,phi_min,phi_max)

                generated_CC1p_in_kin_cuts = generated_CC1p[(np.abs(generated_CC1p['theta_12']-90)<delta_theta_12)
                                                            &(generated_CC1p['Pt']<Pt_max)
                                                            &(np.abs(np.abs(generated_CC1p['delta_phi'])-180.)<delta_Delta_phi)]

                N['gen. in kin. cuts'],List['gen. in kin. cuts'] = len_sam_in_3d_bin(generated_CC1p_in_kin_cuts,
                                                                                     true_1,Pmin,Pmax,
                                                                                     true_2,cos_theta_min,cos_theta_max,
                                                                                     true_3,phi_min,phi_max)
                N['CC1p'],List['CC1p'] = len_sam_in_3d_bin(selected_CC1p,
                                                           reco_1,Pmin,Pmax,
                                                           reco_2,cos_theta_min,cos_theta_max,
                                                           reco_3,phi_min,phi_max)


                N['overlay'],List['overlay'] = len_sam_in_3d_bin(overlay,
                                                                     reco_1,Pmin,Pmax,
                                                                     reco_2,cos_theta_min,cos_theta_max,
                                                                     reco_3,phi_min,phi_max)
                        

                Ngen = N['gen. in kin. cuts'] if do_in_kin_cuts else N['generated']
                if option=="CC1p truth":#{
                    N['CC1p'],List['CC1p'] = len_sam_in_3d_bin(selected_CC1p,
                                                                           true_1,Pmin,Pmax,
                                                                           true_2,cos_theta_min,cos_theta_max,
                                                                           true_3,phi_min,phi_max)
                #}
                eff, eff_err = get_eff(Ngen=Ngen , Nsel=N['CC1p'])
                
                # efficiency weight assigned to the event
                w,werr=0,0
                if eff>eff_cutoff:#{
                    power = float("1.e%d"%power_factor)
                    w = power/(eff*flux*Ntargets)
                    werr = power*eff_err/(eff*eff*flux*Ntargets)
                #}
                if debug:#{
                    print 'phi_min,phi_max:',phi_min,phi_max
                    print N['CC1p'],'CC1p',N['gen. in kin. cuts'],'gen. in kin. cuts'
                    print 'eff=',eff,'+/-',eff_err
                    print 'w=',w,'+/-',werr
                    print "List['on']:",len(List['on'])
                    print "List['off']:",len(List['off'])
                    print "List['CC1p']:",len(List['CC1p'])
                #}
                wname = ob_1 + ' weight'
                
                if option=="CC1p truth": wname = ob_1 + ' weight truth'
                
                beam_on.loc[List['on'],wname] = w
                beam_on.loc[List['on'],wname+' err'] = werr

                beam_off.loc[List['off'],wname] = w
                beam_off.loc[List['off'],wname+' err'] = werr

                selected_CC1p.loc[List['CC1p'],wname] = w
                selected_CC1p.loc[List['CC1p'],wname+' err'] = werr

                overlay.loc[List['overlay'],wname] = w
                overlay.loc[List['overlay'],wname+' err'] = werr
            #} i_phi
        #} i_cos_theta
    #} i_P
    if debug: print 'done.'
    return
# ----------------------------------------------------------




#
#
## ----------------------------------------------------------
## Sep-14, 2018
#def get_Xsec(h = None ,afro_genie_dict=dict(),ylim_P=(0,9),ylim_cos_theta=(0,9),ylim_phi=(0,0.1)
#             ,do_add_genie_models=True
#             ,ob_1='Pmu',ob_2='cos(theta(mu))',ob_3='phi(mu)'
#             ,true_1='truth_Pmu',true_2='truth_Pmu_cos_theta',true_3='truth_Pmu_phi'
#             ):#{
#    get_integrated_Xsec(h,bins1,bins2,bins3,N1,N2,N3);
#    if do_add_genie_models:#{
#        for gname,ls in zip(['nominal','hA2015','hA_SRC'],['-','--','-.']):
#            afro_genie_CC1p = afro_genie_dict[gname]
#            afro_Xsec,afro_Xsec_err = compute_Xsec(Non=len(afro_genie_CC1p), Non_err=np.sqrt(len(afro_genie_CC1p)))
#            afro_Xsec,afro_Xsec_err = afro_Xsec*4.908e19/4.9e20,afro_Xsec_err*4.908e19/4.9e20
#            print gname,'afro genie Xsec: %.2f +/- %.2f'%(afro_Xsec,afro_Xsec_err),'e-38 cm2'
#    #}
#    fig=plt.figure(figsize=(28,8))
#    observable = ob_1
#    bins,mid,bin_width,vlabel,xlabel,units = get_labels(observable=observable)
#    Xsec_1d,Xsec_1d_err,mc_Xsec_1d,mc_Xsec_1d_err = np.zeros(N1),np.zeros(N1),np.zeros(N1),np.zeros(N1)
#    for i_P in range(N1):#{
#        Xsec_1d_err_sq_sum,mc_Xsec_1d_err_sq_sum = 0,0
#        for i_cos_theta in range(N2):#{
#            cos_theta_bin_width = bins2[i_cos_theta+1] - bins2[i_cos_theta]
#            for i_phi in range(N3):#{
#                phi_bin_width = bins3[i_phi+1] - bins3[i_phi]
#                Xsec_1d[i_P] += h['Xsec'][i_P][i_cos_theta][i_phi] * cos_theta_bin_width * phi_bin_width
#                Xsec_1d_err_sq_sum += np.square(h['Xsec err'][i_P][i_cos_theta][i_phi] * cos_theta_bin_width * phi_bin_width)
#                mc_Xsec_1d[i_P] += h['mc-Xsec'][i_P][i_cos_theta][i_phi] * cos_theta_bin_width * phi_bin_width
#                mc_Xsec_1d_err_sq_sum += np.square(h['mc-Xsec err'][i_P][i_cos_theta][i_phi] * cos_theta_bin_width * phi_bin_width)
#        #}
##}
#Xsec_1d_err[i_P] = np.sqrt(Xsec_1d_err_sq_sum)
#    mc_Xsec_1d_err[i_P] = np.sqrt(mc_Xsec_1d_err_sq_sum)
#    #}
#    ax = fig.add_subplot(1,3,1)
#    ax.bar( x=mid , height=2*mc_Xsec_1d_err, bottom=mc_Xsec_1d-mc_Xsec_1d_err, width=bin_width, color=Colors['CC1p'])
#    #     ax.bar( x=mid , height=2*genie_Xsec_err, bottom=genie_Xsec-genie_Xsec_err, width=bin_width, color='black')
#    ax.errorbar( x=mid , xerr=0.5*bin_width, y=Xsec_1d, yerr=Xsec_1d_err , fmt='o', markersize=10
#                , color=Colors['beam on'], capsize=1, capthick=3, label='data')
#                if do_add_genie_models:
#                    for gname,ls in zip(['nominal','hA2015','hA_SRC'],['-','--','-.']):
#                        h_genie,_ = np.histogram(afro_genie_dict[gname][true_1],bins=bins)
#                            h_genie_err = np.sqrt(h_genie)
#                                h_genie,h_genie_err = h_genie*4.908e19/4.9e20, h_genie_err*4.908e19/4.9e20
#                                    genie_Xsec,genie_Xsec_err = np.zeros((len(bins)-1)), np.zeros((len(bins)-1))
#                                        for i in range(len(bins)-1):
#                                            genie_Xsec[i],genie_Xsec_err[i] = compute_Xsec(Non=h_genie[i], Non_err=h_genie_err[i] ,eff=1,eff_err=0,B=0, bin_width=bin_width )
#                                                mystep(x=mid ,dx=bin_width, y=genie_Xsec, y_width=genie_Xsec_err, color='black',linestyle=ls,linewidth=3,label=r'genie ('+gname+')')
#                                            set_axes(ax,xlabel
#                                                     ,y_label=(r'$\frac{d\sigma}{d'+vlabel+'}$' +r'$\left[10^{-38} \frac{cm^{2}}{(%s)}\right]$'%units)
#                                                     ,ylim=ylim_P)
##-------------------------------
#observable = ob_2
#    bins,mid,bin_width,vlabel,xlabel,units = get_labels(observable=observable)
#    Xsec_1d,Xsec_1d_err,mc_Xsec_1d,mc_Xsec_1d_err = np.zeros(N2),np.zeros(N2),np.zeros(N2),np.zeros(N2)
#    for i_cos_theta in range(N2):#{
#        Xsec_1d_err_sq_sum,mc_Xsec_1d_err_sq_sum = 0,0
#        for i_P in range(N1):#{
#            P_bin_width = bins1[i_P+1] - bins1[i_P]
#            for i_phi in range(N3):#{
#                phi_bin_width = bins3[i_phi+1] - bins3[i_phi]
#                Xsec_1d[i_cos_theta] += h['Xsec'][i_P][i_cos_theta][i_phi] * P_bin_width * phi_bin_width
#                Xsec_1d_err_sq_sum += np.square(h['Xsec err'][i_P][i_cos_theta][i_phi] * P_bin_width * phi_bin_width)
#                mc_Xsec_1d[i_cos_theta] += h['mc-Xsec'][i_P][i_cos_theta][i_phi] * P_bin_width * phi_bin_width
#                mc_Xsec_1d_err_sq_sum += np.square(h['mc-Xsec err'][i_P][i_cos_theta][i_phi] * P_bin_width * phi_bin_width)
#        #}
#        #}
#        Xsec_1d_err[i_cos_theta] = np.sqrt(Xsec_1d_err_sq_sum)
#    mc_Xsec_1d_err[i_cos_theta] = np.sqrt(mc_Xsec_1d_err_sq_sum)
##}
#ax = fig.add_subplot(1,3,2)
#    ax.bar( x=mid , height=2*mc_Xsec_1d_err, bottom=mc_Xsec_1d-mc_Xsec_1d_err, width=bin_width, color=Colors['CC1p'])
#    # ax.bar( x=mid , height=2*genie_Xsec_err, bottom=genie_Xsec-genie_Xsec_err, width=bin_width, color='black')
#    ax.errorbar( x=mid , xerr=0.5*bin_width, y=Xsec_1d, yerr=Xsec_1d_err , fmt='o', markersize=10
#                , color=Colors['beam on'], capsize=1, capthick=3, label='data')
#                if do_add_genie_models:
#                    for gname,ls in zip(['nominal','hA2015','hA_SRC'],['-','--','-.']):
#                        h_genie,_ = np.histogram(afro_genie_dict[gname][true_2],bins=bins)
#                            h_genie_err = np.sqrt(h_genie)
#                                h_genie,h_genie_err = h_genie*4.908e19/4.9e20, h_genie_err*4.908e19/4.9e20
#                                    genie_Xsec,genie_Xsec_err = np.zeros((len(bins)-1)), np.zeros((len(bins)-1))
#                                        for i in range(len(bins)-1):
#                                            genie_Xsec[i],genie_Xsec_err[i] = compute_Xsec(Non=h_genie[i], Non_err=h_genie_err[i] ,eff=1,eff_err=0,B=0, bin_width=bin_width )
#                                                mystep(x=mid ,dx=bin_width, y=genie_Xsec, y_width=genie_Xsec_err, color='black',linestyle=ls,linewidth=3,label=r'genie ('+gname+')')
#                                            set_axes(ax,xlabel
#                                                     ,y_label=(r'$\frac{d\sigma}{d'+vlabel+'}$'+r'[$10^{-38}$ cm$^{2}$]'),ylim=ylim_cos_theta)
#    # #-------------------------------
#    observable = ob_3
#    bins,mid,bin_width,vlabel,xlabel,units = get_labels(observable=observable)
#    
#    Xsec_1d,Xsec_1d_err,mc_Xsec_1d,mc_Xsec_1d_err = np.zeros(N3),np.zeros(N3),np.zeros(N3),np.zeros(N3)
#    for i_phi in range(N3):#{
#        Xsec_1d_err_sq_sum,mc_Xsec_1d_err_sq_sum = 0,0
#        for i_P in range(N1):#{
#            P_bin_width = bins1[i_P+1] - bins1[i_P]
#            for i_cos_theta in range(N2):#{
#                cos_theta_bin_width = bins2[i_cos_theta+1] - bins2[i_cos_theta]
#                Xsec_1d[i_phi] += h['Xsec'][i_P][i_cos_theta][i_phi] * P_bin_width * cos_theta_bin_width
#                Xsec_1d_err_sq_sum += np.square(h['Xsec err'][i_P][i_cos_theta][i_phi] * P_bin_width * cos_theta_bin_width)
#                mc_Xsec_1d[i_phi] += h['mc-Xsec'][i_P][i_cos_theta][i_phi] * P_bin_width * cos_theta_bin_width
#                mc_Xsec_1d_err_sq_sum += np.square(h['mc-Xsec err'][i_P][i_cos_theta][i_phi] * P_bin_width * cos_theta_bin_width)
#        #}
#        #}
#        Xsec_1d_err[i_phi] = np.sqrt(Xsec_1d_err_sq_sum)
#    mc_Xsec_1d_err[i_phi] = np.sqrt(mc_Xsec_1d_err_sq_sum)
#    #}
#    ax = fig.add_subplot(1,3,3)
#    ax.bar( x=mid , height=2*mc_Xsec_1d_err, bottom=mc_Xsec_1d-mc_Xsec_1d_err, width=bin_width, color=Colors['CC1p'])
#    #ax.bar( x=mid , height=2*genie_Xsec_err, bottom=genie_Xsec-genie_Xsec_err, width=bin_width, color='black')
#    ax.errorbar( x=mid , xerr=0.5*bin_width, y=Xsec_1d, yerr=Xsec_1d_err , fmt='o', markersize=10
#                , color=Colors['beam on'], capsize=1, capthick=3, label='data')
#    if do_add_genie_models:
#    for gname,ls in zip(['nominal','hA2015','hA_SRC'],['-','--','-.']):
#        h_genie,_ = np.histogram(180./np.pi*afro_genie_dict[gname][true_3],bins=bins)
#        h_genie_err = np.sqrt(h_genie)
#        h_genie,h_genie_err = h_genie*4.908e19/4.9e20, h_genie_err*4.908e19/4.9e20
#        genie_Xsec,genie_Xsec_err = np.zeros((len(bins)-1)), np.zeros((len(bins)-1))
#    for i in range(len(bins)-1):
#    genie_Xsec[i],genie_Xsec_err[i] = compute_Xsec(Non=h_genie[i], Non_err=h_genie_err[i] ,eff=1,eff_err=0,B=0, bin_width=bin_width )
#    mystep(x=mid ,dx=bin_width, y=genie_Xsec, y_width=genie_Xsec_err, color='black',linestyle=ls,linewidth=3,label=r'genie ('+gname+')')
#    set_axes(ax,xlabel
#                                                     ,y_label=(r'$\frac{d\sigma}{d'+vlabel+'}$'+r'$\left[10^{-38} \frac{cm^{2}}{(%s)}\right]$'%units)
#                                                     ,do_add_legend=False,ylim=ylim_phi)
#plt.tight_layout()
#    print 'done.'
##}
## ----------------------------------------------------------




# ----------------------------------------------------------
# Sep-14, 2018
# integrated cross-section
def get_integrated_Xsec(h,bins1,bins2,bins3,N1,N2,N3):
    Xsec_integrated,Xsec_integrated_err_sq_sum = 0,0
    mc_Xsec_integrated,mc_Xsec_integrated_err_sq_sum = 0,0
    for i_P in range(N1):#{
        P_bin_width = bins1[i_P+1] - bins1[i_P]
        for i_cos_theta in range(N2):#{
            cos_theta_bin_width = bins2[i_cos_theta+1] - bins2[i_cos_theta]
            for i_phi in range(N3):#{
                phi_bin_width = bins3[i_phi+1] - bins3[i_phi]
                
                bin_width = P_bin_width * cos_theta_bin_width * phi_bin_width
                Xsec_integrated += h['Xsec'][i_P][i_cos_theta][i_phi]*bin_width
                Xsec_integrated_err_sq_sum += np.square(h['Xsec err'][i_P][i_cos_theta][i_phi]*bin_width)
                mc_Xsec_integrated += h['mc-Xsec'][i_P][i_cos_theta][i_phi]*bin_width
                mc_Xsec_integrated_err_sq_sum += np.square(h['mc-Xsec err'][i_P][i_cos_theta][i_phi]*bin_width)

    Xsec_integrated_err = np.sqrt(Xsec_integrated_err_sq_sum)
    print 'Xsec_integrated: %.2f +/- %.2f'%(Xsec_integrated,Xsec_integrated_err),'e-38 cm2'
    mc_Xsec_integrated_err = np.sqrt(mc_Xsec_integrated_err_sq_sum)
    print 'mc_Xsec_integrated: %.2f +/- %.2f'%(mc_Xsec_integrated,mc_Xsec_integrated_err),'e-38 cm2'
    return Xsec_integrated,Xsec_integrated_err,mc_Xsec_integrated,mc_Xsec_integrated_err
# ----------------------------------------------------------


# ----------------------------------------------------------
# Sep-04, 2018 (edited Sep-22)
def get_labels(observable=''):
    var = Vars[observable]
    bins=Bins[observable]; vlabel=vlabels[observable]; Vlabel=Vlabels[observable]; units=Units[observable]
    xlabel=Vlabel+' ['+units+']' if units is not None else Vlabel
    mid = 0.5*(bins[1:]+bins[:-1]);
#    bin_width = (mid[1]-mid[0])
    bin_width = bins[1:] - bins[:-1]
    return var,bins,mid,bin_width,vlabel,xlabel,units
# ----------------------------------------------------------


# ----------------------------------------------------------
# Sep-02, 2018 (edited Sep-16)
def len_sam_in_3d_bin(sam,
                  Pvar,Pmin,Pmax,
                  cos_theta_var,cos_theta_min,cos_theta_max,
                  phi_var,phi_min,phi_max):
    
    list_of_indices = sam.index[(Pmin<=sam[Pvar])&(sam[Pvar] < Pmax)
                                &
                                (cos_theta_min<=sam[cos_theta_var])&(sam[cos_theta_var] < cos_theta_max)
                                &
                                (phi_min <= 180./np.pi*sam[phi_var])&(180./np.pi*sam[phi_var] < phi_max)
                                ].tolist()
                                
    return len(list_of_indices),list_of_indices
# ----------------------------------------------------------


# ----------------------------------------------------------
# Sep-02, 2018 (edited Sep-16)
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
        Xsec        cross-section in units of (10^-38) cm2 / bin_units
        Xsec_err    cross-section uncertainty in units of (10^-38) cm2 / bin_units
        '''
    
    if eff<eff_cutoff:  return 0,0
    

    num = Non - Noff - B
    den = eff * Ntargets * flux * bin_width
    Xsec = np.max( [0, num/den] )

    num_err = np.sqrt( np.square(Non_err) + np.square(Noff_err) + np.square(B_err) )
    den_err = den * np.sqrt( np.square(eff_err/eff) + np.square(Ntargets_err/Ntargets) + np.square(flux_err/flux) )
    Xsec_err = Xsec * np.max( [0 , np.sqrt(  (np.square(num_err/num) if num>0 else 0) + (np.square(den_err/den) if den>0 else 0) ) ] )
    
    return Xsec*1e38, Xsec_err*1e38
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
    h['Xsec'],herr['Xsec'] = h['Xsec']*1e38,herr['Xsec']*1e38
                           

    h['mc'] = h['mc']/(h['eff']*Ntargets*flux*bin_width)
    
    herr['mc'] = np.sqrt( np.square(herr['mc']/(h['eff']*Ntargets*flux*bin_width))
                                                +np.square(h['mc']*herr['eff']/(h['eff']*h['eff']*Ntargets*flux*bin_width))
                                                +np.square(h['mc']*Ntargets_err/(h['eff']*Ntargets*Ntargets*flux*bin_width))
                                                +np.square(h['mc']*flux_err/(h['eff']*Ntargets*flux*flux*bin_width)))
    h['mc'],herr['mc'] = h['mc']*1e38*mc_scale_factor,herr['mc']*1e38*mc_scale_factor
    
    # draw the cross-section
    if ax is None:#{
        fig,ax = plt.subplots(figsize=(9.708,6))
    #}
    ax.bar( x=mid , height=2*herr['mc'], bottom=h['mc']-herr['mc'], width=2*bin_width, color=Colors['CC1p'], label=r'mc ($\Delta$stat.)')
    ax.errorbar( x=mid , xerr=bin_width, y=h['Xsec'], yerr=herr['Xsec'] , fmt='o', markersize=10
                    , color=Colors['beam on'], capsize=1, capthick=3, label='data')
    set_axes(ax,xlabel
             ,y_label=(r'$\frac{d\sigma}{d'+vlabel+'}$'
                       +r' $\left[\times 10^{-38}\right]$ '
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
    return sam[ (Limits['Pmu'][0] <= sam[varPmu]) & (sam[varPmu] < Limits['Pmu'][1])
               & (Limits['cos(theta(mu))'][0] <= sam[varPmu_cos_theta]) & (sam[varPmu_cos_theta] < Limits['cos(theta(mu))'][1])
               & (Limits['Pp'][0] <= sam[varPp])  & (sam[varPp] < Limits['Pp'][1])
               & (Limits['cos(theta(p))'][0] <= sam[varPp_cos_theta]) & (sam[varPp_cos_theta] < Limits['cos(theta(p))'][1])
               ]
#}
#---------------------------------------------------------------------------------------------



# ----------------------------------------------------------
# Aug-27, 2018 (last edit Oct-09, 2018)
def load_mc_and_data(extra_name='',debug=0
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
    overlay_prefix = Paths['selected events'] + versions['Overlay'] + '_' + versions['overlay date'] + '_' + extra_name
    data_prefix = Paths['selected events'] + versions['beam on'] + '_' + versions['data date'] + '_' + extra_name
    selected_overlay=dict()
    
    cuts_order  = ['no cut','Chi2Proton','Nflashes','MatchedFlash','length'
                   ,'non-collinearity','vertex activity','delta phi','Pt & delta phi']
        
    if os.path.isfile( overlay_prefix + 'selected_CC1p.csv'):#{
        print 'found selected overlay files from '+extra_name+', loading them...'
        for pair_type in pair_types:#{
            selected_overlay[pair_type]=pd.read_csv(overlay_prefix+'selected_'+pair_type+'.csv')
        #}
        selected_overlay_concat = pd.read_csv(overlay_prefix+'selected_overlay.csv')
    #}
    else:#{
        print 'did not find selected overlay files from '+extra_name+', so creating it...'
        OverlaySamples = load_samples(date=versions['overlay date'],filename=versions['Overlay']+'_'+versions['overlay date']+'_vertices')
        reducedOverlay,pureffOverlay,pureffNumbers = apply_cuts_to_overlay(OverlaySamples=OverlaySamples, cuts_order=cuts_order,debug=debug
                                                                           ,minPEcut = minPEcut
                                                                           ,maxdYZcut = maxdYZcut
                                                                           ,delta_theta_12 = delta_theta_12
                                                                           ,r_max_RdQ_CC1p = r_max_RdQ_CC1p
                                                                           ,delta_Delta_phi=delta_Delta_phi
                                                                           ,Pt_max=Pt_max
                                                                           ,Chi2Proton_muCandidate_min=Chi2Proton_muCandidate_min
                                                                           ,Chi2Proton_pCandidate_max=Chi2Proton_pCandidate_max)
        print 'applied cuts to overlay \n'
        for pair_type in pair_types:#{
            selected_overlay[pair_type] = sample_in_limits(sam=reducedOverlay['Pt & delta phi'][pair_type])
            outcsvname = overlay_prefix+'selected_'+pair_type+'.csv'
            selected_overlay[pair_type].to_csv(outcsvname)
            print 'saved',len(selected_overlay[pair_type]),'selected',pair_type,'events to','selected_'+pair_type+'.csv'
        #}
        selected_overlay_concat = pd.concat([selected_overlay['1mu-1p'],selected_overlay['cosmic'],selected_overlay['other-pairs']])
        outcsvname = overlay_prefix+'selected_overlay.csv'
        selected_overlay_concat.to_csv(outcsvname)
        print 'saved',len(selected_overlay_concat),'selected overlay to',outcsvname
        # overlay scaling
        summary = pd.read_csv('/Users/erezcohen/Desktop/uBoone/CCQEanalysis/csvFiles/summary/'
                              +versions['overlay date']+'/'
                              +versions['Overlay']+'_'+versions['overlay date']+'_summary.csv')
        Nevents['OnBeam POT']   = 4.908e+19
        Nevents['overlay']      = np.sum(summary.Nevents)
        Nevents['overlay POT']  = np.sum(summary.POT)
        Nevents['f(POT)']       = Nevents['OnBeam POT']/Nevents['overlay POT']
        print "Nevents['f(POT)']:",Nevents['f(POT)']
    #}
    selected_CC1p = selected_overlay['CC1p']
    print len(selected_CC1p),'selected CC1p events overlay'
    # ----------------------------------------------------------
    ## (2) DATA
    if os.path.isfile(data_prefix+'selected_beam_on.csv'):#{
        selected_beam_on = pd.read_csv(data_prefix+'selected_beam_on.csv')
        selected_beam_off = pd.read_csv(data_prefix+'selected_beam_off.csv')
        print 'found ',len(selected_beam_on),'selected on beam and',len(selected_beam_off),'beam off events...'
    #}
    else:#{
        print 'found selected on beam events and there was no file there...'
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
    if os.path.isfile(overlay_prefix + 'selected_genie_CC1p.csv'):#{
        print 'found selected genie CC1p...'
        genie_CC1p = pd.read_csv(overlay_prefix+'selected_genie_CC1p.csv')
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
        outcsvname = overlay_prefix+'selected_genie_CC1p.csv'
        genie_CC1p.to_csv(outcsvname)
        print 'saved %d'%len(genie_CC1p),'CC1p events in genie_CC1p to',outcsvname
    #}
    return selected_overlay,selected_overlay_concat,selected_CC1p,genie_CC1p,selected_beam_on,selected_beam_off
#}
# ----------------------------------------------------------

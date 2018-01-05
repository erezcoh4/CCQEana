import numpy as np

M = 0.939 # GeV
m = 0.1056 # GeV
m_W = 80.42 # GeV # +/- 0.06
g_A = 1.267 # +/- 0.004
cos_theta_C = 0.9755 # +/- 0.004
G_F = 1.1663787 # GeV**-2 # +/- 0.0000006
mu_p = 2.7928474 # proton dipole moment
mu_n = -1.913 # neutron dipole moment
M2 = M * M
m2 = m * m

def Q2label(subscript=''):
    return r'$Q^2_{'+subscript+'}$ (GeV/c)$^2$'


def plot_Q2_Ev_distribution_in_reQ2_bins(histo_array=None,histtitle=''):
    fig = plt.figure(figsize=(20,20))
    for i,recQ2bin in zip(range(len(Q2bins)-1),Q2bins):
        ax = fig.add_subplot(len(Q2bins)/3,3,i+1)
        sns.heatmap(histo_array[i].T,cmap='hot_r'#                 ,annot=True#,fmt='%.1f'            
                    ,xticklabels=['%.1f'%largeQ2bins[j] if j%3==0 else '' for j in range(len(largeQ2bins)-1)]
                    ,yticklabels=['%.1f'%largeEvbins[j] if j%3==0 else '' for j in range(len(largeEvbins)-1)]
               )
        ax.invert_yaxis()
        set_axes(ax
             ,title=r"%s $%.2f<Q^{2}_{rec.}<%.2f$"%(histtitle,Q2bins[i],Q2bins[i+1])
             ,x_label=r'$Q^2_{gen.}$ (GeV/c)$^2$' if i>=(len(Q2bins)-1-3) else ''
             ,y_label=r'$E^\nu_{gen.}$ (GeV)' if i%3==0 else ''
             ,fontsize=25,do_add_grid=True,do_add_legend=True
            )
        if i%3>0: ax.yaxis.set_major_formatter(ticker.NullFormatter())    
    plt.tight_layout()        
    
    
def fill_2d_histogram( xbins , ybins , samlpe , varx , vary ):
    '''
    return: 
            histo: np.array((len(xbins),len(ybins))) 
    '''
    histo = np.zeros((len(xbins),len(ybins)))
    for i_xbin,xbin in zip(range(len(xbins)-1),xbins):
        xmin,xmax = xbins[i_xbin],xbins[i_xbin+1]
        for j_ybin,ybin in zip(range(len(ybins)-1),ybins):
            ymin,ymax = ybins[j_ybin],ybins[j_ybin+1]

            histo[i_xbin][j_ybin] = len(samlpe[(samlpe[varx]>xmin)
                                                  &(samlpe[varx]<xmax)
                                                  &(samlpe[vary]>ymin)
                                                  &(samlpe[vary]<ymax)]) 
    return histo

def normalize_yield(y,yerr,label=''):
    total_yield = np.sum(y)
    y = np.array([100*np.float(y[i])/total_yield for i in range(len(y))])
    yerr = np.array([100*np.float(yerr[i])/total_yield for i in range(len(yerr))])
    print label
    print ['%.1f+/-%.1f'%(y[i],yerr[i])+'%' for i in range(len(yerr))]
    return y,yerr



def tau( Q2 = 0.0 ):
    return (Q2 / (4.*M2))

def G_d( Q2 = 0.0 , Lambda2 = 0.71 ):    
    return np.power((1. + Q2/Lambda2),-2)

def G_A_dipole( Q2 = 0.0 , mA=1.0):  
    '''
    Here mA is coming to play as the axial mass of the dipole form-factor
    '''
    return np.power((1. + Q2/np.square(mA)),-2)

def G_E_n_Bradford( Q2 = 0.0 ):
    
    a = np.array([0 , 1.25 , 1.30 ])
    a_err = np.array([0 , 0.368 , 1.99 ])
    
    b = np.array([ 0 , -9.86 , 305 , -758 , 802])
    b_err = np.array([ 0 , 6.46 , 28.6 , 77.5 , 156])

    numerator = a[0]*(np.power(tau(Q2),0)) + a[1]*(np.power(tau(Q2),1)) + a[2]*(np.power(tau(Q2),2))
    denominator = 1.0 + b[1]*(np.power(tau(Q2),1)) + b[2]*(np.power(tau(Q2),2)) + b[3]*(np.power(tau(Q2),3)) + b[4]*(np.power(tau(Q2),4))


    return (numerator/denominator)

def G_E_p_Bradford( Q2 = 0.0 ):
    
    a = np.array([1. , -0.0578 ])
    a_err = np.array([0 , 0.166 ])
    
    b = np.array([ 0 , 11.1 , 13.6 , 33.0])
    b_err = np.array([ 0 , 0.217 , 1.39 , 8.95])

    numerator = a[0]*(np.power(tau(Q2),0)) + a[1]*(np.power(tau(Q2),1))
    denominator = 1.0 + b[1]*(np.power(tau(Q2),1)) + b[2]*(np.power(tau(Q2),2)) + b[3]*(np.power(tau(Q2),3))

    return (numerator/denominator)

def G_E_Bradford( Q2 = 0.0 ):
    return 0.5*( G_E_p_Bradford(Q2) - G_E_n_Bradford(Q2) )

def G_M_n_Bradford( Q2 = 0.0 ):

    # Bradford gives G_M_n/mu_n in their paper
    a = np.array([ 1 , 1.81 ])
    a_err = np.array([ 0 , 0.402 ])
    
    b = np.array([ 0 , 14.1 , 20.7 , 68.7])
    b_err = np.array([ 0 , 0.597 , 2.55 , 14.1])

    numerator = a[0]*(np.power(tau(Q2),0)) + a[1]*(np.power(tau(Q2),1))
    denominator = 1.0 + b[1]*(np.power(tau(Q2),1)) + b[2]*(np.power(tau(Q2),2)) + b[3]*(np.power(tau(Q2),3))

    return mu_n*(numerator/denominator)

def G_M_p_Bradford( Q2 = 0.0 ):
    a = np.array([ 1 , 0.150 ])
    a_err = np.array([ 0 , 0.0312 ])
    
    b = np.array([ 0 , 11.1 , 19.6 , 7.54])
    b_err = np.array([ 0 , 0.103 , 0.281 , 0.967])

    numerator = a[0]*(np.power(tau(Q2),0)) + a[1]*(np.power(tau(Q2),1))
    denominator = 1.0 + b[1]*(np.power(tau(Q2),1)) + b[2]*(np.power(tau(Q2),2)) + b[3]*(np.power(tau(Q2),3))

    return mu_p*(numerator/denominator)

def G_M_Bradford( Q2 = 0.0 ):
    return 0.5*( G_M_p_Bradford(Q2) - G_M_n_Bradford(Q2) )

def F1_p_Bradford( Q2 = 0.0 ):
    return (1./(1.0+tau(Q2))) * ( tau(Q2)*G_M_p_Bradford(Q2) + G_E_p_Bradford(Q2) )
def F2_p_Bradford( Q2 = 0.0 ):
    return (1./(1.0+tau(Q2))) * ( G_M_p_Bradford(Q2) - G_E_p_Bradford(Q2) )
def F1_n_Bradford( Q2 = 0.0 ):
    return (1./(1.0+tau(Q2))) * ( tau(Q2)*G_M_n_Bradford(Q2) + G_E_n_Bradford(Q2) )
def F2_n_Bradford( Q2 = 0.0 ):
    return (1./(1.0+tau(Q2))) * ( G_M_n_Bradford(Q2) - G_E_n_Bradford(Q2) )


def F1V_Bradford( Q2 = 0.0 ):
    return (1./(1.0+tau(Q2))) * ( tau(Q2)*G_M_Bradford(Q2) + G_E_Bradford(Q2) )
def F2V_Bradford( Q2 = 0.0 ):
    return (1./(1.0+tau(Q2))) * ( G_M_Bradford(Q2) - G_E_Bradford(Q2) )

def G_E_p_Bernauer( Q2 = 0.0 ):
    a = np.array([ 1.0 , -3.3686, 14.5606 , -88.1912 
                  , 453.6244 , -1638.7911 , 3980.7174 
                  , -6312.6333 , 6222.3646 , -3443.2251 
                  , 814.4112 ])
    Q2_powers = np.array([np.power(Q2,n) for n in range(len(a))])
    return  a[0]*(Q2**0) + a[1]*np.power(Q2,1) + a[2]*np.power(Q2,2) + a[3]*np.power(Q2,3) + a[4]*np.power(Q2,4) + a[5]*np.power(Q2,5)
def G_M_p_Bernauer( Q2 = 0.0 ):
    a = np.array([1.0 ,-2.5952 , 1.0222 , 23.4945, -93.0372
                  ,140.7984, -0.3656 , -305.6759 , 444.6251 
                  , -273.6688 , 64.5811 ])
    a = np.array([ 1.0 ,-2.5952])
    return  a[0]*(Q2**0) + a[1]*np.power(Q2,1) + a[2]*np.power(Q2,2) + a[3]*np.power(Q2,3)
    Q2_powers = np.array([np.power(Q2,n) for n in range(len(a))])
    
    
    
G_F2 = G_F * G_F
cos2_theta_C = cos_theta_C*cos_theta_C
m_W2 = m_W * m_W

# Kinematical Factor
def KinFactor( E=0 , Q2=0.0):
    F = (M2*G_F2*cos2_theta_C)/(8*np.pi)
    #     return ( F/np.square(E) ) * np.square( m_W2/(m_W2+Q2) )  # F is just a Fermi-factor, independent for the experiment
    return ( 1./np.square(E) ) * np.square( m_W2/(m_W2+Q2) ) 


# xi
def xi( E=0 , Q2=0.0 ):
    return ( (4*M*E - Q2 - m2)/(M2))


# A
def A( Q2=0.0 , mA=1.0 ):
    t = tau(Q2)
    G_A = G_A_dipole(Q2=Q2,mA=mA)
    G_E_V = G_E_Bradford(Q2)
    G_M_V = G_M_Bradford(Q2)
    #     I = np.ones(len(t))    
    return (4.*t*( (1+t)*np.square(G_A) - 4*np.square(G_E_V) + 4*t*np.square(G_M_V) ))
# The last two terms in A are usually discarded in experimental analyses as they are suppressed by powers of (m/M)2 \simeq 1% 


# B
def B( Q2=0.0 , mA=1.0 ):
    t = tau(Q2)
    G_A = G_A_dipole(Q2=Q2,mA=mA)
    G_M_V = G_M_Bradford(Q2)
    return (8 * t * G_A * G_M_V)


# C
def C( Q2=0.0 , mA=1.0 ):
    t = tau(Q2)
    G_A = G_A_dipole(Q2=Q2,mA=mA)
    F_1_V = F1V_Bradford(Q2)
    F_2_V = F2V_Bradford(Q2)
    return ( 0.25*np.square(G_A) + np.square(F_1_V) + t*np.square(F_2_V) )

def CCelasticXsec( Q2 , Ev, mA ):
    return (KinFactor(Ev , Q2) * ( A(Q2,mA) + xi(Ev,Q2)*B(Q2,mA) + np.square(xi(Ev,Q2))*C(Q2,mA) ))

def N_CCelasticXsec( Q2 , Ev, mA , NormFact ):
    # normalization factor for the data mutiplies the CCelastic Xsec
    return (NormFact*KinFactor(Ev , Q2) * ( A(Q2,mA) + xi(Ev,Q2)*B(Q2,mA) + np.square(xi(Ev,Q2))*C(Q2,mA) ))    
/**
 * \file GenieFile.h
 *
 * \ingroup CCQEPackage
 *
 * \brief Class def header for a class GenieFile
 *
 * @author erezcohen
 */

/** \addtogroup CCQEPackage
 
 @{*/
#ifndef GENIEFILE_H
#define GENIEFILE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>

#include "TFile.h"
#include "TTree.h"
#include "Rtypes.h"
#include "TRandom3.h"
#include "GENIEinteraction.h"
#include "TGraph.h"
#include "TSpline.h"

#define NMAX 40
#define mAr40 39.96238312251*0.931494
/**
 \class GenieFile
 User defined class GenieFile ... these comments are used to generate
 doxygen documentation!
 */
class GenieFile{
    
public:
    
    GenieFile() = default;
    GenieFile( TString fPath
              ,TString fRootFileName
              ,TString fRootTreeName
              ,int fdebug
              ,TString fAccMapPath
              ,TString fPmuThetaAccMapName,TString fPpThetaAccMapName,TString fQ2ThetaAccMapName);
    

    
    // SETters
    bool                         SetInTree ();
    void         SetPmuThetaAcceptanceMaps (TString fAccMapPath, TString fAccMapName);
    void          SetPpThetaAcceptanceMaps (TString fAccMapPath, TString fAccMapName);
    void               SetQ2AcceptanceMaps (TString fAccMapPath, TString fAccMapName);
    void                 SetQ2_gen_rec_map (TString fMapPath,TString fMapName);
    void        SetQ2_rec_1d_Probabilities ();
    void                 SetVertexPosition (double x,double y, double z)                    {vertex_position = TVector3(x,y,z);};

    // GETters
    int                         GetNevents () const {return (int)GenieTree->GetEntries();};
    bool                  GetCC_Np_200MeVc () const {return IsCC_Np_200MeVc;};
    bool                  GetCC_1p_200MeVc () const {return IsCC_1p_200MeVc;};
    bool              GetCC_1p_200MeVc_0pi () const {return IsCC_1p_200MeVc_0pi;};
    TVector3             GetVertexPosition () const {return vertex_position;};
    
    
    // funcionality
    bool                         ReadEvent (int fi_event);
    bool                         HeaderCSV ();
    bool                       StreamToCSV ();
    bool                             Print ();
    bool                       SetTopology ();
    bool                        Initialize ();
    bool                            EndJob ();
    void               MimicDetectorVolume ();
    void             ProjectMuonTrajectory ();
    void                 CutMuonTrajectory ();
    void           ProjectProtonTrajectory ();
    void               CutProtonTrajectory ();
    double   LinePlaneIntersectionDistance (TVector3 l,TVector3 l0,TVector3 p0,TVector3 n,double epsilon=0.0001);
    void                 SetRecoKinematics ();
    int             SampleFromDistribution ( std::vector<double> probabilities );

    std::vector<double>              Read1dArrayFromFile (TString filename);
    std::vector<std::vector<double>> Read2dArrayFromFile (TString filename);
    void                             SetMicroBooNEWeights ();
    int                                     FindWhichBin ( double x, std::vector<double> bins );
    float                                 LightConeAlpha ( TLorentzVector v ){ return (v.E() - v.Pz())/(mAr40/40.); };
    
    
    
    
    
    
    
    int     debug=0;
    void Debug(Int_t verobosity_level, const char* format) // base function
    {
        if ( debug < verobosity_level ) return;
        std::cout << format << std::endl;
    }
    template<typename T, typename... Targs>
    void Debug(Int_t verobosity_level, const char* format, T value, Targs... Fargs) // recursive variadic function
    {
        if ( debug < verobosity_level ) return;
        for ( ; *format != '\0'; format++ ) {
            if ( *format == '%' ) {
                std::cout << value << " ";
                Debug(verobosity_level, format+1, Fargs...); // recursive call
                return;
            }
            std::cout << *format;
        }
        std::cout << endl;
    }
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    

    
    // stopping range in LAr for muons and protons
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    void SetLArTools(){
        
        cout << "SetLArTools()..." << endl;
        double E;
        // range in gr/cm2
        Float_t Range_grampercm_proton[31] =
        {1.887E-1, 3.823E-1, 6.335E-1, 1.296, 2.159, 7.375, 1.092E1,
            2.215E1, 3.627E1, 5.282E1, 7.144E1, 9.184E1, 1.138E2,
            1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2, 2.681E2,
            2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2,
            7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3
        };
        Float_t KE_MeV_proton[31] =
        {10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300, 350, 400,
            450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
            1500, 2000, 2500, 3000, 4000, 5000
        };
        
        
        double LAr_density = 1.396; //g/cm^3
        
        
        Float_t Range_grampercm_proton_LAr[31];
        Float_t p_MeVc_proton[31];
        for (int i = 0 ; i < 31; i++) {
            Range_grampercm_proton_LAr[i] = Range_grampercm_proton[i] / LAr_density; // range [cm] = range [gr/cm^2] / rho [gr/cm^3]
            E = 938.272 + KE_MeV_proton[i];
            p_MeVc_proton[i] = sqrt(E*E - 938.272*938.272);
        }
        
        KE_vs_range_proton = new TGraph(31, Range_grampercm_proton_LAr, KE_MeV_proton);
        KE_vs_range_proton_s3 = new TSpline3("KE_vs_range_proton", KE_vs_range_proton);
        p_vs_range_proton = new TGraph(31, Range_grampercm_proton_LAr, p_MeVc_proton);
        p_vs_range_proton_s3 = new TSpline3("p_vs_range_proton", p_vs_range_proton);

        range_vs_p_proton = new TGraph(31, p_MeVc_proton, Range_grampercm_proton_LAr);
        range_vs_p_proton_s3 = new TSpline3("range_vs_p_proton", range_vs_p_proton);

        
        // range in gr/cm2
        Float_t Range_grampercm_muon[29] = {9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1,
            6.732E1, 1.063E2, 1.725E2, 2.385E2, 4.934E2,
            6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3,
            4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4,
            1.910E4, 3.558E4, 4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5
        };
        
        Float_t KE_MeV_muon[29] = {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000,
            1400, 2000, 3000, 4000, 8000, 10000, 14000, 20000, 30000,
            40000, 80000, 100000, 140000, 200000, 300000, 400000
        };
        
        Float_t Range_grampercm_muon_LAr[29];
        Float_t p_MeVc_muon[29];
        for (int i = 0 ; i < 29; i++) {
            Range_grampercm_muon_LAr[i] = Range_grampercm_muon[i] / LAr_density;
            E = 105.6 + KE_MeV_muon[i];
            p_MeVc_muon[i] = sqrt(E*E - 105.6*105.6);
        }
        
        KE_vs_range_muon = new TGraph(29, Range_grampercm_muon_LAr, KE_MeV_muon);
        KE_vs_range_muon_s3 = new TSpline3("KE_vs_range_muon", KE_vs_range_muon);
        p_vs_range_muon = new TGraph(29, Range_grampercm_muon_LAr, p_MeVc_muon);
        p_vs_range_muon_s3 = new TSpline3("p_vs_range_muon", p_vs_range_muon);
        
        range_vs_p_muon = new TGraph(29, p_MeVc_muon, Range_grampercm_muon_LAr);
        range_vs_p_muon_s3 = new TSpline3("range_vs_p_muon", range_vs_p_muon);
    }
    TGraph      *KE_vs_range_proton,    *p_vs_range_proton;
    TSpline3    *KE_vs_range_proton_s3, *p_vs_range_proton_s3;
    TGraph      *KE_vs_range_muon,      *p_vs_range_muon;
    TSpline3    *KE_vs_range_muon_s3,   *p_vs_range_muon_s3;
    TGraph      *range_vs_p_muon,       *range_vs_p_proton;
    TSpline3    *range_vs_p_muon_s3,    *range_vs_p_proton_s3;
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    double Get_protonMomentumFromRange( double range ){      // return the result in MeV
        return p_vs_range_proton_s3 -> Eval(range);
    }
    double Get_muonMomentumFromRange( double range ){        // return the result in MeV
        return p_vs_range_muon_s3 -> Eval(range);
    }
    double Get_protonRangeFromMomentum( double momentum ){   // input the result in MeV, return in cm
        return range_vs_p_proton_s3 -> Eval(momentum);
    }
    double Get_muonRangeFromMomentum( double momentum ){     // input the result in MeV, return in cm
        return range_vs_p_muon_s3 -> Eval(momentum);
    }

    
    
    
private:
    
    TString Path="", RootFileName="", RootTreeName="", OutputCSVname="";
    
    
    TFile * InputGenieFile=nullptr;
    TTree * GenieTree=nullptr;
    ofstream csv_file;

    
    bool    cc, qel, res, dis, coh; // Is it a CC /  quasi-elastic / Resonance / DIS / Coherent scattering?
    bool    mec;
    // my variables for CC 1p 0pi
    bool            IsCC_Np_200MeVc=false,IsCC_1p_200MeVc=false,IsCC_1p_200MeVc_0pi=false;
    
    int     i_event=0;
    int     counter=0;
    int     ni; // Number of 'primary' particles in hadronic system.
    int     neu; //Neutrino PDG code.
    int     nf; // Number of final state particles in hadronic system.
    int     i_gen,  i_rec;

    Int_t   pdgi[NMAX]; //PDG code of kth 'primary' particle in hadronic system (before FSI)
    Int_t   pdgf[NMAX]; //PDG code of kth final state particle in hadronic system
    
    double  Q2, W, x;
    double  pxn, pyn, pzn; // Initial state hit nucleon px (in GeV).
    double  pxv, pyv, pzv; // Incoming neutrino px/y/z (in GeV).
    double  pxl, pyl, pzl; // Final state primary lepton px/py/pz (in GeV)
    double  vtxx, vtxy, vtxz; // vertex position
    double  pxf[NMAX],pyf[NMAX],pzf[NMAX]; //Px/y/z of kth final state particle in hadronic system (in GeV).
    double  pxi[NMAX],pyi[NMAX],pzi[NMAX];
    
    // a weight to the event based on MicroBooNE acceptance
    double  uBacc_truth_muon=0      , uBacc_truth_proton=0  , uBacc_truth_Q2=0;
    double  uBacc_reco_muon=0       , uBacc_reco_proton=0   , uBacc_reco_Q2=0;
    double  side_x_right = 257   ,side_x_left = 0;
    double  side_y_up = 116.5    ,side_y_dw = -116.5;
    double  side_z_downstream = 0,side_z_upstream = 1037;
    double  dmu_cm, dp_cm;      // the stopping range of the muon and proton for their simulated "truth" momentum
    double  reco_l_mu, reco_l_p;   // the length of the muon and proton "tracks"
    double  reco_Ev, reco_Q2, reco_delta_phi, reco_Pt;
    double  gen_Q2_gen_rec, rec_Q2_gen_rec;

    
    TVector3        vertex_position;
    TVector3        muonTrajectory_dir;                 // direction of the muon trajectory
    TVector3        muonTrajectory_end, muonTrack_end;  // the "truth" trajectory vs. the track which might be "broken" by a virtual detector box
    TVector3        protonTrajectory_dir;
    TVector3        protonTrajectory_end, protonTrack_end;
    
    TLorentzVector  pmomentum; //tmp...
    TLorentzVector  proton, muon, nu, q, Pmiss;
    TLorentzVector  proton_before_FSI;   // before FSI
    TLorentzVector  reco_Pp, reco_Pmu, reco_q, reco_Pnu;

    std::vector<TLorentzVector> protons, neutrons;
    
    // event weight for detector-effect emulation
    std::vector<double>              Pmu_xbins, Pmu_theta_ybins;
    std::vector<std::vector<double>> Pmu_theta_acceptance, Pmu_theta_acc_err;
    std::vector<double>              Pp_xbins, Pp_theta_ybins;
    std::vector<std::vector<double>> Pp_theta_acceptance, Pp_theta_acc_err;
    std::vector<double>              Q2_bins;
    std::vector<double>              Q2_acceptance, Q2_acc_err;
    std::vector<double>              Q2_gen_rec_bins;
    std::vector<std::vector<double>> Q2_gen_rec_map;
//    std::vector<std::vector<double>> h_Q2_rec;
    
    std::default_random_engine generator;
    std::vector<std::discrete_distribution<double>> Q2_rec_distributions;
    
    TRandom3  rand;
};

#endif
/** @} */ // end of doxygen group


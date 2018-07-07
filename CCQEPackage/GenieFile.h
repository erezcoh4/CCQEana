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

#include "TFile.h"
#include "TTree.h"
#include "Rtypes.h"
#include "TRandom3.h"
#include "GENIEinteraction.h"

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
    
    // GETters
    int                         GetNevents () const {return (int)GenieTree->GetEntries();};
    bool                  GetCC_Np_200MeVc () const {return IsCC_Np_200MeVc;};
    bool                  GetCC_1p_200MeVc () const {return IsCC_1p_200MeVc;};
    bool              GetCC_1p_200MeVc_0pi () const {return IsCC_1p_200MeVc_0pi;};
    
    // funcionality
    bool                         ReadEvent (int fi_event);
    bool                         HeaderCSV ();
    bool                       StreamToCSV ();
    bool                             Print ();
    bool                       SetTopology ();
    bool                        Initialize ();
    bool                            EndJob ();
    
    std::vector<double>              Read1dArrayFromFile (TString filename);
    std::vector<std::vector<double>> Read2dArrayFromFile (TString filename);
    void                             SetMicroBooNEWeight ();
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
    

    
private:
    
    TString Path="", RootFileName="", RootTreeName="", OutputCSVname="";
    
    
    TFile * InputGenieFile=nullptr;
    TTree * GenieTree=nullptr;
    ofstream csv_file;

    
    int     i_event=0;
    int     counter=0;
    
    // GENIE variables
    int     neu; //Neutrino PDG code.
    bool    cc, qel, res, dis, coh; // Is it a CC /  quasi-elastic / Resonance / DIS / Coherent scattering?
    bool    mec;
    double  Q2, W, x;
    double  pxn, pyn, pzn; // Initial state hit nucleon px (in GeV).
    double  pxv, pyv, pzv; // Incoming neutrino px/y/z (in GeV).
    double  pxl, pyl, pzl; // Final state primary lepton px/py/pz (in GeV)
    
    int     nf; // Number of final state particles in hadronic system.
    Int_t   pdgf[NMAX]; //PDG code of kth final state particle in hadronic system
    double  pxf[NMAX],pyf[NMAX],pzf[NMAX]; //Px/y/z of kth final state particle in hadronic system (in GeV).
    
    int     ni; // Number of 'primary' particles in hadronic system.
    Int_t   pdgi[NMAX]; //PDG code of kth 'primary' particle in hadronic system (before FSI)
    double  pxi[NMAX],pyi[NMAX],pzi[NMAX];
    // a weight to the event based on MicroBooNE acceptance
    double  MicroBooNEWeight_Pmu_theta=0, MicroBooNEWeight_Pp_theta=0,MicroBooNEWeight_Q2=0;
    
    // my variables for CC 1p 0pi
    bool            CC1p0pi=false;
    bool            IsCC_Np_200MeVc=false,IsCC_1p_200MeVc=false,IsCC_1p_200MeVc_0pi=false;
    TLorentzVector  pmomentum; //tmp...
    TLorentzVector  proton, muon, nu, q, Pmiss;
    TLorentzVector  proton_before_FSI;   // before FSI
    std::vector<TLorentzVector> protons, neutrons;
    
    // event weight for detector-effect emulation
    std::vector<double>              Pmu_xbins, Pmu_theta_ybins;
    std::vector<std::vector<double>> Pmu_theta_acceptance, Pmu_theta_acc_err;
    std::vector<double>              Pp_xbins, Pp_theta_ybins;
    std::vector<std::vector<double>> Pp_theta_acceptance, Pp_theta_acc_err;
    std::vector<double>              Q2_bins;
    std::vector<double>              Q2_acceptance, Q2_acc_err;
    
    TRandom3  rand;
};

#endif
/** @} */ // end of doxygen group


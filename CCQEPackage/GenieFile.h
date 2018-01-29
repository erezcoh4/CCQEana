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

#include "TFile.h"
#include "TTree.h"
#include "Rtypes.h"
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
              ,int fdebug );
    

    
    // SETters
    bool    SetInTree ();
    
    
    // GETters
    int           GetNevents () const {return (int)GenieTree->GetEntries();};
    bool          GetCC1p0pi () const {return CC1p0pi;};
    
    // funcionality
    bool           ReadEvent (int fi_event);
    bool           HeaderCSV ();
    bool         StreamToCSV ();
    bool               Print ();
    bool         SetTopology ();
    bool          Initialize ();
    bool              EndJob ();
    
    float     LightConeAlpha ( TLorentzVector v ){ return (v.E() - v.Pz())/(mAr40/40.); };
    
    
private:
    
    TString Path="", RootFileName="", RootTreeName="", OutputCSVname="";
    
    
    TFile * InputGenieFile=nullptr;
    TTree * GenieTree=nullptr;
    ofstream csv_file;

    
    int     i_event=0;
    int     counter=0;
    int     debug=0;
    
    // GENIE variables
    int     neu; //Neutrino PDG code.
    bool    cc, qel, res, dis, coh; // Is it a CC /  quasi-elastic / Resonance / DIS / Coherent scattering?
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
    
    
    // my variables for CC 1p 0pi
    bool            CC1p0pi=false;
    TLorentzVector  proton, muon, nu, q, Pmiss;
    TLorentzVector  proton_before_FSI;   // before FSI
};

#endif
/** @} */ // end of doxygen group


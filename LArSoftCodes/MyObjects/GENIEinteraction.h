/**
 * \file GENIEinteraction.h
 *
 * \ingroup GBDTprotonsPackage
 *
 * \brief Class def header for a class GENIEinteraction
 *
 * @author erezcohen
 */

/** \addtogroup GBDTprotonsPackage
 
 @{*/
#ifndef GENIEinteraction_H
#define GENIEinteraction_H

#include "PandoraNuTrack.h"
#include "TLorentzVector.h"
#include "TVector3.h"

/**
 \class GENIEinteraction
 User defined class GENIEinteraction ... these comments are used to generate
 doxygen documentation!
 */
using namespace std;






class GENIEinteraction{
    
public:
    
    /// Default constructor
    GENIEinteraction () = default;
    GENIEinteraction (Int_t frun, Int_t fsubrun, Int_t fevent, Int_t fmcevent_id);
    
    
    
    
    
    
    // running

    vector<size_t> sort_by_momentum_magnitude(const vector<TVector3> &v);
    bool            SortNucleons ();
    bool       ComputeKinematics ();
    bool        ComputePmissPrec ();
    bool      FindCC1p200MeVc0pi ();
    bool CheckContainement (float max_y = 120,
                            float min_z = 0, float max_z = 1050,
                            float min_x = 0, float max_x = 260){
        if( ( vertex_position.x() < min_x )    | ( vertex_position.x() > max_x ) )    IsVertexContained=false;
        if( ( vertex_position.y() < -max_y )   | ( vertex_position.y() > max_y ) )    IsVertexContained=false;
        if( ( vertex_position.z() < min_z )    | ( vertex_position.z() > max_z ) )    IsVertexContained=false;
        IsVertexContained=true;
        return IsVertexContained;
    }
    void                   Print (bool DoPrintTracks=false) const;
    void                AddTrack (PandoraNuTrack ftrack);
   
   
    
    
    
    // SETters
    void                   SetRSE (int frun , int fsubrun , int fevent)   {run=frun; subrun=fsubrun;event=fevent;};
    void                  SetCCNC (int fccnc)                             {ccnc = fccnc;};
    void                  SetMode (int fmode)                             {mode = fmode;};
    void       SetVertexContained (bool fcontained)                       {IsVertexContained = fcontained;};
    void            SetNprimaries (Int_t n)                               {Nprimaries = n;};
    void                 SetNuPDG (Int_t pdg)                             {nuPDG = pdg;};
    void            SetNuMomentum (TLorentzVector fnu)                    {nu = fnu;};
    void        SetLeptonMomentum (TLorentzVector fmuon)                  {muon = fmuon;}; // if the neutrino is v(e) this will be the e!
    void      SetMomentumTransfer (); // should only be called after SetNuMomentum() and SetLeptonMomentum()
    void        SetVertexPosition (TVector3 fpos);
    void            SetKinematics (Float_t fQ2, Float_t fXb, Float_t fW, Float_t fy, Int_t fccnc, Int_t fmode);
    bool               AddPrimary ( Int_t fpdg                  // pdg code
                                   ,TLorentzVector fmomentum    // 4-momentum
                                   ,Int_t fstatus_code          // status code
                                   ,Int_t fmother               // mother
                                   ,std::string fprocess        // process
    );

    
    
    // GETters
    bool                        GetVertexContained () const {return IsVertexContained;};
    bool                              AskIfCC1p0pi () const {return IsCC_1p_200MeVc_0pi;};
    
    Int_t                             GetMCeventID () const {return mcevent_id;};
    Int_t                                  GetCCNC () const {return ccnc;};
    Int_t                                  GetMode () const {return mode;};
    Int_t                            GetNprimaries () const {return Nprimaries;};
    Int_t                             GetPrimaries () const {return Nprimaries;};
    Int_t                              GetNprotons () const {return protons.size();};

    Float_t                                  GetPt () const {return (protons.size()) ? (muon+protons.at(0)).Pt() : -1;};
    Float_t                           Get_theta_pq () const {return theta_pq;};
    Float_t                                  GetEv () const {return nu.E();};
    Float_t                                  GetQ2 () const {return Q2;};
    Float_t                                  GetXb () const {return Xb;};
    Float_t                                   GetY () const {return y;};
    Float_t                                  GetW2 () const {return W*W;};

    TVector3                     GetVertexPosition () const {return vertex_position;};
    TLorentzVector               GetLeptonMomentum () const {return muon;}; // if the neutrino is v(e) this will be the e!
    TLorentzVector             GetMomentumTransfer () const {return q;}; // if the neutrino is v(e) this will be the e!
    TLorentzVector                          GetPmu () const {return GetLeptonMomentum();};
    TLorentzVector                           GetPp () const {return ((protons.size()>0) ? protons.at(0) : TLorentzVector());};
    
    std::vector<TLorentzVector>  GetProtonsMomenta () const {return protons;}; // when taking this, check if the vector is not empty...
    
    
    

private:
    
    // booleans on the genie interaction
    bool                    IsVertexContained=false;
    bool                    IsCC_1p_200MeVc_0pi=false; // an interaction with at least 1 muon and 1 proton > 200 MeV/c and no pions
    bool                    Is_mu_TrackReconstructed=false, Is_p_TrackReconstructed=false, IsVertexReconstructed=false;

    // Int_t
    Int_t                   nuPDG=-9999;
    Int_t                   run=-9999, subrun=-9999, event=-9999, ccnc=-9999, mode=-9999, mcevent_id=-9999;
    Int_t                   Nprimaries=0;
    Int_t                   Nnu=0, Nnu_e=0, Nnu_mu=0;
    Int_t                   Np=0, Nn=0, Npi=0;
    Int_t                   Nmu=0, Nel=0, Ngamma=0;
    Int_t                   Nmu_minus=0, Nmu_plus=0;
    Int_t                   Npi_minus=0, Npi_plus=0, Npi_0=0;
    Int_t                   Ne_plus=0, Ne_minus=0;
    Int_t                   pos_wire_u=-9999, pos_wire_v=-9999, pos_wire_y=-9999;
    Int_t                   pos_time_u=-9999, pos_time_v=-9999, pos_time_y=-9999;
    
    // Flaot_t
    Float_t                 Xb=-9999 , Q2=-9999, W=-9999, y=-9999;
    Float_t                 theta_pq=-9999, p_over_q=-9999, Mmiss=-9999;
    Float_t                 reco_mu_p_distance=-9999;
    
    // TVector3
    TVector3                vertex_position=TVector3();
    
    // TLorentzVector
    TLorentzVector          momentum=TLorentzVector(), Plead=TLorentzVector() ;
    TLorentzVector          nu=TLorentzVector(-1,-1,-1,-1), muon=TLorentzVector()  , q=TLorentzVector() ;
    TLorentzVector          n_miss=TLorentzVector()  , Pcm=TLorentzVector()   , Prec=TLorentzVector();
    
    // std::vector-s
    std::vector<Int_t>      pdg;           //particle type (pdg) of the GENIE particle
    std::vector<Int_t>      trackID;       //trackID of the GENIE particle (different from the GEANT-assigned track ID)
    std::vector<Int_t>      ND;            //number of daughters of the GENIE particle
    std::vector<Int_t>      mother;        //mother trackID of the GENIE particle
    std::vector<Int_t>      status_code;   //particle status code of the GENIE particle
    
    std::vector<Float_t>    Eng;           //Energy of the GENIE particle in GeV
    std::vector<Float_t>    Px;            //Px of the GENIE particle in GeV
    std::vector<Float_t>    Py;            //Py of the GENIE particle in GeV
    std::vector<Float_t>    Pz;            //Pz of the GENIE particle in GeV
    std::vector<Float_t>    P;             //Magnitude of the momentum vector of the GENIE particle in GeV
    std::vector<Float_t>    mass;          //mass of the GENIE particle in GeV

    std::vector<TVector3>       p3vect  , n3vect;
    std::vector<TLorentzVector> protons , neutrons;

    
    std::vector<std::string>    process;   
    
    // PandoraNuTrack
    PandoraNuTrack              muonTrack=PandoraNuTrack(), protonTrack=PandoraNuTrack();
    std::vector<PandoraNuTrack> tracks; // pandoraNu tracks that are associated with the genie interacion

};
#endif
/** @} */ // end of doxygen group


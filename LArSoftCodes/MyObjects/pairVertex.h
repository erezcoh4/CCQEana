/**
 * \file pairVertex.h
 *
 * \ingroup CCQEPackage
 * 
 * \brief Class def header for a class pairVertex
 *
 * @author erezcohen
 */

/** \addtogroup CCQEPackage

    @{*/
#ifndef PAIRVERTEX_H
#define PAIRVERTEX_H

#include <iostream>
#include "PandoraNuTrack.h"
#include "hit.h"
#include "box.h"
#include "flash.h"
#include "GENIEinteraction.h"




/**
   \class pairVertex
   User defined class pairVertex ... these comments are used to generate
   doxygen documentation!
 */
class pairVertex {

public:

    
    pairVertex() = default;
    pairVertex (Int_t frun,Int_t fsubrun,Int_t fevent,Int_t fID);

    
   
    
    
    
    // functionallity
    bool                           MatchGenieInteraction ( std::vector<GENIEinteraction> , PandoraNuTrack );
    vector<size_t>                                sort_l (const std::vector<PandoraNuTrack> &v);
    bool                              SortTracksByLength ();
    vector<size_t>                             sort_pida (const std::vector<PandoraNuTrack> &v);
    bool                                SortTracksByPIDA ();
    vector<size_t>                       sort_chi2proton (const std::vector<PandoraNuTrack> &v);
    bool                          SortTracksByChi2Proton ();
    
    void                             FixTracksDirections ();
    void                                        AddTrack (PandoraNuTrack ftrack);
    void                                      AddTrackID (Int_t ftrack_id)               {track_id.push_back(ftrack_id);};
    
    void                                          Print (bool DoPrintTracks=false, bool DoPrintFull=true, bool DoPrintGENIE=false) const;
    bool                                  IncludesTrack (Int_t ftrack_id);
    bool                                RemoveFarTracks (float max_mu_p_distance );
    
    std::vector<PandoraNuTrack>  CloseSemiContainedTracks ( std::vector<PandoraNuTrack> AllTracksInTheEvent , float fmax_distance=100 );
    std::vector<PandoraNuTrack>     RemoveTrackFromVector ( std::vector<PandoraNuTrack> AllTracks , PandoraNuTrack TrackToBeRemoved );
    std::vector<PandoraNuTrack>    RemoveTracksFromVector ( std::vector<PandoraNuTrack> AllTracks , std::vector<PandoraNuTrack> TracksToBeRemoved );
    
    void                                ReconstructBeam ();
    void                          ReconstructKinematics ();
    
    void                          AssociateHitsToTracks (std::vector<hit> hits);
    void                        BuildVertexIDFromTracks ();
    
    // @brief check if the recostructed vertex is inside the TPC
    bool                                  IsVertexInTPC (float max_FV_y = 116.5,
                                                         float min_FV_z = 0, float max_FV_z = 1037,
                                                         float min_FV_x = 0, float max_FV_x = 257)  const{
        if( ( position.x() < min_FV_x )    | ( position.x() > max_FV_x ) )    return false;
        if( ( position.y() < -max_FV_y )   | ( position.y() > max_FV_y ) )    return false;
        if( ( position.z() < min_FV_z )    | ( position.z() > max_FV_z ) )    return false;
        return true;
    };
    void                             CheckIsVertexInTPC () { IsRecoVertexInTPC = IsVertexInTPC(); };
    
    
    
    
    
    
    // SETters
    void                SetVertexID (Int_t fvertex_id)                              {vertex_id = fvertex_id;};
    void                     SetRSE (Int_t r, Int_t s, Int_t e)                     {run=r; subrun=s; event=e;};
    void                SetPosition (TVector3 fposition)                            {position = fposition;};
    void                 SetAs1mu1p ();
    void                  SetAsCC1p ();
    void               SetAsCC1p0pi ();
    void              SetAsNon1mu1p ();
    void                SetAsCosmic ();
    void    SetTrueMuonProtonTracks (PandoraNuTrack , PandoraNuTrack );
    void          SetTrueMuonProton (PandoraNuTrack , PandoraNuTrack );
    void         SetTracksRelations ();
    bool         SetIsReconstructed (float fmax_mu_p_distance = 11 );
    void               SetGENIEinfo (GENIEinteraction fgenie_interaction)           {genie_interaction = fgenie_interaction; };
    void            SetClosestGENIE (GENIEinteraction fgenie_interaction)           {closest_genie_interaction = fgenie_interaction; };
    void       SetReconstructedInfo ();
    void            AssignMuonTrack (PandoraNuTrack ftrack)                         {Track_muCandidate = ftrack; };
    void          AssignProtonTrack (PandoraNuTrack ftrack)                         {Track_pCandidate = ftrack; };
    void    SetReconstructedMomenta (float PmuFromRange = 0, float PpFromRange = 0 );
    void   SetReconstructedFeatures (float PmuFromRange = 0, float PpFromRange = 0 );
    void           SetMCSMuMomentum (float fMCSMuMomentum = 0);
    void         SetPlaneProjection (int plane , float _wire , float _time )        {vertex_wire[plane]=_wire; vertex_time[plane]=_time;};
    void            SetClosestFlash (flash _flash)                                  {ClosestFlash = _flash;};
    void SetIsVertexReconstructable (bool fIsRec)                                   {IsVertexReconstructable = fIsRec;};
    void   SetIsVertexReconstructed (bool fIsRec)                                   {IsVertexReconstructed = fIsRec;};
    void       SetIsVertexContained (bool fIsCon)                                   {IsVertexContained = fIsCon;};
    void      SetIsBrokenTrajectory (bool fIsBroken)                                {IsBrokenTrajectory = fIsBroken;};
    void         SetClosestDistance (float fdistances_ij)                           {reco_mu_p_distance=fdistances_ij;};

    void            SetMatchedFlash ( flash _flash, float _fscore )                 {MatchedFlash = _flash; MatchedFlashScore = _fscore;};
    
    
    
    // GETters
    TString         GetTruthTopologyString () const {return TruthTopologyString;}
    
    bool                        GetIs1mu1p () const {return Is1mu1p;};
    bool           GetIsGENIECC_1p_200MeVc () const {return IsGENIECC_1p_200MeVc;};
    bool       GetIsGENIECC_1p_200MeVc_0pi () const {return IsGENIECC_1p_200MeVc_0pi;};
    bool                     GetIsNon1mu1p () const {return IsNon1mu1p;};
    bool                       GetIsCosmic () const {return IsCosmic;};
    bool              GetIsVertexContained () const {return IsVertexContained;};
    bool       GetIs_mu_TrackReconstructed () const {return Is_mu_TrackReconstructed;};
    bool        GetIs_p_TrackReconstructed () const {return Is_p_TrackReconstructed;};
    bool          GetIsVertexReconstructed () const {return IsVertexReconstructed;};
    bool        GetIsVertexReconstructable () const {return IsVertexReconstructable;};
    bool             GetIsBrokenTrajectory () const {return IsBrokenTrajectory;};
    bool              GetIsRecoVertexInTPC () const {return IsRecoVertexInTPC;};
    
    Int_t                           GetRun () const {return run;};
    Int_t                        GetSubrun () const {return subrun;};
    Int_t                         GetEvent () const {return event;};
    Int_t                      GetVertexID () const {return vertex_id;};
    Int_t                       GetNtracks () const {return (Int_t)tracks.size();};
 
    Int_t                    GetVertexWire (int plane) const {return vertex_wire[plane];};
    Int_t                    GetVertexTime (int plane) const {return vertex_time[plane];};
    
    float           GetAngleBetween2tracks () const; // return the angle between the muon and proton candidates, in degrees (!)
    float               GetClosestDistance () const {return reco_mu_p_distance;};
    float                        GetRecoEv () const {return reco_Pnu.E();};
    float                    GetRecoEv_mcs () const {return reco_Pnu_mcs.E();};
    float                        GetReco_q () const {return reco_q.P();};
    float                    GetReco_q_mcs () const {return reco_q_mcs.P();};
    float                    GetReco_omega () const {return reco_omega;};
    float                        GetRecoQ2 () const {return reco_Q2;};
    float                    GetRecoQ2_mcs () const {return reco_Q2_mcs;};
    float                        GetRecoXb () const {return reco_Xb;};
    float                         GetRecoY () const {return reco_y;};
    float                        GetRecoW2 () const {return reco_W2;};
    float                        GetRecoPt () const {return (reco_Pmu + reco_Pp).Pt();};
    float                    GetRecoPt_mcs () const {return (reco_Pmu_mcs + reco_Pp).Pt();};
    float                 GetReco_theta_pq () const {return reco_theta_pq;};
    float                 GetTruthDeltaPhi () const;
    float               GetDistanceToGENIE () const {return (genie_interaction.GetVertexPosition()-position).Mag();};
    float        GetDistanceToClosestGENIE () const {return (closest_genie_interaction.GetVertexPosition()-position).Mag();};

    
    float                       GetRecoEmu () const {return reco_Pmu.E();};
    float                        GetRecoEp () const {return reco_Pp.E();};
    float                        GetRecoTp () const {return reco_Pp.E()-reco_Pp.M();};

    // get the ratio of tracks-charge deposited to total-charge deposited
    // in a box of N(wires) x N(time-ticks) around the vertex in plane i=0,1,2
    // input: plane, N(wires) & N(time-ticks) for the box, hits in event
    // see docdb-10958
    float               GetRdQaroundVertex (int plane, int Nwires, int Nticks, std::vector<hit> hits) const ;
    float                   GetChargeInBox (int plane, std::vector<hit> hits, box VertexBox) const;
    
    // get the ratio of tracks-charge deposited to total-charge deposited
    // in a Sphere of r [cm] around the vertex in plane i=0,1,2
    // input: plane, radius in [cm]
    float         GetRdQSphereAroundVertex (int plane, float r, std::vector<hit> hits ) const ;
    float                GetChargeInSphere (int plane, std::vector<hit> hits, float r) const;
    float             GetMatchedFlashScore () const {return MatchedFlashScore;};
    
    Float_t                  GetYDis2Flash (flash) const;
    Float_t                  GetZDis2Flash (flash) const;
    Float_t                   GetDis2Flash (flash) const;
    Float_t            GetDis2ClosestFlash () const {return GetDis2Flash( ClosestFlash );} ;
    Float_t            GetDis2MatchedFlash () const {return GetDis2Flash( MatchedFlash );};
    Float_t           GetYDis2MatchedFlash () const {return GetYDis2Flash( MatchedFlash );};
    Float_t           GetZDis2MatchedFlash () const {return GetZDis2Flash( MatchedFlash );};
    


    std::vector<float>    Get_delta_phi_ij () const {return delta_phi_ij;};
    std::vector<float>    Get_distances_ij () const {return distances_ij;};
    std::vector<float>  Get_delta_theta_ij () const {return delta_theta_ij;};

    
    
    TVector3                   GetPosition () const {return position;};
    
    TLorentzVector              GetRecoPnu () const {return reco_Pnu;};
    TLorentzVector              GetRecoPmu () const {return reco_Pmu;};
    TLorentzVector          GetRecoPnu_mcs () const {return reco_Pnu_mcs;};
    TLorentzVector          GetRecoPmu_mcs () const {return reco_Pmu_mcs;};
    TLorentzVector               GetRecoPp () const {return reco_Pp;};
    TLorentzVector            GetRecoPmiss () const {return (reco_Pp - reco_q);};
    
    
    PandoraNuTrack        GetShortestTrack () const {return ShortestTrack;};
    PandoraNuTrack         GetLongestTrack () const {return LongestTrack;};
    PandoraNuTrack       GetSmallPIDaTrack () const {return SmallPIDaTrack;};
    PandoraNuTrack       GetLargePIDaTrack () const {return LargePIDaTrack;};
    PandoraNuTrack GetSmallChi2ProtonTrack () const {return SmallChi2ProtonTrack;};
    PandoraNuTrack GetLargeChi2ProtonTrack () const {return LargeChi2ProtonTrack;};
    PandoraNuTrack    GetTrack_muCandidate () const {return Track_muCandidate;};
    PandoraNuTrack     GetTrack_pCandidate () const {return Track_pCandidate;};
    
    std::vector<PandoraNuTrack>  GetTracks () const {return tracks;};
    
    GENIEinteraction          GetGENIEinfo () const {return genie_interaction;};
    GENIEinteraction       GetClosestGENIE () const {return closest_genie_interaction;};
    flash                  GetClosestFlash () const {return ClosestFlash;};
    flash                  GetMatchedFlash () const {return MatchedFlash;};

    std::vector<hit>           GetMuonHits (int plane) const {return hits_muon[plane];};
    std::vector<hit>         GetProtonHits (int plane) const {return hits_proton[plane];};
    
    
    // operators
    inline bool operator==(const pairVertex & v) {
        return  (vertex_id == v.GetVertexID());
    }
    inline bool operator!=(const pairVertex & v) {
        return  (vertex_id != v.GetVertexID());
    }
    
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    // debug
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    Int_t debug=0;
    void Debug (Int_t verobosity_level, std::string text){
        if ( debug > verobosity_level ) cout << text << endl;
    }
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
   
private:
    
    
    // -------------------------------------------------------
    
    // variables
    TString             TruthTopologyString="unknown truth topology";
    
    bool                Is1mu1p=false,   IsNon1mu1p=false,   IsCosmic=false;
    bool                IsGENIECC_1p_200MeVc=false; // CC event with a single proton above 200 MeV/c and no charged pions above 70 MeV/c
    bool                IsGENIECC_1p_200MeVc_0pi=false; // CC with a single proton above 200 MeV/c and no pions or photons or electrons
    bool                IsVertexContained=false, Is_mu_TrackReconstructed=false, Is_p_TrackReconstructed=false;
    bool                IsVertexReconstructed=false, IsVertexReconstructable=false;
    bool                IsBrokenTrajectory=false, IsRecoVertexInTPC=false;

    Int_t               run=-1 , subrun=-1 , event=-1, vertex_id=-1;
    Int_t               Ntracks=-1;
    
    // location in each plane
    float               vertex_wire[3]={0,0,0} , vertex_time[3]={0,0,0}; // these are floating point numbers since they are projections
    
    // reconstructed features
    // calorimentric reconstruction:
    // Ev = Eµ + Tp + Sn + T(A-1)
    float               reco_mu_p_distance=-1;
    float               reco_BeamE=-1,   reco_theta_pq=-1, reco_Pp_3momentum=-1, reco_Pmu_3momentum=-1;
    float               reco_p_over_q=-1, reco_Q2=-1, reco_Q2_mcs=-1;
    float               reco_omega=-1;
    float               reco_Xb=-1, reco_y=-1, reco_W2=-1, reco_s=-1;
    float               reco_alpha_p=-1 , reco_alpha_q=-1 , reco_alpha_mu=-1, reco_alpha_miss=-1;
    // --- - - --- -- - -- -- -- -- --- -- - --- - -- - -- -- -- --- - -- - --- - - -- - -- -
    
    float               truth_alpha_q, truth_alpha_p, truth_alpha_mu, truth_alpha_miss;
    float               MatchedFlashScore=-1;   // score of flash-matchig by Marco
    
    TVector3            position=TVector3();
    TVector3            reco_Pp_3vect=TVector3(), reco_Pmu_3vect=TVector3(), reco_Pmu_3vect_mcs=TVector3();
    
    // Tp + Eµ
    TLorentzVector      reco_Pnu=TLorentzVector(-1,-1,-1,-1),  reco_Pp=TLorentzVector(-1,-1,-1,-1);
    TLorentzVector      reco_Pmu=TLorentzVector(-1,-1,-1,-1),  reco_q=TLorentzVector(-1,-1,-1,-1);
    TLorentzVector      reco_Pnu_mcs=TLorentzVector(-1,-1,-1,-1),   reco_Pmu_mcs=TLorentzVector(-1,-1,-1,-1), reco_q_mcs=TLorentzVector(-1,-1,-1,-1);
    TLorentzVector      reco_n_miss=TLorentzVector(-1,-1,-1,-1);
    
    PandoraNuTrack      muonTrueTrack,  protonTrueTrack;
    PandoraNuTrack      ShortestTrack,  LongestTrack;
    PandoraNuTrack      LargePIDaTrack, SmallPIDaTrack;
    PandoraNuTrack      LargeChi2ProtonTrack, SmallChi2ProtonTrack;
    PandoraNuTrack      Track_muCandidate, Track_pCandidate;
    
    GENIEinteraction    genie_interaction=GENIEinteraction();
    GENIEinteraction    closest_genie_interaction=GENIEinteraction();

    std::vector<Int_t>  track_id;

    std::vector <std::vector<float> >   tracks_distances;
    std::vector <std::vector<float> >   tracks_delta_phi;
    std::vector <std::vector<float> >   tracks_delta_theta;
    std::vector<float>                  tracks_dis_from_vertex, delta_phi_ij,    distances_ij , delta_theta_ij;
    std::vector<PandoraNuTrack>         tracks, tracks_lengthsorted,  tracks_pidasorted, tracks_chi2protonsorted ;
    
    std::vector<hit>    hits_muon[3], hits_proton[3]; // in 3 wire planes
    flash               ClosestFlash=flash();
    flash               MatchedFlash=flash();           // Flash-matchig by Marco
    
};

#endif
/** @} */ // end of doxygen group 


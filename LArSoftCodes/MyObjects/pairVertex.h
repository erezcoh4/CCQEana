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
#include "GENIEinteraction.h"

#define r2d TMath::RadToDeg()
#define d2r TMath::DegToRad()
#define PI TMath::Pi()



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
    
    void                             FixTracksDirections ();
    void                                        AddTrack (PandoraNuTrack ftrack);
    void                                      AddTrackID (Int_t ftrack_id)               {track_id.push_back(ftrack_id);};
    
    void                                          Print (bool DoPrintTracks=false) const;
    bool                                  IncludesTrack (Int_t ftrack_id);
    bool                                RemoveFarTracks (float max_mu_p_distance );
    
    std::vector<PandoraNuTrack>  CloseSemiContainedTracks ( std::vector<PandoraNuTrack> AllTracksInTheEvent , float fmax_distance=100 );
    std::vector<PandoraNuTrack>     RemoveTrackFromVector ( std::vector<PandoraNuTrack> AllTracks , PandoraNuTrack TrackToBeRemoved );
    std::vector<PandoraNuTrack>    RemoveTracksFromVector ( std::vector<PandoraNuTrack> AllTracks , std::vector<PandoraNuTrack> TracksToBeRemoved );
    
    void                                ReconstructBeam ();
    void                          ReconstructKinematics ();
    
    void                          AssociateHitsToTracks (std::vector<hit> hits);
    
    
    // SETters
    void                SetVertexID (Int_t fvertex_id)                              {vertex_id = fvertex_id;};
    void                     SetRSE (Int_t r, Int_t s, Int_t e)                     {run=r; subrun=s; event=e;};
    void                SetPosition (TVector3 fposition)                            {position = fposition;};
    void                 SetAs1mu1p ();
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
    void            AssignMuonTrack (PandoraNuTrack ftrack)                         {AssignedMuonTrack = ftrack; };
    void          AssignProtonTrack (PandoraNuTrack ftrack)                         {AssignedProtonTrack = ftrack; };
    void    SetReconstructedMomenta (float PmuFromRange = 0, float PpFromRange = 0 );
    void   SetReconstructedFeatures (float PmuFromRange = 0, float PpFromRange = 0 );
    void         SetPlaneProjection (int plane , float _wire , float _time )        {vertex_wire[plane]=_wire; vertex_time[plane]=_time;};

    // GETters
    
    bool                        GetIs1mu1p () const {return Is1mu1p;};
    bool       GetIsGENIECC_1p_200MeVc_0pi () const {return IsGENIECC_1p_200MeVc_0pi;};
    bool                     GetIsNon1mu1p () const {return IsNon1mu1p;};
    bool                       GetIsCosmic () const {return IsCosmic;};
    bool              GetIsVertexContained () const {return IsVertexContained;};
    bool       GetIs_mu_TrackReconstructed () const {return Is_mu_TrackReconstructed;};
    bool        GetIs_p_TrackReconstructed () const {return Is_p_TrackReconstructed;};
    bool          GetIsVertexReconstructed () const {return IsVertexReconstructed;};

    
    Int_t                           GetRun () const {return run;};
    Int_t                        GetSubrun () const {return subrun;};
    Int_t                         GetEvent () const {return event;};
    Int_t                      GetVertexID () const {return vertex_id;};
    Int_t                       GetNtracks () const {return (Int_t)tracks.size();};
 
    
    float           GetAngleBetween2tracks () const; // return the angle between the muon and proton candidates, in degrees (!)
    float                        GetRecoEv () const {return reco_Pnu.E();};
    float                        GetRecoQ2 () const {return reco_Q2;};
    float                        GetRecoXb () const {return reco_Xb;};
    float                         GetRecoY () const {return reco_y;};
    float                        GetRecoW2 () const {return reco_W2;};
    float                        GetRecoPt () const {return (IsVertexReconstructed) ? (reco_Pmu + reco_Pp).Pt() : -1;};
    float                 GetReco_theta_pq () const {return reco_theta_pq;};
    float                 GetTruthDeltaPhi () const;
    
    // get the ratio of tracks-charge deposited to total-charge deposited
    // in a box of N(wires) x N(time-ticks) around the vertex in plane i=0,1,2
    // input: plane, N(wires) & N(time-ticks) for the box, hits in event
    float               GetRdQaroundVertex (int plane, int Nwires, int Nticks, std::vector<hit> hits) const ;
    
    std::vector<float>    Get_delta_phi_ij () const {return delta_phi_ij;};
    std::vector<float>    Get_distances_ij () const {return distances_ij;};
    std::vector<float>  Get_delta_theta_ij () const {return delta_theta_ij;};

    
    
    TVector3                   GetPosition () const {return position;};
    
    TLorentzVector              GetRecoPnu () const {return reco_Pnu;};
    TLorentzVector              GetRecoPmu () const {return reco_Pmu;};
    TLorentzVector               GetRecoPp () const {return ((IsVertexReconstructed) ? (reco_Pp) : TLorentzVector() );};
    
    
    PandoraNuTrack        GetShortestTrack () const {return ShortestTrack;};
    PandoraNuTrack         GetLongestTrack () const {return LongestTrack;};
    PandoraNuTrack       GetSmallPIDaTrack () const {return SmallPIDaTrack;};
    PandoraNuTrack       GetLargePIDaTrack () const {return LargePIDaTrack;};
    PandoraNuTrack    GetAssignedMuonTrack () const {return AssignedMuonTrack;};
    PandoraNuTrack  GetAssignedProtonTrack () const {return AssignedProtonTrack;};
    
    std::vector<PandoraNuTrack>  GetTracks () const {return tracks;};
    
    GENIEinteraction          GetGENIEinfo () const {return genie_interaction;};
    GENIEinteraction       GetClosestGENIE () const {return closest_genie_interaction;};

    
    
    // operators
    inline bool operator==(const pairVertex & v) {
        return  (vertex_id == v.GetVertexID());
    }
    inline bool operator!=(const pairVertex & v) {
        return  (vertex_id != v.GetVertexID());
    }
   
private:
    
    
//    Float_t             AllChargeInVertexROI[3], AllChargeInVertexROI_enlarged_20_100[3], AllChargeInVertexROI_enlarged_40_200[3]; // sum of charge of all hits in the vertex-roi per plane
//    Float_t             dQtotROI_20x40_AroundVertex[3], dQassociatedROI_20x40_AroundVertex[3];
//    Float_t             TracksAssociatedCharge[3]; // sum of charge of all hits that are associated with my-tracks
//    Float_t             ratio_associated_hit_charge_to_total[3] , average_ratio_associated_hit_charge_to_total , max_ratio_associated_hit_charge_to_total;
//    Float_t             ratio_associated_hit_charge_to_total_enlarged_20_100[3];
//    Float_t             ratio_associated_hit_charge_to_total_enlarged_40_200[3];
//    Float_t             ratio_dQassociated_dQtot_ROI_20x40_AroundVertex[3];
//    
//    MyTrack             MyTrackMuonTrack[3] , MyTrackProtonTrack[3];
//    MyTrack             MyTrackMuon_u, MyTrackMuon_v, MyTrackMuon_y, MyTrackProton_u, MyTrackProton_v, MyTrackProton_y;
//    std::vector<MyTrack> my_tracks;
//    
//    hit                 ClosestHitToVertex[3];
//    std::vector<hit>    HitsInPlane_u, HitsInPlane_v, HitsInPlane_y;
//    std::vector<hit>    AllHitsInROI[3], AllHitsInROI_u,AllHitsInROI_v,AllHitsInROI_y;
    // -------------------------------------------------------
    
    
    
    
    
    

    
    // variables
    TString             TopologyString="" , TruthTopologyString="unknown truth topology";
    
    bool                Is1mu1p=false,    IsGENIECC_1p_200MeVc_0pi=false,   IsNon1mu1p=false,   IsCosmic=false;
    bool                IsVertexContained=false, Is_mu_TrackReconstructed=false, Is_p_TrackReconstructed=false, IsVertexReconstructed=false;

    Int_t               run=-1 , subrun=-1 , event=-1, vertex_id=-1;
    Int_t               Ntracks=-1;
    
    // location in each plane
    float               vertex_wire[3]={0,0,0} , vertex_time[3]={0,0,0}; // these are floating point numbers since they are projections
    
    // reconstructed features
    // calorimentric reconstruction:
    // Ev = Eµ + Tp + Sn + T(A-1)
    float               reco_mu_p_distance=-1;
    float               reco_BeamE=-1,   reco_theta_pq=-1, reco_Pp_3momentum=-1, reco_Pmu_3momentum=-1;
    float               reco_p_over_q=-1, reco_Q2=-1;
    float               reco_omega=-1;
    float               reco_Xb=-1, reco_y=-1, reco_W2=-1, reco_s=-1;
    float               reco_alpha_p=-1 , reco_alpha_q=-1 , reco_alpha_mu=-1, reco_alpha_miss=-1;
    // --- - - --- -- - -- -- -- -- --- -- - --- - -- - -- -- -- --- - -- - --- - - -- - -- -
    
    float               truth_alpha_q, truth_alpha_p, truth_alpha_mu, truth_alpha_miss;
    
    //    float               dqdx_around_vertex,   dqdx_around_vertex_tracks_associated, dqdx_around_vertex_non_tracks_associated;

    TVector3            position=TVector3();
    TVector3            reco_Pp_3vect=TVector3(), reco_Pmu_3vect=TVector3();
//
//    
//    TLorentzVector      reconstructed_nu, reconstructed_muon, reconstructed_q ;
    
    // Tp + Eµ
    TLorentzVector      reco_Pnu=TLorentzVector(-1,-1,-1,-1),  reco_Pp=TLorentzVector(-1,-1,-1,-1);
    TLorentzVector      reco_Pmu=TLorentzVector(-1,-1,-1,-1),  reco_q=TLorentzVector(-1,-1,-1,-1);
    TLorentzVector      reco_n_miss=TLorentzVector(-1,-1,-1,-1);
    
//    // momentum correction from p(mu)/theta(mu) and p(p) / theta(p) correlations
//    TLorentzVector      reco_Pp_corrected,  reco_Pmu_corrected,  reco_q_corrected;
//    TLorentzVector      reco_Pnu_corrected, reco_n_miss_corrected;
//    // --- - - --- -- - -- -- -- -- --- -- - --- - -- - -- -- -- --- - -- - --- - - -- - -- -
//    
//    
//    // from MCS LLHD
//    TVector3            reco_Pmu_3vect_mcsllhd;
//    TLorentzVector      reco_Pmu_mcsllhd, reco_Pnu_mcsllhd, reco_q_mcsllhd;
//    float               reco_theta_pq_mcsllhd, reco_BeamPz_mcsllhd;
//    float               reco_p_over_q_mcsllhd, reco_Q2_mcsllhd, reco_Q2_from_angles_mcsllhd;
//    // --------------------------------------------------------------------------------------------------------
//    
//    box                 roi[3] , roi_u , roi_v , roi_y, Roi_20x40_AroundVertex[3];
    
    PandoraNuTrack      muonTrueTrack,  protonTrueTrack;
    PandoraNuTrack      ShortestTrack,  LongestTrack;
    PandoraNuTrack      LargePIDaTrack, SmallPIDaTrack;
    PandoraNuTrack      AssignedMuonTrack, AssignedProtonTrack;
    
    GENIEinteraction    genie_interaction=GENIEinteraction();
    GENIEinteraction    closest_genie_interaction=GENIEinteraction();

    std::vector<Int_t>  track_id;

    std::vector <std::vector<float> >   tracks_distances;
    std::vector <std::vector<float> >   tracks_delta_phi;
    std::vector <std::vector<float> >   tracks_delta_theta;
    std::vector<float>                  tracks_dis_from_vertex, delta_phi_ij,    distances_ij , delta_theta_ij;
    std::vector<PandoraNuTrack>         tracks, tracks_lengthsorted,  tracks_pidasorted ;
    
    std::vector<hit>    hits_muon[3], hits_proton[3]; // in 3 wire planes
    
};

#endif
/** @} */ // end of doxygen group 

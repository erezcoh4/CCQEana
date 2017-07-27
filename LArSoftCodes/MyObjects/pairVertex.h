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
    bool                                        Initiate ();
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
    
    std::vector<PandoraNuTrack>   NearUncontainedTracks ( std::vector<PandoraNuTrack> AllTracksInTheEvent , float fmax_distance=100 );
    std::vector<PandoraNuTrack>   RemoveTrackFromVector ( std::vector<PandoraNuTrack> AllTracks , PandoraNuTrack TrackToBeRemoved );
    std::vector<PandoraNuTrack>  RemoveTracksFromVector ( std::vector<PandoraNuTrack> AllTracks , std::vector<PandoraNuTrack> TracksToBeRemoved );
    
    
    
    
    // SETters
    void                SetVertexID (Int_t fvertex_id)              {vertex_id = fvertex_id;};
    void                     SetRSE (Int_t r, Int_t s, Int_t e)     {run=r; subrun=s; event=e;};
    void                SetPosition (TVector3 fposition)            {position = fposition;};
    void                 SetAs1mu1p ();
    void               SetAsCC1p0pi ();
    void              SetAsNon1mu1p ();
    void                SetAsCosmic ();
    void        SetMuonProtonTracks ( PandoraNuTrack , PandoraNuTrack );
    void              SetMuonProton ( PandoraNuTrack , PandoraNuTrack );
    void         SetTracksRelations ();
    bool         SetIsReconstructed ( float fmax_mu_p_distance = 11 );
    void               SetGENIEinfo (GENIEinteraction fgenie_interaction){ genie_interaction = fgenie_interaction; };
    void            SetClosestGENIE (GENIEinteraction fgenie_interaction){ closest_genie_interaction = fgenie_interaction; };
    void       SetReconstructedInfo ();
    void            SetAssignTracks (PandoraNuTrack fAssignedMuonTrack, PandoraNuTrack fAssignedProtonTrack,
                                     float PmuonFromRange = 0, float PprotonFromRange = 0 );
    void    SetReconstructedMomenta ( float PmuonFromRange = 0, float PprotonFromRange = 0 );
    void     SetReconstructedBeamPz ();
    void         SetReconstructed_q ();
    void   SetReconstructedFeatures ( float PmuonFromRange = 0, float PprotonFromRange = 0 );
    

    // GETters
    Int_t                           GetRun () const {return run;};
    Int_t                        GetSubrun () const {return subrun;};
    Int_t                         GetEvent () const {return event;};
    Int_t                      GetVertexID () const {return vertex_id;};
    
    TVector3                   GetPosition () const {return position;};
    
    std::vector<PandoraNuTrack>  GetTracks () const {return tracks;};
    GENIEinteraction          GetGENIEinfo () const {return genie_interaction;};
    GENIEinteraction       GetClosestGENIE () const {return closest_genie_interaction;};

    
    
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
    
    bool                Is1mu1p=false,    IsGENIECC_1p_200MeVc_0pi=false,   Non1mu1p=false,   IsCosmic=false;
    bool                IsVertexContained=false, IsmuTrackReconstructed=false, IspTrackReconstructed=false, IsVertexReconstructed=false;

    Int_t               run=-1 , subrun=-1 , event=-1, vertex_id=-1;
    Int_t               Ntracks=-1;
//    Int_t               reconstructed_Np, reconstructed_Nn, reconstructed_Npi, reconstructed_Nmu, reconstructed_Nel;
//    
//    // location in each plane
//    float               vertex_wire[3] , vertex_time[3];
//    float               delta_phi_LongestShortestTracks;
//    float               reconstructed_Xb, reconstructed_Q2 ;
//    float               reconstructed_theta_pq, reconstructed_p_over_q, reconstructed_Mmiss;
//    
//    // CC1p reconstructed features
//    float               reco_mu_p_distance;
//    float               reco_BeamPz,   reco_theta_pq, reco_Pp_3momentum, reco_Pmu_3momentum;
//    float               reco_p_over_q, reco_Q2, reco_Q2_from_angles;
//    float               reco_omega;
//    float               reco_Xb, reco_y, reco_W2, reco_s;
//    float               reco_Ev_from_angles, reco_Ev_from_angles_Ev_from_mu_p_diff, reco_Ev_from_angles_Ev_from_mu_p_ratio;
//    float               reco_alpha_p , reco_alpha_q , reco_alpha_mu, reco_alpha_miss;
//    float               reco_Q2_from_angles_diff, reco_Q2_from_angles_ratio;
//    float               reco_Ev_with_binding, reco_Ev_with_binding_diff, reco_Ev_with_binding_ratio;
//
//    
//    // momentum correction from p(mu)/theta(mu) and p(p) / theta(p) correlations
//    float               reco_Ev_corrected;
//    float               reco_theta_pq_corrected,   reco_p_over_q_corrected, reco_Q2_corrected;
//    float               reco_omega_corrected;
//    float               reco_Xb_corrected, reco_y_corrected, reco_W2_corrected, reco_s_corrected;
//    float               reco_alpha_q_corrected, reco_alpha_p_corrected, reco_alpha_miss_corrected, reco_alpha_mu_corrected;
//    // --- - - --- -- - -- -- -- -- --- -- - --- - -- - -- -- -- --- - -- - --- - - -- - -- -
//
//    // Tp + Eµ
//    float               reco_BeamE, reco_Ev_fromE;
//    float               reco_theta_pq_fromE,   reco_p_over_q_fromE, reco_Q2_fromE;
//    float               reco_omega_fromE;
//    float               reco_Xb_fromE, reco_y_fromE, reco_W2_fromE, reco_s_fromE;
//    float               reco_alpha_p_fromE, reco_alpha_q_fromE, reco_alpha_miss_fromE;
//    // --- - - --- -- - -- -- -- -- --- -- - --- - -- - -- -- -- --- - -- - --- - - -- - -- -
//
//    
//    float               dqdx_around_vertex,   dqdx_around_vertex_tracks_associated, dqdx_around_vertex_non_tracks_associated;
//    float               truth_alpha_q, truth_alpha_p, truth_alpha_mu, truth_alpha_miss;
//    
    TVector3            position=TVector3()    ;
//    TVector3            reco_Pp_3vect, reco_Pmu_3vect, reco_Pp_3vect_corrected, reco_Pmu_3vect_corrected;
//    
//    
//    TLorentzVector      reconstructed_nu, reconstructed_muon, reconstructed_q ;
//    TLorentzVector      reco_Pnu,  reco_Pp,   reco_Pmu,  reco_q;
//    TLorentzVector      reco_n_miss;    
//    // Tp + Eµ
//    TLorentzVector      reco_Pnu_fromE, reco_q_fromE, reco_n_miss_fromE;
//
//    
//    
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
//    
//
//    PandoraNuTrack      muonTrueTrack,  protonTrueTrack;
    PandoraNuTrack      ShortestTrack,  LongestTrack;
    PandoraNuTrack      LargePIDATrack, SmallPIDATrack;
    PandoraNuTrack      AssignedMuonTrack, AssignedProtonTrack;
    
    GENIEinteraction    genie_interaction=GENIEinteraction();
    GENIEinteraction    closest_genie_interaction=GENIEinteraction();

    std::vector<Int_t>  track_id;

    std::vector <std::vector<float> >   tracks_distances;
    std::vector <std::vector<float> >   tracks_delta_phi;
    std::vector <std::vector<float> >   tracks_delta_theta;
    std::vector<float>                  tracks_dis_from_vertex, delta_phi_ij,    distances_ij , delta_theta_ij;
//
//    
//    std::vector<TLorentzVector> reconstructed_protons;
//    
    std::vector<PandoraNuTrack>         tracks, tracks_lengthsorted,  tracks_pidasorted ;
    
    
    
};

#endif
/** @} */ // end of doxygen group 


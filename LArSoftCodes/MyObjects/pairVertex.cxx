#ifndef PAIRVERTEX_CXX
#define PAIRVERTEX_CXX

#include "pairVertex.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
pairVertex::pairVertex(Int_t frun,Int_t fsubrun,Int_t fevent,Int_t fID):
run(frun),
subrun(fsubrun),
event(fevent),
vertex_id(fID)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::AddTrack (PandoraNuTrack ftrack){
    
    tracks.push_back(ftrack);
    Ntracks=(int)tracks.size();
    AddTrackID( ftrack.GetTrackID() );
    
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::SetAs1mu1p(){
    Is1mu1p                     = true;
    IsGENIECC_1p_200MeVc_0pi    = false;
    IsNon1mu1p                  = false;
    IsCosmic                    = false;
    TruthTopologyString = "true µp pair";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::SetAsCC1p0pi(){
    Is1mu1p                     = true;
    IsGENIECC_1p_200MeVc_0pi    = true;
    IsNon1mu1p                  = false;
    IsCosmic                    = false;
    TruthTopologyString = "true CC 1p 0π";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::SetAsNon1mu1p(){
    Is1mu1p                     = false;
    IsGENIECC_1p_200MeVc_0pi    = false;
    IsNon1mu1p                  = true;
    IsCosmic                    = false;
    TruthTopologyString = "true non-µp pair";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::SetAsCosmic(){
    IsGENIECC_1p_200MeVc_0pi    = false;
    IsNon1mu1p                  = false;
    IsCosmic                    = true;
    TruthTopologyString = "cosmic pair";
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool pairVertex::IncludesTrack (Int_t ftrack_id){
    for (auto track:tracks){
        if (track.GetTrackID() == ftrack_id) return true;
    }
    return false;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool pairVertex::RemoveFarTracks(float max_mu_p_distance ){
    
    // after we fixed the vertex position,
    // narrow down the set of tracks associated to the vertex by looking only
    // at those tracks that are close enough to the vertex
    std::vector<Int_t>          CloseEnoughTracksID;
    std::vector<PandoraNuTrack> CloseEnoughTracks;
    
    for (auto t: tracks) {
        if ( t.DistanceFromPoint( position ) < max_mu_p_distance ){
            CloseEnoughTracks.push_back( t );
            CloseEnoughTracksID.push_back( t.GetTrackID() );
        }
    }
    tracks = CloseEnoughTracks;
    track_id = CloseEnoughTracksID;
    Ntracks = tracks.size();
    return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<size_t> pairVertex::sort_l(const vector<PandoraNuTrack> &v) {
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1].GetLength() > v[i2].GetLength();});
    return idx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool pairVertex::SortTracksByLength(){
    std::vector<PandoraNuTrack> tmp_tracks = tracks;
    for (auto i:sort_l( tmp_tracks )){
        tracks_lengthsorted.push_back( tmp_tracks.at(i) );
    }
    ShortestTrack = tracks_lengthsorted.at( tmp_tracks.size() - 1 );
    LongestTrack = tracks_lengthsorted.at( 0 );
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<size_t> pairVertex::sort_pida(const vector<PandoraNuTrack> &v) {
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1].GetPIDa() > v[i2].GetPIDa();});
    return idx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool pairVertex::SortTracksByPIDA(){
    std::vector<PandoraNuTrack> tmp_tracks = tracks;
    for (auto i:sort_pida( tmp_tracks )){
        tracks_pidasorted.push_back( tmp_tracks.at(i) );
    }
    SmallPIDaTrack = tracks_pidasorted.at( tmp_tracks.size() - 1 );
    LargePIDaTrack = tracks_pidasorted.at( 0 );
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::SetTracksRelations(){
    // July 2017
    // the angles (delta_phi and delta_theta) are calculated in degrees (!)
    for(int i = 0; i < Ntracks ; i++){
        
        tracks_dis_from_vertex.push_back( tracks[i].DistanceFromPoint( position ) );
        
        std::vector<float> distances_track_i;
        std::vector<float> delta_phi_track_i;
        std::vector<float> delta_theta_track_i;
        
        for(int j = 0; j < Ntracks ; j++){
            distances_track_i.push_back( tracks[i].ClosestDistanceToOtherTrack(tracks[j]) );
            delta_phi_track_i.push_back( r2d*fabs(tracks[i].GetPhi() - tracks[j].GetPhi()) );
            delta_theta_track_i.push_back( r2d*fabs(tracks[i].GetTheta() - tracks[j].GetTheta() ) );
            
            
            if (i==0 && j!=i) {
                distances_ij.push_back( distances_track_i.back() );
                delta_phi_ij.push_back( delta_phi_track_i.back() );
                delta_theta_ij.push_back( delta_theta_track_i.back() );
            }
        }
        tracks_distances.push_back( distances_track_i );
        tracks_delta_phi.push_back( delta_phi_track_i );
        tracks_delta_theta.push_back( delta_theta_track_i );
        
        
        distances_track_i.clear();
        delta_phi_track_i.clear();
        delta_theta_track_i.clear();
    }
    
    
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<PandoraNuTrack> pairVertex::RemoveTrackFromVector( vector<PandoraNuTrack> TracksVector , PandoraNuTrack TrackToBeRemoved ){
    std::vector<PandoraNuTrack> tmp_tracks = TracksVector;
    TracksVector.clear();
    for (auto t : tmp_tracks){
        if ( t != TrackToBeRemoved ){
            TracksVector.push_back(t);
        }
    }
    return TracksVector;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<PandoraNuTrack> pairVertex::RemoveTracksFromVector( vector<PandoraNuTrack> TracksVector , vector<PandoraNuTrack> TracksToBeRemoved ){
    // loop over tracks to be removed
    // and remove each of them
    for (auto t:TracksToBeRemoved) TracksVector = RemoveTrackFromVector( TracksVector , t );
    return TracksVector;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<PandoraNuTrack> pairVertex::CloseSemiContainedTracks( vector<PandoraNuTrack> AllTracksInTheEvent , float fMaxDistance ){
    
    vector<PandoraNuTrack> NonVertexTracks = RemoveTracksFromVector( AllTracksInTheEvent , tracks );
    std::vector<PandoraNuTrack> CloseTracks;
    for (auto track:NonVertexTracks) {
        if ( track.DistanceFromPoint(position) < fMaxDistance ){
            CloseTracks.push_back(track);
        }
    }
    return CloseTracks;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::SetTrueMuonProtonTracks( PandoraNuTrack track_a , PandoraNuTrack track_b ){
    protonTrueTrack = track_a;
    muonTrueTrack = track_b;
    Is_mu_TrackReconstructed = Is_p_TrackReconstructed = IsVertexReconstructed = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::SetTrueMuonProton( PandoraNuTrack t1 , PandoraNuTrack t2 ){
    if ( t1.GetMCpdgCode()==2212 && t2.GetMCpdgCode()==13 ){
        SetTrueMuonProtonTracks( t1 , t2 );
    }
    else if ( t1.GetMCpdgCode()==13 && t2.GetMCpdgCode()==2212 ){
        SetTrueMuonProtonTracks( t2 , t1 );
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float pairVertex::GetAngleBetween2tracks() const{
    // July-30 2017
    // return the angle between the two tracks in the vertex, in degrees (!)
    TVector3 t1_dir, t2_dir;
    t1_dir.SetMagThetaPhi ( 1 , AssignedMuonTrack.GetTheta(), AssignedMuonTrack.GetPhi() );
    t2_dir.SetMagThetaPhi ( 1 , AssignedProtonTrack.GetTheta(), AssignedProtonTrack.GetPhi() );
    return r2d*(t1_dir.Angle( t2_dir ));
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::FixTracksDirections(){
    // for CC1p events, we can fix the directions of the track
    // by looking at the reconstructed vertex position
    // and comparing the start/end point of the track to the vertex position
    float start_start_distance  = (AssignedMuonTrack.GetStartPos()  - AssignedProtonTrack.GetStartPos()).Mag();
    float end_start_distance    = (AssignedMuonTrack.GetEndPos()    - AssignedProtonTrack.GetStartPos()).Mag();
    float start_end_distance    = (AssignedMuonTrack.GetStartPos()  - AssignedProtonTrack.GetEndPos()).Mag();
    float end_end_distance      = (AssignedMuonTrack.GetEndPos()    - AssignedProtonTrack.GetEndPos()).Mag();
    float min_distance = std::min({start_start_distance, end_start_distance, start_end_distance, end_end_distance});
    
    // first fix the position of the vertex
    if (min_distance == start_start_distance){
        position = 0.5*(AssignedMuonTrack.GetStartPos() + AssignedProtonTrack.GetStartPos());
    }
    else if (min_distance == end_start_distance){
        position = 0.5*(AssignedMuonTrack.GetEndPos() + AssignedProtonTrack.GetStartPos());
    }
    else if (min_distance == start_end_distance){
        position = 0.5*(AssignedMuonTrack.GetStartPos() + AssignedProtonTrack.GetEndPos());
    }
    else if (min_distance == end_end_distance){
        position = 0.5*(AssignedMuonTrack.GetEndPos() + AssignedProtonTrack.GetEndPos());
    }
    
    // then, flip the tracks accordingly
    if ( (AssignedMuonTrack.GetEndPos() - position).Mag() < (AssignedMuonTrack.GetStartPos() - position).Mag() ){
        // Debug(4,"Flipping muon track");
        AssignedMuonTrack.FlipTrack();
    }
    if ( (AssignedProtonTrack.GetEndPos() - position).Mag() < (AssignedProtonTrack.GetStartPos() - position).Mag() ){
        // Debug(4,"Flipping proton track");
        AssignedProtonTrack.FlipTrack();
    }
    
    // -- - --- -- -- --- - - -- -- -- - -- - -- -- -- - -- -- - - -- - - -- - -- - -- - - - -- - - -- - - -
    // muon angle is better reconstructed than proton angle
    // since the muon is longer and more 'pronounced' in the detector.
    // hence, we can use the muon angle to correct the proton angle:
    // flip proton track - based on \theta_muon-\theta_proton correlation
    // the MC correlation is a band around
    // 𝜽(p) = -𝜽(µ)/𝛑 + 1
    // so if 𝜽(p) is too far from this correlation we can flip the p-track
    if (fabs( AssignedProtonTrack.GetTheta() - (-AssignedMuonTrack.GetTheta()/PI + 1.)) > 1.){
        AssignedProtonTrack.FlipTrack();
    }
    // -- - --- -- -- --- - - -- -- -- - -- - -- -- -- - -- -- - - -- - - -- - -- - -- - - - -- - - -- - - -
}

//
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void pairVertex::SetEDepAroundVertex(){
//    
//    dqdx_around_vertex_non_tracks_associated = dqdx_around_vertex_tracks_associated = dqdx_around_vertex = -10000;
//    
//    dqdx_around_vertex = AssignedMuonTrack.dqdx_around_start_total + AssignedProtonTrack.dqdx_around_start_total;
//    
//    dqdx_around_vertex_tracks_associated = AssignedMuonTrack.dqdx_around_start_track_associated_total + AssignedProtonTrack.dqdx_around_start_track_associated_total;
//    
//    if(dqdx_around_vertex_tracks_associated!=0 && dqdx_around_vertex!=0
//       && dqdx_around_vertex_tracks_associated!=-10000 && dqdx_around_vertex!=-10000){
//        dqdx_around_vertex_non_tracks_associated = dqdx_around_vertex - dqdx_around_vertex_tracks_associated;
//    }
//    else {
//        dqdx_around_vertex_non_tracks_associated = -10000;
//    }
//}
//
//
//



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::SetReconstructedMomenta( float PmuFromRange, float PpFromRange ){
    // reconstruct the µ and p momenta by using the minimal features possible
    // theta / phi of the reconstructed track
    // and the momentum from stopping range
    // at a later stage we can maybe use calorimetery or multiple Coulomb scattering
    // p
    reco_Pp_3momentum = PpFromRange;
    reco_Pp_3vect.SetMagThetaPhi( reco_Pp_3momentum , AssignedProtonTrack.GetTheta() , AssignedProtonTrack.GetPhi() );
    reco_Pp.SetVectMag ( reco_Pp_3vect , 0.9385 );
    
    // µ
    reco_Pmu_3momentum = PmuFromRange;
    reco_Pmu_3vect.SetMagThetaPhi( reco_Pmu_3momentum , AssignedMuonTrack.GetTheta() , AssignedMuonTrack.GetPhi() );
    reco_Pmu.SetVectMag ( reco_Pmu_3vect , 0.1056 );
    
    //    reco_Pmu_3vect_mcsllhd.SetMagThetaPhi( AssignedMuonTrack.mommsllhd , AssignedMuonTrack.theta , AssignedMuonTrack.phi );
    //    reco_Pmu_mcsllhd.SetVectMag ( reco_Pmu_3vect_mcsllhd , 0.1056 );
    
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::ReconstructBeam(){
    // reconstruct the beam by using the reconstructed Pµ and Pp
    reco_BeamE = reco_Pmu.E() + (reco_Pp.E() - reco_Pp.Mag()) + 0.040 ; // Eµ + Tp + Sn + T(A-1)
    reco_Pnu = TLorentzVector( 0 , 0 , reco_BeamE , reco_BeamE );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::ReconstructKinematics(){
    // [http://pdg.lbl.gov/2014/reviews/rpp2014-rev-structure-functions.pdf]
    
    
    
    // reconstruct the momentum transfer from minimal features of CC1p and stopping range
    reco_q = reco_Pnu - reco_Pmu;
    reco_omega = reco_q.E();
    
    // reconstructed 𝜃(p,q) based on these minimal features
    reco_theta_pq = r2d * reco_Pp.Vect().Angle( reco_q.Vect() );
    reco_p_over_q = reco_Pp.P()/reco_q.P();
    reco_Q2 = - reco_q.Mag2();
    
    // kinematics
    reco_Xb = reco_Q2 / (2*0.939*reco_q.E());
    reco_n_miss = reco_Pp - reco_q;
    reco_y = reco_omega/reco_Pnu.E();
    reco_s = reco_Q2/(reco_Xb*reco_y) + 0.939*0.939 + 0.106*0.106;
    
    // LC momentum fraction
    reco_alpha_p = (reco_Pp.E()-reco_Pp.Pz())/0.931;
    reco_alpha_mu = (reco_Pmu.E()-reco_Pmu.Pz())/0.931;
    reco_alpha_q = (reco_q.E() - reco_q.Pz())/0.931;
    reco_alpha_miss = reco_alpha_p - reco_alpha_q;

    // invariant mass of the interaction
    reco_W2 = 0.939*(0.939 + 2*(reco_Pnu.E() - reco_Pmu.E())) - 4*reco_Pnu.E()*reco_Pmu.E()*(1.-cos(reco_Pmu.Theta()));

    
    // truth information for MC
    if (genie_interaction.GetNprotons() > 0){
        TLorentzVector truth_Plead = genie_interaction.GetProtonsMomenta().at(0);
        truth_alpha_p = (truth_Plead.E() - truth_Plead.Pz())/0.931;
    }
    TLorentzVector truth_muon = genie_interaction.GetLeptonMomentum();
    truth_alpha_mu = (truth_muon.E() - truth_muon.Pz())/0.931;
    
    TLorentzVector truth_q = genie_interaction.GetMomentumTransfer();
    truth_alpha_q = (truth_q.E() - truth_q.Pz())/0.931;
    truth_alpha_miss = truth_alpha_p - truth_alpha_q;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::SetReconstructedFeatures( float PmuFromRange, float PpFromRange ){
    
    // reconstructed distance between µ and p
    reco_mu_p_distance = AssignedMuonTrack.ClosestDistanceToOtherTrack( AssignedProtonTrack );
    
    SetReconstructedMomenta( PmuFromRange, PpFromRange );
    ReconstructBeam();
    ReconstructKinematics();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float pairVertex::GetTruthDeltaPhi () const {
    // return the truth \Delta \phi between tracks in degrees
    float muon_truth_phi = AssignedMuonTrack.GetTruthMomentum().Phi();
    float proton_truth_phi = AssignedProtonTrack.GetTruthMomentum().Phi();
    return r2d*fabs( muon_truth_phi - proton_truth_phi );
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::AssociateHitsToTracks (std::vector<hit> hits) {
    // Aug-3,2017
    // plug the proton hits to the proton (candidate) and the muon (candidate) hits
    
    for (auto hit:hits){
        
        int plane = hit.GetPlane();
        
        if (hit.GetTrackKey()== AssignedMuonTrack.GetTrackID() ){
            hits_muon[plane].push_back(hit);
        }
        else if (hit.GetTrackKey()== AssignedProtonTrack.GetTrackID() ){
            hits_proton[plane].push_back(hit);
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
float pairVertex::GetRdQaroundVertex (int plane, int Nwires, int Nticks , std::vector<hit> hits) const {
    // Aug-3,2017
    // get the ratio of tracks-charge deposited to total-charge deposited
    // in a box of N(wires) x N(time-ticks) around the vertex in plane i=0,1,2
    // a mal-function wil return -1
    // input:
    // plane, N(wires) & N(time-ticks) for the box, hits in event
    if (plane<0 || plane>2) return -1;
    
    box VertexBox( (Int_t)vertex_wire[plane] - Nwires , (Int_t)vertex_time[plane] - Nticks,
                   (Int_t)vertex_wire[plane] + Nwires , (Int_t)vertex_time[plane] + Nticks );
    
    
    float Qtotal = 0.0 , Qtracks = 0.0;
    for (auto hit:hits){
        
        if (hit.InPlane(plane) && hit.InBox(VertexBox)){
            
            Qtotal += hit.GetCharge();
            
            // the hit belongs to the muon/proton track add its charge to the tracks charge
            if ( hit.GetTrackKey()== AssignedMuonTrack.GetTrackID()  || hit.GetTrackKey()== AssignedProtonTrack.GetTrackID() ){
                Qtracks += hit.GetCharge();
            }
        }
    }
    if (fabs(Qtracks)>0) return (Qtracks/Qtotal);
    else return -1;
}




//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//bool pairVertex::BuildROI(int plane){
//    
//    // for GENIE vertices in which only one track was reconstructed, build the roi based on this track
//    if (muonTrackReconstructed && !protonTrackReconstructed) {AssignedProtonTrack=AssignedMuonTrack;}
//    if (!muonTrackReconstructed && protonTrackReconstructed) {AssignedMuonTrack=AssignedProtonTrack;}
//    // -------------------------------------------------------------------------------------------------
//    
//    int wire_min = std::min( {AssignedMuonTrack.roi[plane].start_wire , AssignedMuonTrack.roi[plane].end_wire , AssignedProtonTrack.roi[plane].start_wire, AssignedProtonTrack.roi[plane].end_wire} );
//    int wire_max = std::max( {AssignedMuonTrack.roi[plane].start_wire , AssignedMuonTrack.roi[plane].end_wire , AssignedProtonTrack.roi[plane].start_wire, AssignedProtonTrack.roi[plane].end_wire} );
//    int time_min = std::min( {AssignedMuonTrack.roi[plane].start_time , AssignedMuonTrack.roi[plane].end_time , AssignedProtonTrack.roi[plane].start_time, AssignedProtonTrack.roi[plane].end_time} );
//    int time_max = std::max( {AssignedMuonTrack.roi[plane].start_time , AssignedMuonTrack.roi[plane].end_time , AssignedProtonTrack.roi[plane].start_time, AssignedProtonTrack.roi[plane].end_time} );
//    
//    roi[plane] = box( wire_min - 10 , time_min - 30 , wire_max + 10 , time_max + 30 );
//
//    switch (plane) {
//        case 0:
//            roi_u = roi[plane];
//            break;
//        case 1:
//            roi_v = roi[plane];
//            break;
//        case 2:
//        default:
//            roi_y = roi[plane];
//            break;
//    }
//    return true;
//}

//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//bool pairVertex::BuildLocationInPlane(int plane){
//    
//    mu_start_wire[plane] = mu_start_time[plane] = p_start_wire[plane] = p_start_time[plane] = 0;
//    mu_end_wire[plane] = mu_end_time[plane] = p_end_wire[plane] = p_end_time[plane] = 0;
//
//    
//    switch (plane) {
//        case 0:
//            mu_start_wire[plane] = AssignedMuonTrack.start_wire_u;
//            mu_start_time[plane] = AssignedMuonTrack.start_time_u;
//            mu_end_wire[plane] = AssignedMuonTrack.end_wire_u;
//            mu_end_time[plane] = AssignedMuonTrack.end_time_u;
//            
//            p_start_wire[plane] = AssignedProtonTrack.start_wire_u;
//            p_start_time[plane] = AssignedProtonTrack.start_time_u;
//            p_end_wire[plane] = AssignedProtonTrack.end_wire_u;
//            p_end_time[plane] = AssignedProtonTrack.end_time_u;
//            
//            break;
//        case 1:
//            mu_start_wire[plane] = AssignedMuonTrack.start_wire_v;
//            mu_start_time[plane] = AssignedMuonTrack.start_time_v;
//            mu_end_wire[plane] = AssignedMuonTrack.end_wire_v;
//            mu_end_time[plane] = AssignedMuonTrack.end_time_v;
//            
//            p_start_wire[plane] = AssignedProtonTrack.start_wire_v;
//            p_start_time[plane] = AssignedProtonTrack.start_time_v;
//            p_end_wire[plane] = AssignedProtonTrack.end_wire_v;
//            p_end_time[plane] = AssignedProtonTrack.end_time_v;
//
//            break;
//        case 2:
//        default:
//            mu_start_wire[plane] = AssignedMuonTrack.start_wire_y;
//            mu_start_time[plane] = AssignedMuonTrack.start_time_y;
//            mu_end_wire[plane] = AssignedMuonTrack.end_wire_y;
//            mu_end_time[plane] = AssignedMuonTrack.end_time_y;
//            
//            p_start_wire[plane] = AssignedProtonTrack.start_wire_y;
//            p_start_time[plane] = AssignedProtonTrack.start_time_y;
//            p_end_wire[plane] = AssignedProtonTrack.end_wire_y;
//            p_end_time[plane] = AssignedProtonTrack.end_time_y;
//            
//            break;
//    }
//    
//    // first fix the position of the vertex
//    TVector2 mu_start( mu_start_wire[plane] , mu_start_time[plane] );
//    TVector2  mu_end( mu_end_wire[plane]    , mu_end_time[plane] );
//    TVector2 p_start(  p_start_wire[plane]  , p_start_time[plane] );
//    TVector2   p_end( p_end_wire[plane]     , p_end_time[plane] );
//    
//    float d_start_start = WireTimeDistance( mu_start.X() , mu_start.Y() , p_start.X() , p_start.Y() );
//    float   d_end_start = WireTimeDistance( mu_end.X()   , mu_end.Y()   , p_start.X() , p_start.Y() );
//    float   d_start_end = WireTimeDistance( mu_start.X() , mu_start.Y() , p_end.X()   , p_end.Y() );
//    float     d_end_end = WireTimeDistance( mu_end.X()   , mu_end.Y()   , p_end.X()   , p_end.Y() );
//    float d_min = std::min({ d_start_start , d_end_start , d_start_end , d_end_end });
//    
//    if (d_min==d_start_start){
//        vertex_wire[plane] = 0.5*(mu_start.X() + p_start.X());
//        vertex_time[plane] = 0.5*(mu_start.Y() + p_start.Y());
//    }
//    else if (d_min==d_end_start){
//        vertex_wire[plane] = 0.5*(mu_end.X() + p_start.X());
//        vertex_time[plane] = 0.5*(mu_end.Y() + p_start.Y());
//    }
//    else if (d_min==d_start_end){
//        vertex_wire[plane] = 0.5*(mu_start.X() + p_end.X());
//        vertex_time[plane] = 0.5*(mu_start.Y() + p_end.Y());
//    }
//    else if (d_min==d_end_end){
//        vertex_wire[plane] = 0.5*(mu_end.X() + p_end.X());
//        vertex_time[plane] = 0.5*(mu_end.Y() + p_end.Y());
//    }
//    // build a small ROI of 20 wires x 40 time ticks around the vertex
//    Roi_20x40_AroundVertex[plane] = box( vertex_wire[plane]-20, vertex_time[plane]-40, vertex_wire[plane]+20, vertex_time[plane]+40  );
//    
//    // now flip the tracks that need to be flipped
//    if ( WireTimeDistance( vertex_wire[plane] , vertex_time[plane] , mu_start_wire[plane] , mu_start_time[plane] ) > WireTimeDistance( vertex_wire[plane] , vertex_time[plane] , mu_end_wire[plane] , mu_end_time[plane] ) ){
//        int tmp_wire = mu_start_wire[plane];
//        mu_start_wire[plane] = mu_end_wire[plane];
//        mu_end_wire[plane] = tmp_wire;
//        int tmp_time = mu_start_time[plane];
//        mu_start_time[plane] = mu_end_time[plane];
//        mu_end_time[plane] = tmp_time;
//    }
//    if ( WireTimeDistance( vertex_wire[plane] , vertex_time[plane] , p_start_wire[plane] , p_start_time[plane] ) > WireTimeDistance( vertex_wire[plane] , vertex_time[plane] , p_end_wire[plane] , p_end_time[plane] ) ){
//        int tmp_wire = p_start_wire[plane];
//        p_start_wire[plane] = p_end_wire[plane];
//        p_end_wire[plane] = tmp_wire;
//        int tmp_time = p_start_time[plane];
//        p_start_time[plane] = p_end_time[plane];
//        p_end_time[plane] = tmp_time;
//    }
//    
//    mu_angle[plane] = WireTimeAngle( (float)mu_start_wire[plane] ,  (float)mu_start_time[plane]
//                                    , (float)mu_end_wire[plane] , (float)mu_end_time[plane] );
//    AssignedMuonTrack.WireTimeAngle[plane] = mu_angle[plane];
//    
//    p_angle[plane] = WireTimeAngle( (float)p_start_wire[plane] ,  (float)p_start_time[plane]
//                                    , (float)p_end_wire[plane] , (float)p_end_time[plane] );
//    AssignedProtonTrack.WireTimeAngle[plane] = p_angle[plane];
//    
//    return true;
//    
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void pairVertex::Print(bool DoPrintTracks) const {
    
    cout << "\033[35m" << "~~~~~~~~~~~~~~\n vertex " << vertex_id << "\n~~~~~~~~~~~~~~"<< "\033[0m" << endl;
    if (Is1mu1p) Printf("1mu-1p vertex!");
    // SHOW3( run , subrun , event );
    SHOWTVector3( position );
    
    for (auto t: tracks) {
        Printf("track %d (pdg %d), vertex distance %.1f cm",t.GetTrackID(),t.GetMCpdgCode(),t.DistanceFromPoint( position ));
    }
    
    // tracks
    if (!tracks.empty()){
        
        if (DoPrintTracks) {
            cout << "\033[33m" << tracks.size() << " tracks in vertex " << "\033[31m" << endl;
            for (auto t: tracks) {
                if (t.GetTrackID()!=-100)   t.Print( true );
                else                    Printf("unidentified/non-reconstructed track. continuing...");
            }
        }
        
        // inter-tracks distances
        cout << "\033[33m" << tracks.size()*tracks.size() << " inter-tracks distances " << "\033[31m" << endl;
        for(auto vec : tracks_distances) {
            for(auto x : vec) std::cout << setprecision(2) << x << "\t";
            std::cout << std::endl;
        }
    }
    else {
        cout << "\033[33m" << "no reconstructed tracks in vertex " << "\033[31m" << endl;
    }
    
    SHOW( IsGENIECC_1p_200MeVc_0pi );
    if (IsGENIECC_1p_200MeVc_0pi){
        cout << "This vertex is a GENIE true CC1p0π" << endl;
        SHOW3( Is_mu_TrackReconstructed, Is_p_TrackReconstructed, IsVertexReconstructed );
        
        // neutrino
        SHOWTLorentzVector( reco_Pnu );
        
        // muon
        SHOWTLorentzVector( genie_interaction.GetLeptonMomentum() );
        if (!tracks.empty()) SHOWTLorentzVector( reco_Pmu );
        
        // proton
        if (genie_interaction.GetNprotons()) SHOWTLorentzVector( genie_interaction.GetProtonsMomenta().at(0) );
        if (!tracks.empty()) SHOWTLorentzVector( reco_Pp );
        
        // theta (p,q)
        SHOW(genie_interaction.Get_theta_pq());
        if (!tracks.empty()) SHOW(reco_theta_pq);
    }
    
    // GENIE interaction features
    PrintLine();
    cout << "MC information: " << endl;
    cout << "truth topology: " << TruthTopologyString << endl;
    if (Is1mu1p){
        if ( genie_interaction.GetMCeventID() == -9999 ) Printf("vertex does not have a matching GENIE interaction");
        else genie_interaction.Print();
        PrintLine();
        if ( closest_genie_interaction.GetMCeventID() == -9999 ) Printf("vertex does not have a matching closest GENIE interaction around");
        else {
            Printf("closest GENIE interaction is %.1f cm away from reconstructed vertex:",(closest_genie_interaction.GetVertexPosition() - position).Mag());
            closest_genie_interaction.Print();
        }
        PrintLine();
    }
    
//    if (!tracks.empty()) for (int plane = 0 ; plane<3 ; plane ++ ) {printf("ROI in plane %d: ",plane); roi[plane].Print();}
}


#endif

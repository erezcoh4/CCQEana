////////////////////////////////////////////////////////////////////////
// Class:       ErezCCQEAnalyzer
// Plugin Type: analyzer (art v2_05_00)
// File:        ErezCCQEAnalyzer_module.cc
//
// Generated at Wed Jul 12 16:10:16 2017 by Erez Cohen using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larsim/MCCheater/BackTracker.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCFlux.h"


// ROOT includes
#include "TTree.h"

// C/C++ libraries
#include <memory>
#include <utility>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>


// my pandoraNu track...
#include "uboone/ErezCCQEana/MyObjects/PandoraNuTrack.h"
#include "uboone/ErezCCQEana/MyObjects/hit.h"
#include "uboone/ErezCCQEana/MyObjects/box.h"
#include "uboone/ErezCCQEana/MyObjects/GENIEinteraction.h"
#include "uboone/ErezCCQEana/MyObjects/pairVertex.h"

// constants
constexpr int debug          = 1;
constexpr int kMaxTrack      = 1000;  //maximum number of tracks
constexpr int kMaxHits       = 40000; //maximum number of hits;
constexpr int kMaxTruth      = 100;
constexpr int kMaxNgenie     = 100;
constexpr float kMaxInterTrackDistance = 11; // 11 cm between tracks - maximal distance for clustering
constexpr float EPSILON      = 0.1;   // tollerance for equations


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace ub { class ErezCCQEAnalyzer; }





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ub::ErezCCQEAnalyzer : public art::EDAnalyzer {
public:
    explicit ErezCCQEAnalyzer(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    ErezCCQEAnalyzer(ErezCCQEAnalyzer const &) = delete;
    ErezCCQEAnalyzer(ErezCCQEAnalyzer &&) = delete;
    ErezCCQEAnalyzer & operator = (ErezCCQEAnalyzer const &) = delete;
    ErezCCQEAnalyzer & operator = (ErezCCQEAnalyzer &&) = delete;
    
    // Required functions.
    void                   analyze (art::Event const & e) override;
    
    // Selected optional functions.
    void                  beginJob () override;
    void               reconfigure (fhicl::ParameterSet const& p) override;
    
    // functionallity
    void         ConstructVertices ();
    void   ClusterTracksToVertices ();
    void           AnalyzeVertices ();
    void          FindPairVertices ();
    void               TagVertices ();
    void          PrintInformation ();
    bool    TrackAlreadyInVertices (int ftrack_id);
    void       HeaderVerticesInCSV ();
    void       StreamVerticesToCSV ();

    
    
    
    
private:
    
    // Declare member data here.
    void ResetVars();
    
    // Declare member data here.
    TTree *fTree;
    Int_t run, subrun, event;
    
    
    short   isdata;
    
    bool    MCmode;
    
    int     Ntracks;                // number of reconstructed tracks
    int     Nhits , Nhits_stored;   // number of recorded hits in the event
    int     Nvertices;
    int     vertices_ctr;
    
    // my objects
    std::vector<PandoraNuTrack>     tracks;
    std::vector<hit>                hits;
    std::vector<GENIEinteraction>   genie_interactions;
    std::vector<pairVertex>         vertices;
    
    //Module labels to get data products
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fFlashModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fGenieGenModuleLabel;

    //mctruth information
    Int_t    mcevts_truth;    //number of neutrino Int_teractions in the spill
    
    
    // time stamp
    std::chrono::time_point<std::chrono::system_clock> start_ana_time, end_ana_time;
    
    // output csv file of vertices
    ofstream vertices_file;

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ub::ErezCCQEAnalyzer::ErezCCQEAnalyzer(fhicl::ParameterSet const & p):EDAnalyzer(p){
    reconfigure(p);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::analyze(art::Event const & evt){
    
    ResetVars();
    
    art::ServiceHandle<geo::Geometry> geom;
    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<cheat::BackTracker> bt;
    //    const sim::ParticleList& plist = bt->ParticleList();
    isdata = evt.isRealData();
    run = evt.run(); subrun = evt.subRun(); event = evt.id().event();
    
    // * hits
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;
    if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

    // * tracks
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
    
    // * associations
    art::FindManyP<recob::Hit> fmth(trackListHandle, evt, fTrackModuleLabel);
    art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

    
    // ----------------------------------------
    // hits information
    // ----------------------------------------
    Nhits = hitlist.size();
    Nhits_stored = std::min(Nhits, kMaxHits);
    for (int i = 0; i < Nhits_stored ; ++i){//loop over hits
        
        hit fhit(
                  hitlist[i]->WireID().Plane    // hit plane
                  ,hitlist[i]->WireID().Wire    // hit wire
                  ,i                            // hit id
                  ,hitlist[i]->PeakTime()       // hit peak time
                  ,hitlist[i]->Integral()       // hit charge
                  );
        
        hits.push_back( fhit );
        
    }

    // ----------------------------------------
    // tracks information
    // ----------------------------------------
    Ntracks = tracklist.size();
    for(int i=0; i < std::min(int(tracklist.size()),kMaxTrack); ++i ){
        recob::Track::Point_t start_pos, end_pos;
        std::tie( start_pos, end_pos ) = tracklist[i]->Extent();
        
        PandoraNuTrack track(
                             run , subrun, event        // r/s/e
                             ,tracklist[i]->ID()        // track id
                             ,tracklist[i]->Length()    // length
                             ,tracklist[i]->Theta()     // polar angle
                             ,tracklist[i]->Phi()       // azimuthal angle
                             ,TVector3(start_pos.X(),start_pos.Y(),start_pos.Z())   // start position
                             ,TVector3(end_pos.X(),end_pos.Y(),end_pos.Z())         // end position
                            );
        
        double StartLoc[3] = {start_pos.X(), start_pos.Y(), start_pos.Z()};
        
        // U / V / Y coordinates
        for (int plane = 0; plane < 3; plane++){
            // Get the PlaneID in (geo::PlaneID) type
            raw::ChannelID_t channel = geom -> NearestChannel( StartLoc , plane );
            std::vector<geo::WireID> wires = geom -> ChannelToWire(channel);
            const geo::PlaneID & planeID = wires[0].planeID();
            // get start point
            Int_t start_wire = (Int_t) ( geom->WireCoordinate( start_pos.Y() , start_pos.Z() ,  planeID ) );
            Int_t start_time = (Int_t) ( detprop->ConvertXToTicks( start_pos.X() , planeID ) ) ;
            // get end point
            Int_t end_wire = (Int_t) ( geom->WireCoordinate( end_pos.Y() , end_pos.Z() ,  planeID ) );
            Int_t end_time = (Int_t) ( detprop->ConvertXToTicks( end_pos.X() , planeID ) ) ;
            // plug into the track
            track.SetStartEndPlane( plane , start_wire , start_time , end_wire , end_time );
        }

        
        // Hits-Tracks association
        if (fmth.isValid()){
            std::vector< art::Ptr<recob::Hit> > vhit = fmth.at(i);
            for (size_t h = 0; h < vhit.size(); ++h){
                if (vhit[h].key()<kMaxHits){
                    hits.at( vhit[h].key() ).SetTrackKey( tracklist[i].key() );
                }
            }
        }
        
        // PIDa and calorimetric KE
        if (fmcal.isValid()){
            unsigned maxnumhits = 0;
            std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
            for (auto const& calo : calos){
                if (calo->PlaneID().isValid){
                    int plane = calo->PlaneID().Plane;
                    
                    // get the calorimetric kinetic energy of the track
                    track.SetCaloKEPerPlane( plane , calo->KineticEnergy() );
                    
                    // select the best plane as the one with the maximal number of charge deposition points
                    if (calo->dEdx().size() > maxnumhits){
                        maxnumhits = calo->dEdx().size();
                        track.SetBestPlane ( plane );
                        track.SetMaxNHits ( maxnumhits );
                    }
                    // build the PIDa as a fit the reduced Bethe Bloch
                    // dE/dx = A * R^{0.42}
                    double pida = 0;
                    int used_trkres = 0;
                    for (size_t ip = 0; ip < calo->dEdx().size(); ++ip){
                        if (calo->ResidualRange()[ip]<30){
                            pida += calo->dEdx()[ip]*pow( calo->ResidualRange()[ip],0.42);
                            ++used_trkres;
                        }
                    }
                    if (used_trkres) pida /= used_trkres;
                    track.SetPIDaPerPlane( plane , pida );
                }
            }
            track.SetPIDa();
        }
        
        
        
        // MC information
        if (!isdata&&fmth.isValid()){
            // Find true track for each reconstructed track
            int TrackID = 0;
            std::vector< art::Ptr<recob::Hit> > allHits = fmth.at(i);
            
            std::map<int,double> trkide;
            for(size_t h = 0; h < allHits.size(); ++h){
                art::Ptr<recob::Hit> hit = allHits[h];
                std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
                for( size_t e = 0; e < TrackIDs.size(); ++e ){
                    trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
                }
            }
            // Work out which IDE despoited the most charge in the hit if there was more than one.
            double max_Edep = -1;
            for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
                if ((ii->second)>max_Edep){
                    max_Edep = ii->second;
                    TrackID = ii->first;
                }
            }
            // Now have trackID, so get PDG code and T0 etc.
            const simb::MCParticle *particle = bt->TrackIDToParticle( TrackID );
            if (particle){
                track.SetMCpdgCode( particle->PdgCode() );
                track.SetTruthStartPos( TVector3(particle->Vx() , particle->Vy() , particle->Vz()) );
                track.SetTruthEndPos( TVector3(particle->EndX() , particle->EndY() , particle->EndZ()) );
                track.SetTruthMomentum( particle -> Momentum() );
            }//if (particle)
        }//MC
        
        tracks.push_back( track );
    }
    
    
    
    // ----------------------------------------
    // MC truth information
    // ----------------------------------------
    if (!isdata){
        
        // * MC truth information
        art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
        std::vector<art::Ptr<simb::MCTruth> > mclist;
        if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);
        
        art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
        std::vector<art::Ptr<simb::MCFlux> > fluxlist;
        if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
        art::fill_ptr_vector(fluxlist, mcfluxListHandle);
        
        
        mcevts_truth = mclist.size();
        if (mcevts_truth){
            MCmode = true;
            
            for( int mc_evend_id = 0; (mc_evend_id < mcevts_truth) && (mc_evend_id < kMaxTruth) ; mc_evend_id++ ){
                art::Ptr<simb::MCTruth> mctruth = mclist[mc_evend_id];
                
                if (mctruth->Origin() == simb::kBeamNeutrino){
                    
                    GENIEinteraction genie_interaction( run , subrun , event , mc_evend_id );
                    genie_interaction.SetNuPDG( mctruth->GetNeutrino().Nu().PdgCode() );
                    genie_interaction.SetKinematics(
                                                    mctruth->GetNeutrino().QSqr() // Q2
                                                    ,mctruth->GetNeutrino().W()     // W
                                                    ,mctruth->GetNeutrino().X()     // Bjorken x
                                                    ,mctruth->GetNeutrino().Y()     // y
                                                    ,mctruth->GetNeutrino().CCNC()  // CC=0/NC=1
                                                    ,mctruth->GetNeutrino().Mode()  // QE=0
                                                    );
                    genie_interaction.SetVertexPosition( TVector3(mctruth->GetNeutrino().Nu().Vx()
                                                                  ,mctruth->GetNeutrino().Nu().Vy()
                                                                  ,mctruth->GetNeutrino().Nu().Vz()));
                    
                    genie_interaction.SetNuMomentum( mctruth->GetNeutrino().Nu().Momentum() );
                    genie_interaction.SetLeptonMomentum( mctruth->GetNeutrino().Lepton().Momentum() );
                    genie_interaction.SetMomentumTransfer();
                    
                    // all other particles in the interaction
                    Int_t NgenieParticles = (Int_t)mctruth->NParticles();
                    if ( NgenieParticles ){
                        for( int iPart = 0; iPart < std::min( NgenieParticles , kMaxNgenie ); iPart++ ){
                            const simb::MCParticle& part( mctruth->GetParticle(iPart) );
                            // add a primary to genie-interaction
                            genie_interaction.AddPrimary( part.PdgCode()    // pdg code
                                                         ,part.Momentum()   // 4-momentum
                                                         ,part.StatusCode() // status code
                                                         ,part.Mother()     // mother
                                                         );
                            
                            // match the primary particle with a track
                            for (auto & track : tracks){
                                if (
                                    ( part.PdgCode() == track.GetMCpdgCode() )
                                    &&
                                    ( (part.Momentum() - track.GetTruthMomentum()).Mag() < EPSILON )
                                    &&
                                    ( (genie_interaction.GetVertexPosition() - track.GetTruthStartPos()).Mag() < EPSILON )
                                    ) {
                                    // Printf("found a primary-track match! plugging mc_evend_id=%d into track %d",mc_evend_id,track.GetTrackID());
                                    genie_interaction.AddTrack ( track );
                                    track.SetMCeventID( mc_evend_id );
                                }
                            }

                        } // for particle
                    }
                    genie_interaction.SortNucleons();
                    genie_interaction.ComputePmissPrec();
                    genie_interaction.FindCC1p200MeVc0pi();

                    genie_interactions.push_back( genie_interaction );
                    
                }//mctruth->Origin()
            }
        }
    }//is neutrino


    // ----------------------------------------
    // event topology (my-vertex....)
    // ----------------------------------------
    ConstructVertices();
    
    
    PrintInformation();
    fTree -> Fill();
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::ConstructVertices(){
    
    // cluster all tracks at close proximity to vertices
    ClusterTracksToVertices();
    // analyze these vertices: inter-tracks distances, angles...
    AnalyzeVertices();
    // retain only vertices with pairs of 2-tracks at close proximity
    FindPairVertices();
    // if its a MC event, tag the vertex by their MC information
    TagVertices();
    // output to csv file
    StreamVerticesToCSV();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::ClusterTracksToVertices(){
    // July-25, 2017
    // cluster all tracks at close proximity to vertices
    bool    FoundCloseTracks , AlreadySetPosition;
    float   closest_distance_ij;
    TVector3 vertex_position;
    
    for (int i=0; i < Ntracks; i++){
        
        // if (!tracks[i].IsFullyContained) continue;
        if (!tracks[i].IsTrackContainedSoft()) continue;
        
        // skip if track was clustered to a vertex by in one of the previous loop steps
        if ( TrackAlreadyInVertices( tracks[i].GetTrackID() )) continue;
        
        pairVertex vertex( run, subrun, event , vertices.size() );
        vertex.AddTrack( tracks[i] );
        
        FoundCloseTracks = AlreadySetPosition = false;
        
        for ( int j=0 ; j < Ntracks ; j++ ){ // i+1
            
            // if (!tracks[j].IsFullyContained) continue;
            if (tracks[j].IsTrackContainedSoft() && j!=i){
                
                // if this is the first time we go over these two tracks
                // and they are close enough to define a vertex,
                // we also define the position of their mutual vertex
                if (!AlreadySetPosition){
                    
                    // two close tracks (at a separation distance smaller that max_mu_p_distance)
                    std::string StartOrEnd = "None";
                    closest_distance_ij = tracks[i].ClosestDistanceToOtherTrack(tracks[j],&StartOrEnd);
                    
                    if ( closest_distance_ij < kMaxInterTrackDistance ){
                        
                        vertex.AddTrack( tracks[j] );
                        FoundCloseTracks = true;
                        
                        if (StartOrEnd.compare("Start")==0)     vertex_position = tracks[i].GetStartPos();
                        else if (StartOrEnd.compare("End")==0)  vertex_position = tracks[i].GetEndPos() ;
                        else                                    vertex_position = TVector3(-1000,-1000,-1000) ;
                        
                        vertex.SetPosition( vertex_position );
                        AlreadySetPosition = true;
                    }
                }
                
                // else, namely if we have already clustered a vertex
                // and positioned it in space,
                // we only need to check wether the new track (j) is close enough to this vertex
                // to be also associated with it
                else {
                    // Debug(3, Form("cheking track %d distance from vertex %d ",tracks[j].track_id,c_vertex.vertex_id));
                    if ( tracks[j].DistanceFromPoint(vertex.GetPosition()) < kMaxInterTrackDistance ){
                        
                        // SHOWTVector3(c_vertex.position);
                        // Printf("track %d close enough...",tracks[j].track_id);
                        // SHOW(tracks[j].DistanceFromPoint(c_vertex.position));
                        
                        vertex.AddTrack( tracks[j] );
                    }
                }
            }
        }
        
        
        // if this is an MC event,
        // match a GENIE interaction to the vertex
        if (MCmode){
            if (FoundCloseTracks) {

                // 1st (and best) method: match the GENIE interaction by the mc event-id
                // use the mc event-id of the first/second track. This is the mc event-id of the proper GENIE interaction
                if (vertex.GetTracks().size()>1){
                    if ( vertex.GetTracks().at(0).GetMCeventID() == vertex.GetTracks().at(1).GetMCeventID() ){
                        vertex.SetGENIEinfo( genie_interactions.at( vertex.GetTracks().at(0).GetMCeventID() ) );
                    }
                }
                // the problem with this method is that not allways the two tracks came from
                // the same GENIE interaction
                // for example, happenings in which one track came from a GENIE int. and the other from cosmic...
                // so,
                // 2nd method: look for the closest GENIE interaction in the event, spatially
                GENIEinteraction closest_genie;
                bool MatchedGENIEinteraction = false;
                float closest_genie_interaction_distance = 10000; // [cm]
                for (auto genie_interaction : genie_interactions){
                    float genie_distance = (genie_interaction.GetVertexPosition() - vertex.GetPosition()).Mag();
                    if ( genie_distance < closest_genie_interaction_distance ){
                        closest_genie_interaction_distance = genie_distance;
                        closest_genie = genie_interaction;
                        MatchedGENIEinteraction = true;
                        break;
                    }
                }
                if (MatchedGENIEinteraction){
                    vertex.SetClosestGENIE( closest_genie );
                }
            }
        }
        // plug into vertices list
        vertices.push_back( vertex );
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool ub::ErezCCQEAnalyzer::TrackAlreadyInVertices(int ftrack_id){
    for (auto v:vertices){
        if ( v.IncludesTrack( ftrack_id ) ) return true;
    }
    return false;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::AnalyzeVertices(){
    for (auto & v:vertices){
        // after fixing the vertext position, remove far tracks
        v.RemoveFarTracks( kMaxInterTrackDistance );
        // now sort the tracks
        v.SortTracksByPIDA ();
        v.SortTracksByLength ();
        // and the relations between the tracks
        // inter-track distances, delta-theta, delta-phi...
        v.SetTracksRelations ();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::FindPairVertices(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::TagVertices(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::HeaderVerticesInCSV(){
    
    vertices_ctr = 0;
    
    vertices_file
    << "run" << "," << "subrun" << "," << "event" << "," << "vertex_id" << ","
    << "x" << "," << "y" << "," << "z" << ","
    << "tracks" << ","
    << endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::StreamVerticesToCSV(){
    // July-25, 2017
    // whatever you add here - must add also in header - ub::ErezCCQEAnalyzer::HeaderVerticesInCSV()
    for (auto v:vertices){
        
        vertices_ctr++;
        
        vertices_file
        << v.GetRun() << "," << v.GetSubrun() << "," << v.GetEvent() << "," << v.GetVertexID() << ","
        << v.GetPosition().x() << "," << v.GetPosition().y() << "," << v.GetPosition().z() << ",";
        
        for ( auto track : v.GetTracks()){
            vertices_file << track.GetTrackID() << ";" ;
        }
        
        
        vertices_file << "," << endl;
        
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::PrintInformation(){
    
    SHOW( fTree->GetEntries() );
    SHOW3( run , subrun , event );
    
    if (MCmode){
        Printf("this is an MC event, with an MC information");
        if(!genie_interactions.empty()){
            cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << genie_interactions.size() << " genie interactions\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
            for (auto g: genie_interactions) {
                g.Print();
            }
        }
    }
    
    if(!tracks.empty()){
        cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << tracks.size() << " pandoraNu tracks\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
        for (auto t: tracks) {
            t.Print( true );
        }
    }
    
    
    if(!vertices.empty()){
        cout << "\033[36m" << "xxxxxxxxxxxxxx\n\n" << vertices.size() << " vertices\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
        for (auto v: vertices) {
            v.Print( );
        }
    }

    // time stamp
    PrintLine();
    end_ana_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_ana_time - start_ana_time;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end_ana_time);
    std::cout << "\033[33m"
    << "finished analysis of this event" << std::ctime(&end_time)
    << "time elapsed: " << elapsed_seconds.count() << "s"
    << "\033[31m" << endl;

    EndEventBlock();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::beginJob(){
    
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("eventsTree","analysis tree of events");
    fTree->Branch("run"     ,&run           ,"run/I");
    fTree->Branch("subrun"  ,&subrun        ,"subrun/I");
    fTree->Branch("event"   ,&event         ,"event/I");
    fTree->Branch("Ntracks" ,&Ntracks       ,"Ntracks/I");
    fTree->Branch("Nhits"   ,&Nhits_stored  ,"Nhits/I");
    fTree->Branch("Nvertices",&Nvertices    ,"Nvertices/I");
    
    // my objects
    fTree->Branch("tracks"              ,&tracks);
    fTree->Branch("hits"                ,&hits);
    fTree->Branch("genie_interactions"  ,&genie_interactions);
    fTree->Branch("vertices"            ,&vertices);

    
    // output csv file
    vertices_file.open("/uboone/data/users/ecohen/CCQEanalysis/csvFiles/ccqe_candidates/vertices.csv");
    HeaderVerticesInCSV();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::reconfigure(fhicl::ParameterSet const & p){
    fTrackModuleLabel = p.get< std::string >("TrackModuleLabel");
    fHitsModuleLabel  = p.get< std::string >("HitsModuleLabel");
    fGenieGenModuleLabel =  p.get< std::string >("GenieGenModuleLabel");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::ResetVars(){
    
    MCmode = false;
    run = subrun = event = -9999;
    Ntracks = Nhits = Nhits_stored = Nvertices = 0 ;
    tracks.clear();
    hits.clear();
    genie_interactions.clear();
    vertices.clear();
    
    // time-stamp
    start_ana_time = std::chrono::system_clock::now();
}



// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
DEFINE_ART_MODULE(ub::ErezCCQEAnalyzer)
// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
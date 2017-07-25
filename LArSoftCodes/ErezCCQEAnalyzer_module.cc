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


// my pandoraNu track...
#include "uboone/ErezCCQEana/MyObjects/PandoraNuTrack.h"
#include "uboone/ErezCCQEana/MyObjects/hit.h"
#include "uboone/ErezCCQEana/MyObjects/box.h"
#include "uboone/ErezCCQEana/MyObjects/GENIEinteraction.h"

// constants
constexpr int kMaxTrack      = 1000;  //maximum number of tracks
constexpr int kMaxHits       = 40000; //maximum number of hits;
constexpr int kMaxTruth      = 100;
constexpr int kMaxNgenie     = 100;

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
    void analyze(art::Event const & e) override;
    
    // Selected optional functions.
    void beginJob() override;
    void reconfigure (fhicl::ParameterSet const& p) override;
    
    // functionallity
    void PrintInformation ();
    
private:
    
    // Declare member data here.
    void ResetVars();
    
    // Declare member data here.
    TTree *fTree;
    Int_t run, subrun, event;
    
    
    short   isdata;
    
    int     Ntracks;                // number of reconstructed tracks
    int     Nhits , Nhits_stored;   // number of recorded hits in the event
    
    
    
    // my objects
    std::vector<PandoraNuTrack>     tracks;
    std::vector<hit>                hits;
    std::vector<GENIEinteraction>   genie_interactions;
    
    
    //Module labels to get data products
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fFlashModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fGenieGenModuleLabel;

    //mctruth information
    Int_t     mcevts_truth;    //number of neutrino Int_teractions in the spill
    Int_t     nuPDG_truth;     //neutrino PDG code
    Int_t     ccnc_truth;      //0=CC 1=NC
    Int_t     mode_truth;      //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
    Float_t  enu_truth;       //true neutrino energy
    Float_t  Q2_truth;        //Momentum transfer squared
    Float_t  W_truth;         //hadronic invariant mass
    Float_t  X_truth;
    Float_t  Y_truth;
    Int_t    hitnuc_truth;    //hit nucleon
    Int_t    target_truth;    //hit nucleus
    Float_t  nuvtxx_truth;    //neutrino vertex x
    Float_t  nuvtxy_truth;    //neutrino vertex y
    Float_t  nuvtxz_truth;    //neutrino vertex z
    Float_t  nu_dcosx_truth;  //neutrino dcos x
    Float_t  nu_dcosy_truth;  //neutrino dcos y
    Float_t  nu_dcosz_truth;  //neutrino dcos z
    Float_t  lep_mom_truth;   //lepton momentum
    Float_t  lep_dcosx_truth; //lepton dcos x
    Float_t  lep_dcosy_truth; //lepton dcos y
    Float_t  lep_dcosz_truth; //lepton dcos z
    Float_t  t0_truth;        // t0

    
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
                track.SetTruthEndPos( TVector3(particle->EndX() , particle->Endy() , particle->Endz()) );
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
                            genie_interaction.AddPrimary( part.PdgCode()    // pdg code
                                                         ,part.Momentum()   // 4-momentum
                                                         ,part.StatusCode() // status code
                                                         ,part.Mother()     // mother
                                                         );
                        } // for particle
                    }
                    genie_interaction.SortNucleons();
                    genie_interaction.ComputePmissPrec();
                    genie_interaction.FindCC1p200MeVc0pi();
                    
                    // loop over all tracks, and find if any of the tracks
                    // belongs to one of the GENIE particles
                    for (auto track : tracks){
                        if ( track.GetTruthStartPos() == genie_interaction.GetVertexPosition() ){
                            // match GENIE muon to muon track
                            if (track.GetMCpdgCode()==13){
                                if ( track. )
                            }
                            // match GENIE protons to proton tracks
                        }
                    }
                    
                    
                    genie_interactions.push_back( genie_interaction );
                    
                    //                float mindist2 = 9999; // cm;
                    //                TVector3 nuvtx(nuvtxx_truth, nuvtxy_truth, nuvtxz_truth);
                    //                infidvol = insideFidVol(nuvtxx_truth, nuvtxy_truth, nuvtxz_truth);
                    //                //find the closest reco vertex to the neutrino mc truth
                    //                if (infidvol)
                    //                {
                    //                    // vertex is when at least two tracks meet
                    //                    for(size_t i = 0; i < vtxlist.size(); ++i){ // loop over vertices
                    //                        Double_t xyz[3] = {};
                    //                        vtxlist[i]->XYZ(xyz);
                    //                        TVector3 vtxreco(xyz);
                    //                        float dist2 = pma::Dist2(vtxreco, nuvtx);
                    //                        if (dist2 < mindist2)
                    //                        {
                    //                            mindist2 = dist2;
                    //                            vtxrecomc = std::sqrt(dist2);
                    //                            vtxrecomcx = vtxreco.X() - nuvtxx_truth;
                    //                            vtxrecomcy = vtxreco.Y() - nuvtxy_truth;
                    //                            vtxrecomcz = vtxreco.Z() - nuvtxz_truth;
                    //                        }
                    //                    }
                    //
                    //                    // two endpoints of tracks are somehow also vertices...
                    //                    for (size_t i = 0; i < tracklist.size(); ++i){ // loop over tracks
                    //                        float dist2 = pma::Dist2(tracklist[i]->Vertex(), nuvtx);
                    //                        if (dist2 < mindist2)
                    //                        {
                    //                            mindist2 = dist2;
                    //                            vtxrecomc = std::sqrt(dist2);
                    //                            vtxrecomcx = tracklist[i]->Vertex().X() - nuvtxx_truth;
                    //                            vtxrecomcy = tracklist[i]->Vertex().Y() - nuvtxy_truth;
                    //                            vtxrecomcz = tracklist[i]->Vertex().Z() - nuvtxz_truth;
                    //                            
                    //                        }
                    //                        dist2 = pma::Dist2(tracklist[i]->End(), nuvtx);
                    //                        if (dist2 < mindist2)
                    //                        {
                    //                            mindist2 = dist2;
                    //                            vtxrecomc = std::sqrt(dist2);
                    //                            vtxrecomcx = tracklist[i]->End().X() - nuvtxx_truth;
                    //                            vtxrecomcy = tracklist[i]->End().Y() - nuvtxy_truth;
                    //                            vtxrecomcz = tracklist[i]->End().Z() - nuvtxz_truth;
                    //                            
                    //                        }
                    //                    }
                    //                }
                }//is neutrino
            }
        }
        
//        //save g4 particle information
//        std::vector<const simb::MCParticle* > geant_part;
//        
//        // ### Looping over all the Geant4 particles from the BackTracker ###
//        for(size_t p = 0; p < plist.size(); ++p)
//        {
//            // ### Filling the vector with MC Particles ###
//            geant_part.push_back(plist.Particle(p));
//        }
//        
//        //std::cout<<"No of geant part= "<<geant_part.size()<<std::endl;
//        
//        // ### Setting a string for primary ###
//        std::string pri("primary");
//        
//        int primary=0;
//        int geant_particle=0;
//        
//        // ############################################################
//        // ### Determine the number of primary particles from geant ###
//        // ############################################################
//        for( unsigned int i = 0; i < geant_part.size(); ++i ){
//            geant_particle++;
//            // ### Counting the number of primary particles ###
//            if(geant_part[i]->Process()==pri)
//            { primary++;}
//        }//<---End i loop
//        
//        
//        // ### Saving the number of primary particles ###
//        no_primaries=primary;
//        // ### Saving the number of Geant4 particles ###
//        geant_list_size=geant_particle;
//        
//        // ### Looping over all the Geant4 particles ###
//        for( unsigned int i = 0; i < geant_part.size(); ++i ){
//            
//            // ### If this particle is primary, set = 1 ###
//            if(geant_part[i]->Process()==pri)
//            {process_primary[i]=1;}
//            // ### If this particle is not-primary, set = 0 ###
//            else
//            {process_primary[i]=0;}
//            
//            // ### Saving the particles mother TrackID ###
//            Mother[i]=geant_part[i]->Mother();
//            // ### Saving the particles TrackID ###
//            TrackId[i]=geant_part[i]->TrackId();
//            // ### Saving the PDG Code ###
//            pdg[i]=geant_part[i]->PdgCode();
//            // ### Saving the particles Energy ###
//            Eng[i]=geant_part[i]->E();
//            
//            // ### Saving the Px, Py, Pz info ###
//            Px[i]=geant_part[i]->Px();
//            Py[i]=geant_part[i]->Py();
//            Pz[i]=geant_part[i]->Pz();
//            
//            // ### Saving the Start and End Point for this particle ###
//            StartPointx[i]=geant_part[i]->Vx();
//            StartPointy[i]=geant_part[i]->Vy();
//            StartPointz[i]=geant_part[i]->Vz();
//            EndPointx[i]=geant_part[i]->EndPosition()[0];
//            EndPointy[i]=geant_part[i]->EndPosition()[1];
//            EndPointz[i]=geant_part[i]->EndPosition()[2];
//            
//            // ### Saving the processes for this particle ###
//            //std::cout<<"finding proc"<<std::endl;
//            G4Process.push_back( geant_part[i]->Process() );
//            G4FinalProcess.push_back( geant_part[i]->EndProcess() );
//            //std::cout<<"found proc"<<std::endl;
//            //      std::cout << "ID " << TrackId[i] << ", pdg " << pdg[i] << ", Start X,Y,Z " << StartPointx[i] << ", " << StartPointy[i] << ", " << StartPointz[i]
//            //		<< ", End XYZ " << EndPointx[i] << ", " << EndPointy[i] << ", " << EndPointz[i] << ", Start Proc " << G4Process[i] << ", End Proc " << G4FinalProcess[i]
//            //		<< std::endl;
//            
//            // ### Saving the Start direction cosines for this particle ###
//            Startdcosx[i] = geant_part[i]->Momentum(0).Px() / geant_part[i]->Momentum(0).P();
//            Startdcosy[i] = geant_part[i]->Momentum(0).Py() / geant_part[i]->Momentum(0).P();
//            Startdcosz[i] = geant_part[i]->Momentum(0).Pz() / geant_part[i]->Momentum(0).P();
//            // ### Saving the number of Daughters for this particle ###
//            NumberDaughters[i]=geant_part[i]->NumberDaughters();
//            
//        } //geant particles
//        
        
    }//is neutrino

    
    PrintInformation();
    fTree -> Fill();
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::PrintInformation(){
    
    SHOW( fTree->GetEntries() );
    SHOW3( run , subrun , event );
    
    if(!genie_interactions.empty()){
        cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << genie_interactions.size() << " genie interactions\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
        for (auto g: genie_interactions) {
            g.Print();
        }
    }
    
    if(!tracks.empty()){
        cout << "\033[33m" << "xxxxxxxxxxxxxx\n\n" << tracks.size() << " pandoraNu tracks\n\n" << "xxxxxxxxxxxxxx"<< "\033[31m" << endl;
        for (auto t: tracks) {
            t.Print( true );
        }
    }
    
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
    
    // my objects
    fTree->Branch("tracks"              ,&tracks);
    fTree->Branch("hits"                ,&hits);
    fTree->Branch("genie_interactions"  ,&genie_interactions);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::reconfigure(fhicl::ParameterSet const & p){
    fTrackModuleLabel = p.get< std::string >("TrackModuleLabel");
    fHitsModuleLabel  = p.get< std::string >("HitsModuleLabel");
    fGenieGenModuleLabel =  p.get< std::string >("GenieGenModuleLabel");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ub::ErezCCQEAnalyzer::ResetVars(){
    run = subrun = event = -9999;
    Ntracks = Nhits = Nhits_stored = 0;
    tracks.clear();
    hits.clear();
    genie_interactions.clear();
}



// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
DEFINE_ART_MODULE(ub::ErezCCQEAnalyzer)
// - -- - -- - - --- -- - - --- -- - -- - -- -- -- -- - ---- -- - -- -- -- -- -
#ifndef GENIEinteraction_CXX
#define GENIEinteraction_CXX

#include "GENIEinteraction.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GENIEinteraction::GENIEinteraction (Int_t frun, Int_t fsubrun, Int_t fevent, Int_t fmcevent_id):
run(frun),
subrun(fsubrun),
event(fevent),
mcevent_id(fmcevent_id)
{
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::SetKinematics (Float_t fQ2, Float_t fXb, Float_t fW, Float_t fy, Int_t fccnc, Int_t fmode){
    Q2 = fQ2;
    Xb = fXb;
    W = fW;
    y = fy;
    ccnc = fccnc;
    mode = fmode;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::SetMomentumTransfer (){
    // July-2017
    // should only be called after GENIEinteraction::SetNuMomentum() and GENIEinteraction::SetLeptonMomentum()
    q = nu - muon;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GENIEinteraction::AddPrimary ( Int_t fpdg                  // pdg code
                                   ,TLorentzVector fmomentum    // 4-momentum
                                   ,Int_t fstatus_code          // status code
                                   ,Int_t fmother               // mother
){
    
    momentum = fmomentum;
     // interaction neutrino
    switch (fpdg) {
            
        case 12: // ν(e)
            Nnu ++;
            Nnu_e ++;
            break;
            
        case 14: // ν(µ)
            Nnu ++;
            Nnu_mu ++;
            break;
            
        default:
            break;
    }
    
    if (fstatus_code==1) { // status code 0 particles are unstable or do not exit the nucleus and are thus irrelevant
        
        pdg.push_back(fpdg);
        mother.push_back(fmother);
        status_code.push_back(fstatus_code);
        
        Nprimaries++;
        
        switch (pdg.back()) {
                
            case 13: // µ-
                Nmu++;
                Nmu_minus++;
                break;
                
            case -13: // µ+
                Nmu++;
                Nmu_plus++;
                break;
                
            case 2212: // p
                p3vect.push_back( momentum.Vect() ) ;
                Np++;
                break;
                
                
            case 2112: // n
                n3vect.push_back( momentum.Vect() ) ;
                Nn++;
                break;
                
            case 211: // π+
                Npi++;
                Npi_plus++;
                break;
                
            case -211: // π-
                Npi++;
                Npi_minus++;
                break;
                
            case 111: // π0
                Npi++;
                Npi_0++;
                break;
                
            case 11: // e-
                Ne_minus++;
                Nel++;
                break;
                
            case -11: // e+
                Ne_plus++;
                Nel++;
                break;
                
            case 22: // photon
                Ngamma++;
                break;
                
            default:
                break;
        }
    }
    return true;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::AddTrack (PandoraNuTrack ftrack){
    
    tracks.push_back(ftrack);
}


//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//bool GENIEinteraction::AddPrimary ( // GENIE information is for outside of the nucleus (after the reaction)
//                                   Int_t fpdg
//                                   ,Float_t fEng
//                                   ,Float_t fPx
//                                   ,Float_t fPy
//                                   ,Float_t fPz
//                                   ,Float_t fP
//                                   ,Int_t fstatus_code
//                                   ,Float_t fmass
//                                   ,Int_t fND
//                                   ,Int_t fmother
//                                   ,bool track_reconstructed
//                                   ,Int_t ftrackID
//                                   ,PandoraNuTrack fprimarPandoraNuTrack
//                                   ){
//    
//    momentum.SetVectM( TVector3 ( fPx , fPy , fPz ) , fmass );
//    
//    
//    switch (fpdg) { // interaction neutrino
//            
//        case 12: // ν(e)
//            nu = momentum;
//            Nnu ++;
//            Nnu_e ++;
//            break;
//            
//        case 14: // ν(µ)
//            nu = momentum;
//            Nnu ++;
//            Nnu_mu ++;
//            break;
//            
//        default:
//            break;
//    }
//    
//    
//    if (fstatus_code==1) { // status code 0 particles are unstable or do not exit the nucleus and are thus irrelevant
//        
//        pdg.push_back(fpdg);
//        Eng.push_back(fEng);
//        ND.push_back(fND);
//        mother.push_back(fmother);
//        status_code.push_back(fstatus_code);
//        
//        Nprimaries++;
//        if (track_reconstructed){
//            tracks.push_back(fprimarPandoraNuTrack);
//            trackID.push_back(ftrackID);
//        }
//        
//        switch (pdg.back()) {
//                
//            case 13: // µ-
//                muon = momentum;
//                muonTrack = fprimarPandoraNuTrack;
//                muonTrackReconstructed = track_reconstructed;
//                Nmu++;
//                Nmu_minus++;
//                break;
//                
//            case -13: // µ+
//                Nmu++;
//                Nmu_plus++;
//                break;
//                
//            case 2212: // p
//                p3vect.push_back( momentum.Vect() ) ;
//                protonTracks.push_back(fprimarPandoraNuTrack);
//                protonTrack = fprimarPandoraNuTrack; // for CC1p, there is only one proton track, protonTrack.
//                if (protonTrackReconstructed==false) protonTrackReconstructed = track_reconstructed;
//                Np++;
//                break;
//                
//                
//            case 2112: // n
//                n3vect.push_back( momentum.Vect() ) ;
//                Nn++;
//                break;
//                
//            case 211: // π+
//                Npi++;
//                Npi_plus++;
//                break;
//                
//            case -211: // π-
//                Npi++;
//                Npi_minus++;
//                break;
//
//            case 111: // π0
//                Npi++;
//                Npi_0++;
//                break;
//                
//            case 11: // e-
//                Ne_minus++;
//                Nel++;
//                break;
//                
//            case -11: // e+
//                Ne_plus++;
//                Nel++;
//                break;
//
//            case 22: // photon
//                Ngamma++;
//                break;
//                
//            default:
//                break;
//        }
//    }
//    return true;
//}
//
//
//
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GENIEinteraction::FindCC1p200MeVc0pi(){
    // tag the vertex as
    // IsCC_1p_200MeVc_0pi
    // which means that the final state includes
    // 1 protons with momentum >= 200 MeV/c
    // 0 neutrons with momentum >= 200 MeV/c
    // 0 pions
    // 0 electrons
    // 0 photons
    // flagged by GENIE as CC interaction
    Int_t Np_200MeVc = 0 , Nn_200MeVc = 0;
    for (auto Pp : p3vect){
        if (Pp.Mag()>=0.2) Np_200MeVc++;
    }
    for (auto Pn : n3vect){
        if (Pn.Mag()>=0.2) Nn_200MeVc++;
    }
    
    if ( ccnc==0 && Nmu==1 && Np_200MeVc==1 && Npi==0 && Nn_200MeVc==0 && Nel==0 && Ngamma==0 ){
        IsCC_1p_200MeVc_0pi = true;
        for (auto & t : tracks){
            t.SetCC_1p_200MeVc_0pi();
        }
        return true;
    }
    IsCC_1p_200MeVc_0pi = false;
    return false;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GENIEinteraction::SortNucleons(){
    

    for (auto i: sort_by_momentum_magnitude( p3vect )){
        protons.push_back( TLorentzVector( p3vect.at(i) , sqrt( p3vect.at(i).Mag2() + 0.938*0.938 ) ) );
    }
    
    for (auto i: sort_by_momentum_magnitude( n3vect )){
        neutrons.push_back( TLorentzVector( n3vect.at(i) , sqrt( n3vect.at(i).Mag2() + 0.939*0.939 ) ) );
    }
    return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OåOooo........oooOO0OOooo......
vector<size_t> GENIEinteraction::sort_by_momentum_magnitude(const vector<TVector3> &v) {
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1].Mag() > v[i2].Mag();});
    return idx;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GENIEinteraction::ComputePmissPrec(){
    
    //    SHOW(Np);
    //    SHOW(protons.size());
    if(Np>0){
        
        Plead = protons.at(0);
        
        // The leading proton was a neutron before the reaction!
        n_miss = Plead - q;
        
        for (auto proton: protons) {
            Pcm += proton;
        }
        Pcm  -= q;
        Prec  = Pcm - n_miss;
     
        
        // for SRC
        theta_pq    = (Plead.P()>0 && q.P()>0) ? TMath::RadToDeg() * Plead.Vect().Angle(q.Vect()) : -9999;
        p_over_q    = Plead.P()/q.P();
        
        // A(ν,𝝁p) missing mass M²(miss) = (q + (Mp+Mn) - Plead)² , all 4-vectors
        Mmiss       = (q + TLorentzVector( TVector3() , 0.938+0.939 ) - Plead).Mag();

    }
    return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::SetVertexPosition (TVector3 fpos){
    vertex_position = fpos;
    CheckContainement ();
//
//    // load GeometryHelper utility
//    auto geomHelper = ::larutil::GeometryHelper::GetME();
//    double pos_xyz[3] = {fpos.x() , fpos.y() , fpos.z() };
//    // shift in time-axis due to the truncation of the waveforms
//    // (the first 2400 ADCs are removed from the waveform, The extra couple ticks could be due to a shift in the signal deconvolution)
//    double time_shift =  802;
//    
//    
//    for (int plane = 0; plane < 3; plane++){
//        
//        // geoHelper is a set of utility functions that help with geometric stuff..
//        auto const& projection2D = geomHelper->Point_3Dto2D(pos_xyz, plane);
//        int wire = (int) ( projection2D.w / geomHelper->WireToCm() );
//        int time = (int) ( projection2D.t / geomHelper->TimeToCm() ) + time_shift;
//        
//        switch (plane) {
//            case 0:
//                pos_wire_u = wire;
//                pos_time_u = time;
//                break;
//            case 1:
//                pos_wire_v = wire;
//                pos_time_v = time;
//                break;
//            case 2:
//                pos_wire_y = wire;
//                pos_time_y = time;
//                break;
//                
//            default:
//                break;
//        }
//    }
//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GENIEinteraction::Print(bool DoPrintTracks) const{

    cout << "\033[31m" << "GENIE interaction " << mcevent_id << "\n~~~~~~~~~~~~~~~~~~~~~ "<< "\033[0m" << endl;
    SHOW(mcevent_id);
    SHOWTVector3(vertex_position);
    if (!IsVertexContained) Printf("vertex not contained!");
    SHOWstdVector( pdg );
    SHOWTLorentzVector( nu );
    SHOWTLorentzVector( muon );
    SHOWTLorentzVector( q );
    if(!protons.empty()){
        cout << "\033[35m" << protons.size() << " protons:" << endl;
        for (auto proton: protons) {
            SHOWTLorentzVector( proton );
        }
        cout << "\033[31m";
    }
    if(!neutrons.empty()){
        cout << "\033[35m" << neutrons.size() << " neutrons:" << endl;
        for (auto neutron: neutrons) {
            SHOWTLorentzVector( neutron );
        }
        cout << "\033[31m";
    }
    SHOW2( Xb , Q2 );
    SHOW3( theta_pq , p_over_q , Mmiss );
    
    
    SHOW( Nprimaries );
    SHOW3( Np , Nn , Npi );
    SHOW3( Nmu , Nel , Ngamma );
    SHOW2( ccnc , mode );
    SHOW( IsVertexContained )
    
    SHOW( IsCC_1p_200MeVc_0pi );
    if ( IsCC_1p_200MeVc_0pi ) {
        SHOW2( muonTrackReconstructed, protonTrackReconstructed );
        SHOW2( muonTrack.IsTrackContainedSoft() , protonTrack.IsTrackContainedSoft() );
    }

    if(DoPrintTracks && !tracks.empty()){
        cout << "\033[33m" << "--------------\n"
        << tracks.size() << " pandoraNu tracks in GENIE interaction "
        << mcevent_id << "\n--------------" <<  "\033[37m" << endl;
        SHOWstdVector( trackID );
        for (auto t: tracks) {
            t.Print(true);
        }
    }


}


#endif

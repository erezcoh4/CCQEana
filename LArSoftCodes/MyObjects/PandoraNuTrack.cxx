#ifndef PANDORANUTRACK_CXX
#define PANDORANUTRACK_CXX

#include "PandoraNuTrack.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PandoraNuTrack::PandoraNuTrack( Int_t frun, Int_t fsubrun, Int_t fevent
                               ,Int_t ftrack_id
                               ,Float_t flength
                               ,Float_t ftheta, Float_t fphi
                               ,TVector3 fstart_pos, TVector3 fend_pos):
run(frun),
subrun(fsubrun),
event(fevent),
track_id(ftrack_id),
length(flength),
theta(ftheta),
phi(fphi),
start_pos(fstart_pos),
end_pos(fend_pos)
{

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PandoraNuTrack::Print( bool DoPrintPandoraNuFeatures ) const{
    
    cout << "\033[31m" << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl
    << "track " << track_id << endl << "-------------------"    << "\033[0m" << endl;
    
    SHOWTVector3(start_pos);
    SHOWTVector3(end_pos);
    
    if (DoPrintPandoraNuFeatures){
        
        
        SHOW( PIDa ); // SHOW3( PIDaPerPlane[0] , PIDaPerPlane[1] , PIDaPerPlane[2] );
        PrintPhys(length,"cm");
        PrintPhys(theta,"rad");
        PrintPhys(phi,"rad");
        
    }
    if (MCpdgCode!=-9999){
        cout << "........................" << endl << "MC information " << endl ;
        SHOW(MCpdgCode);
        SHOW(mcevent_id);
        SHOWTVector3(truth_start_pos);
        SHOWTVector3(truth_end_pos);
        SHOWTLorentzVector(truth_momentum);
        SHOW(truth_mother);
        SHOW(truth_process);
        SHOW(truth_origin);
        cout << "........................" << endl;
    }
    cout << "closest-flash:" << endl;
    ClosestFlash.Print();
    //        SHOW2( truth_ccnc, IsGENIECC1p );
    //        SHOW( IsGENIECC_1p_200MeVc_0pi );
    //    }
    cout << "\033[31m" << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << "\033[0m" << endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PandoraNuTrack::SetStartEndPlane(Int_t plane ,
                                      Int_t start_wire, Int_t start_time ,
                                      Int_t end_wire, Int_t end_time ){
    
    roi[plane] = box( start_wire , start_time , end_wire , end_time );
    switch (plane) {
        case 0:
            start_wire_u = start_wire;
            start_time_u = start_time;
            end_wire_u = end_wire;
            end_time_u = end_time;
            break;
        case 1:
            start_wire_v = start_wire;
            start_time_v = start_time;
            end_wire_v = end_wire;
            end_time_v = end_time;
            break;
        case 2:
            start_wire_y = start_wire;
            start_time_y = start_time;
            end_wire_y = end_wire;
            end_time_y = end_time;
            break;
            
        default:
            break;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetStartWire (int plane) const{
    switch (plane) {
        case 0:
            return start_wire_u;
            break;
        case 1:
            return start_wire_v;
            break;
        case 2:
        default:
            return start_wire_y;
            break;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetStartTime (int plane) const{
    switch (plane) {
        case 0:
            return start_time_u;
            break;
        case 1:
            return start_time_v;
            break;
        case 2:
        default:
            return start_time_y;
            break;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetEndWire (int plane) const{
    switch (plane) {
        case 0:
            return end_wire_u;
            break;
        case 1:
            return end_wire_v;
            break;
        case 2:
        default:
            return end_wire_y;
            break;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetEndTime (int plane) const{
    switch (plane) {
        case 0:
            return end_time_u;
            break;
        case 1:
            return end_time_v;
            break;
        case 2:
        default:
            return end_time_y;
            break;
    }
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::DistanceFromPoint( TVector3 position, std::string * fStartOrEnd  ){
    Float_t DistanceStart, DistanceEnd , distance = 1000;
    std::string StartOrEnd = "None";
    
    DistanceStart = ( start_pos - position).Mag();
    DistanceEnd = ( end_pos - position).Mag();
    if ( DistanceStart < DistanceEnd ){
        StartOrEnd = "Start";
        distance = DistanceStart;
    }
    else{
        StartOrEnd = "End";
        distance = DistanceEnd;
    }
    
    return distance;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::ClosestDistanceToOtherTrack( PandoraNuTrack other_track, std::string * fStartOrEnd ){
    Float_t MinDistanceToOtherTrack = 10000;
    std::string StartOrEnd = "None";
    Float_t DistanceStartStart = (start_pos - other_track.start_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceStartStart)     {MinDistanceToOtherTrack = DistanceStartStart; StartOrEnd = "Start";}
    
    Float_t DistanceStartEnd = (start_pos - other_track.end_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceStartEnd)       {MinDistanceToOtherTrack = DistanceStartEnd; StartOrEnd = "Start";}
    
    Float_t DistanceEndStart = (end_pos - other_track.start_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceEndStart)       {MinDistanceToOtherTrack = DistanceEndStart; StartOrEnd = "End";}
    
    Float_t DistanceEndEnd = (end_pos - other_track.end_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceEndEnd)         {MinDistanceToOtherTrack = DistanceEndEnd; StartOrEnd = "End";}
    
    
    if (fStartOrEnd!=nullptr) *fStartOrEnd = StartOrEnd;
    
    return MinDistanceToOtherTrack;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PandoraNuTrack::FlipTrack(){
    
    // flip start and end positions
    TVector3 tmp_pos = start_pos;
    start_pos   = end_pos;
    end_pos     = tmp_pos;
    
    // change angles
    theta       = 3.1416 - theta;
    phi         = (phi > 0 ? phi-3.1416 : phi+3.1416);
    
    //    Float_t tmp_dqdx = start_dqdx;
    //    start_dqdx  = end_dqdx;
    //    end_dqdx    = tmp_dqdx;
    
    //    for (int plane = 0 ; plane < 3 ; plane++ ){
    //
    //        Float_t tmp1 = dqdx_around_start[plane];
    //        dqdx_around_start[plane] = dqdx_around_end[plane];
    //        dqdx_around_end[plane] = tmp1;
    //
    //        Float_t tmp2 = dqdx_around_start_track_associated[plane];
    //        dqdx_around_start_track_associated[plane] = dqdx_around_end_track_associated[plane];
    //        dqdx_around_end_track_associated[plane] = tmp2;
    //    }
    //
    //    Float_t     tmp3 = dqdx_around_start_total;
    //    dqdx_around_start_total = dqdx_around_end_total;
    //    dqdx_around_end_total = tmp3;
    //
    //    Float_t     tmp4 = dqdx_around_start_track_associated_total;
    //    dqdx_around_start_track_associated_total = dqdx_around_end_track_associated_total;
    //    dqdx_around_end_track_associated_total = tmp4;
}


#endif

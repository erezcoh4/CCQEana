#ifndef hit_CXX
#define hit_CXX

#include "hit.h"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
hit::hit(Int_t fplane, Int_t fwire, Int_t fid, Float_t fpeakT, Float_t fcharge):
hit_plane(fplane),
hit_wire(fwire),
hit_id(fid),
hit_peakT(fpeakT),
hit_charge(fcharge)
{}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void hit::Print() const{
    std::cout << "hit " << hit_id << ", " << hit_plane << "/" << hit_wire << "/" << hit_peakT << std::endl;
    std::cout << "\t trkKey: " << hit_trkKey << std::endl;

}


#endif

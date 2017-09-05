/**
 * \file PandoraNuTrack.h
 *
 * \ingroup GBDTprotonsPackage
 *
 * \brief Class def header for a class PandoraNuTrack
 *
 * @author erezcohen
 */

/** \addtogroup GBDTprotonsPackage
 
 @{*/
#ifndef PANDORANUTRACK_H
#define PANDORANUTRACK_H

#include "Rtypes.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "box.h"
#include "flash.h"
#include <iostream>
#include <iomanip>

using namespace std;

// prints....
#define EndEventBlock() std::cout << "\033[32m"<< "....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......" << "\033[0m"<< endl
#define PrintLine() std::cout << "--------------------------------------------------------------" << std::endl
#define PrintXLine() std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl
#define SHOW(a) std::cout << setprecision(2) << fixed << #a << ": " << (a) << std::endl;
#define SHOW2(a,b) std::cout <<"\033[34m"<<#a<<": "<<(a)<<"," << #b <<": "<<(b)<< "\033[0m"<< std::endl;
#define SHOW3(a,b,c) std::cout <<"\033[36m"<<#a<<": "<<(a)<<"," << #b <<": "<<(b)<<","<<#c<<": "<<(c)<< "\033[0m"<< std::endl;
#define SHOW4(a,b,c,d) std::cout <<"\033[31m"<<#a<<": "<<(a)<<"," << #b <<": "<<(b)<<","<<#c<<": "<<(c)<<","<<#d<<": "<<(d)<< "\033[0m"<< std::endl;
#define SHOWstdVector(v){ if (v.size()<1) {std::cout << #v << " is empty" << std::endl;} else {std::cout << #v << "( " << v.size() << " entries):\t"; for (auto it:v) std::cout << it << "\t"; std::cout << endl;}}
#define SHOWTVector3(v){ std::cout << #v << ": (" << v.X() << "," << v.Y() << "," << v.Z() << ")" << endl;}
#define SHOWTLorentzVector(v) std::cout << #v << ": " << "\t(" << setprecision(2) << fixed << v.Px() << ","  << v.Py() << "," << v.Pz()  << "," << v.E() << ")" << ", P = " << v.P() << ", M = " << v.M() << std::endl

#define PrintPhys(a,units) std::cout  << setprecision(2) << fixed << #a << ": " << (a) <<  " " << (units) << std::endl



/**
 \class PandoraNuTrack
 User defined class PandoraNuTrack ... these comments are used to generate
 doxygen documentation!
 */





class PandoraNuTrack{
    
public:
    
    /// Default constructor
    PandoraNuTrack () = default;
    PandoraNuTrack (Int_t frun, Int_t fsubrun, Int_t fevent
                    ,Int_t ftrack_id
                    ,Float_t flength
                    ,Float_t ftheta, Float_t fphi
                    ,TVector3 fstart_pos, TVector3 fend_pos );
    
    
    
    /// SETters
    void      SetTruthProcess (std::string _process)        {truth_process = _process;};
    void            SetOrigin (std::string _origin)         {truth_origin = _origin;};

    void SetCC_1p_200MeVc_0pi ()                            {IsGENIECC_1p_200MeVc_0pi=true;};
    
    void               SetRun (Int_t _run)                  {run = _run;};
    void            SetSubrun (Int_t _subrun)               {subrun = _subrun;};
    void             SetEvent (Int_t _event)                {event = _event;};
    void         SetMCeventID (Int_t fmcevent_id)           {mcevent_id = fmcevent_id;};
    void              SetCCNC (Int_t fccnc)                 {truth_ccnc = fccnc;};
    void              SetMode (Int_t fmode)                 {truth_mode = fmode;};
    void         SetMCpdgCode (Int_t _mcpdg)                {MCpdgCode = _mcpdg;};
    void         SetBestPlane (Int_t best)                  {BestPlane = best;};
    void          SetMaxNHits (Int_t maxnhits)              {MaxNHits = maxnhits;};
    void       SetTruthMother (Int_t _mother)               {truth_mother = _mother;};
    void    SetCaloKEPerPlane (Int_t plane , Float_t ke )   {CaloKEPerPlane[plane] = ke;};
    void      SetPIDaPerPlane (Int_t plane , Float_t pida ) {PIDaPerPlane[plane] = pida;};
    
    
    void            SetLength (Float_t l)                   {length = l;};
    void             SetTheta (Float_t t)                   {theta = t;};
    void               SetPhi (Float_t ph)                  {phi = ph;};
    void       SetTruthLength ()                            {truth_length = (truth_start_pos-truth_end_pos).Mag();};
    void        SetTruthTheta (Float_t f)                   {truth_theta = f;};
    void          SetTruthPhi (Float_t f)                   {truth_phi = f;};

    void          SetStartPos (TVector3 pos)                {start_pos = pos;};
    void            SetEndPos (TVector3 pos)                {end_pos = pos;};
    void     SetTruthStartPos (TVector3 pos)                {truth_start_pos = pos;};
    void       SetTruthEndPos (TVector3 pos)                {truth_end_pos = pos;};
    void              SetPIDa ()                            {PIDa = PIDaPerPlane[BestPlane]; };
    void     SetTruthMomentum (TLorentzVector fmomentum)    {truth_momentum=fmomentum; };
    
    void     SetStartEndPlane (Int_t plane ,
                               Int_t start_wire, Int_t start_time ,
                               Int_t end_wire, Int_t end_time );
    void      SetClosestFlash (flash _flash)                {ClosestFlash = _flash;};
    
    
    
    
    // GETters

    std::string     GetTruthProcess () const {return truth_process;};
    std::string           GetOrigin () const {return truth_origin;};    // "unknown origin" / "beam neutrino" / "cosmic ray"

    
    Int_t                    GetRun () const {return run;};
    Int_t                 GetSubrun () const {return subrun;};
    Int_t                  GetEvent () const {return event;};
    Int_t                GetTrackID () const {return track_id;};
    Int_t              GetMCeventID () const {return mcevent_id;};
    Int_t            GetTruthMother () const {return truth_mother;};
    
    Int_t              GetMCpdgCode () const {return MCpdgCode;};
    Int_t                   GetCCNC () const {return truth_ccnc;};
    Int_t              GetBestPlane () const {return BestPlane;};
    Int_t               GetMaxNHits () const {return MaxNHits;};

    Int_t              GetStartWire (int plane) const;
    Int_t              GetStartTime (int plane) const;
    Int_t                GetEndWire (int plane) const;
    Int_t                GetEndTime (int plane) const;

    
    Float_t               GetLength () const {return length;};
    Float_t                GetTheta () const {return theta;};
    Float_t                  GetPhi () const {return phi;};
    Float_t          GetTruthLength () const {return truth_length;};
    Float_t           GetTruthTheta () const {return truth_theta;};
    Float_t             GetTruthPhi () const {return truth_phi;};
    
    Float_t       GetCaloKEPerPlane ( Int_t plane ) const { return CaloKEPerPlane[plane];};
    Float_t         GetPIDaPerPlane ( Int_t plane ) const { return PIDaPerPlane[plane];};
    Float_t                 GetPIDa () const {return PIDa; };

    TVector3            GetStartPos () const {return start_pos;};
    TVector3              GetEndPos () const {return end_pos;};
    TVector3       GetTruthStartPos () const {return truth_start_pos;};
    TVector3         GetTruthEndPos () const {return truth_end_pos;};
    TLorentzVector GetTruthMomentum () const {return truth_momentum; };
    flash           GetClosestFlash () const {return ClosestFlash;};
    
    
    
    
    
    
    
    
    // functionallity
    void            FlipTrack ();
    void           CreateROIs ();
    void                Print (bool DoPrintPandoraNuFeatures = true ) const;

    
    // operators
    bool              AskIfCC1p0pi () const {return IsGENIECC_1p_200MeVc_0pi;};
    bool IsTrackContainedSoft (float max_FV_y = 115,
                              float min_FV_z = 5, float max_FV_z = 1045,
                              float min_FV_x = 3, float max_FV_x = 257)  const{
        if( ( start_pos.x() < min_FV_x )    | ( start_pos.x() > max_FV_x ) )    return false;
        if( ( start_pos.y() < -max_FV_y )   | ( start_pos.y() > max_FV_y ) )    return false;
        if( ( start_pos.z() < min_FV_z )    | ( start_pos.z() > max_FV_z ) )    return false;
        if( ( end_pos.x() < min_FV_x )      | ( end_pos.x() > max_FV_x ) )      return false;
        if( ( end_pos.y() < -max_FV_y )     | ( end_pos.y() > max_FV_y ) )      return false;
        if( ( end_pos.z() < min_FV_z )      | ( end_pos.z() > max_FV_z ) )      return false;
        return true;
    };
    
    Float_t           DistanceFromPoint ( TVector3 position , std::string * StartOrEnd=nullptr );
    Float_t ClosestDistanceToOtherTrack ( PandoraNuTrack other_track , std::string * StartOrEnd=nullptr );
    
    inline bool operator==(const PandoraNuTrack & t) {
        return  (run==t.GetRun() && subrun==t.GetSubrun() && event==t.GetEvent() && track_id==t.GetTrackID());
    }
    inline bool operator!=(const PandoraNuTrack & t) {
        return  (run!=t.GetRun() || subrun!=t.GetSubrun() || event!=t.GetEvent() || track_id!=t.GetTrackID());
    }

    
private:

    // std::string
    std::string truth_process="unknown process";
    std::string truth_origin="unknown origin";
    
    // bool
    bool        IsGENIECC_1p_200MeVc_0pi=false;
    
    // Int_t
    Int_t       run=0, subrun=0, event=0, track_id=-9999;
    Int_t       truth_ccnc=0, truth_mode=0;
    Int_t       MCpdgCode=-9999;
    Int_t       start_wire_u=0, start_wire_v=0, start_wire_y=0;
    Int_t       start_time_u=0, start_time_v=0, start_time_y=0;
    Int_t       end_wire_u=0, end_wire_v=0, end_wire_y=0;
    Int_t       end_time_u=0, end_time_v=0, end_time_y=0;
    Int_t       BestPlane=2;
    Int_t       MaxNHits=0;
    Int_t       mcevent_id=-1;
    Int_t       truth_mother=-1;
    
    // Float_t
    Float_t     length=0, theta=0, phi=0;
    Float_t     CaloKEPerPlane[3]={0,0,0};
    Float_t     PIDaPerPlane[3]={0,0,0};
    Float_t     PIDa=-1;

    // truth information - only valid for MC data
    Float_t     truth_theta=-9999, truth_phi=-9999, truth_length=-9999;
    
    
    // TVector3
    TVector3    start_pos=TVector3(), end_pos=TVector3();
    TVector3    truth_start_pos=TVector3(), truth_end_pos=TVector3();
    
    TLorentzVector truth_momentum=TLorentzVector();
    
    box         roi[3]={box(),box(),box()};
    
    flash       ClosestFlash=flash();
};
#endif
/** @} */ // end of doxygen group


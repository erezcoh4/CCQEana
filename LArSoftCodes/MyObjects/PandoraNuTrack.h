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

//#include "hit.h"
//#include "box.h"


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
//    void          SetTruthEng (Float_t f)                   {truth_Eng = f;};
//    void            SetTruthP (Float_t f)                   {truth_P = f;};
//    void         SetTruthMass (Float_t f)                   {truth_Mass = f;};
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
    
    
    
    
    
    // GETters

    std::string     GetTruthProcess () const {return truth_process;};

    
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
//    Int_t              GetStartWire (int plane) const {return roi[plane].GetStartWire();};
//    Int_t              GetStartTime (int plane) const {return roi[plane].GetStartTime();};
//    Int_t                GetEndWire (int plane) const {return roi[plane].GetEndWire();};
//    Int_t                GetEndTime (int plane) const {return roi[plane].GetEndTime();};

    Int_t              GetStartWire (int plane) const;
    Int_t              GetStartTime (int plane) const;
    Int_t                GetEndWire (int plane) const;
    Int_t                GetEndTime (int plane) const;

    
    Float_t               GetLength () const {return length;};
    Float_t                GetTheta () const {return theta;};
    Float_t                  GetPhi () const {return phi;};
//    Float_t             GetTruthEng () const {return truth_Eng;};
//    Float_t               GetTruthP () const {return truth_P;};
//    Float_t            GetTruthMass () const {return truth_Mass;};
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
    Float_t     PIDa;
    // truth information - only valid for MC data
    //    Float_t     truth_Eng=-1, truth_P=-1 , truth_Mass=-1;
    Float_t     truth_theta=-9999, truth_phi=-9999, truth_length=-9999;
    
    
    // TVector3
    TVector3    start_pos=TVector3(), end_pos=TVector3();
    TVector3    truth_start_pos=TVector3(), truth_end_pos=TVector3();
    
    TLorentzVector truth_momentum=TLorentzVector();
    
    box         roi[3]={box(),box(),box()};
    
    /* other features of a pandoraNu track that might be usefull in the future
     
     
     
     
     other features of a pandoraNu track that might be usefull in the future */
    
    
    /*
     void      Calorimetry ();
     void     Straightness ();
     void      SetMomentum (Float_t,Float_t);
     
     
     
     // setters
     void     Set_start_dqdx (Float_t dqdx)  {start_dqdx = dqdx;};
     void       Set_end_dqdx (Float_t dqdx)  {end_dqdx = dqdx;};
     void       Set_tot_dqdx (Float_t dqdx)  {tot_dqdx = dqdx;};
     void       Set_avg_dqdx (Float_t dqdx)  {avg_dqdx = dqdx;};
     void          Set_nhits (Int_t n)       {nhits = n;};
     void  SetCalorimetryPDG (Int_t _pdg[3]) {for (int i=0 ; i < 3 ; i++ ) CalorimetryPDG[i] = _pdg[i];};
     void  SetProcessPrimary (Int_t fpp)     {process_primary = fpp;};
     void        SetCF2Start (Float_t fc)    {cfdistance_start = fc;}; // set the closest flash distance to the start point of the track
     
     
     
     void       SetCosScores (Float_t fcscore, Float_t fccscore)
     {cosmicscore = fcscore; coscontscore = fccscore;};
     
     void       Set_pid_info (Float_t fpida, Float_t fchi)
     {pidpida = fpida; pidchi = fchi;};
     
     void     SetTrackPurity (Float_t fpurtruth_U, Float_t fpurtruth_V, Float_t fpurtruth_Y)
     { purtruth_U = fpurtruth_U; purtruth_V = fpurtruth_V; purtruth_Y = fpurtruth_Y ;};
     
     void           Set_dqdx (Float_t, Float_t, Float_t, Int_t);
     void       SetFlashInfo (Float_t fcftime, Float_t fcftimewidth, Float_t fcfzcenter, Float_t fcfzwidth, Float_t fcfycenter, Float_t fcfywidth, Float_t fcftotalpe, Float_t fcfdistance);
     
     
     
     void           Set_dEdx (vector<Float_t> , vector<Float_t> , vector<Float_t> , vector<Float_t> , vector<Float_t> ,
     vector<Float_t> , vector<Float_t> , vector<Float_t> , vector<Float_t> , vector<Float_t> ,
     vector<Float_t> , vector<Float_t> , vector<Float_t> , vector<Float_t> , vector<Float_t>  );
     
     void    SetCalorimetry_Y (vector<Float_t> , vector<Float_t> , vector<Float_t> , vector<Float_t>   );
     
     //    void   SetSWtrigger (std::string * fswtrigger_name, bool * fswtrigger_triggered){
     //        for (size_t i=0 ; i < (int)(sizeof(fswtrigger_triggered)/sizeof(fswtrigger_triggered[0])) ; i++ ) {
     //            swtrigger_name.push_back(fswtrigger_name[i]);
     //            swtrigger_triggered.push_back(fswtrigger_triggered[i]);
     //        }
     //    };
     
     void    SetSlopeIntercept ( int plane = 0 , float fslope = -1000 , float fintercept = -1000 ){
     slope[plane] = fslope;
     intercept[plane] = fintercept;
     };
     void SetX1Y1X2Y2forTrack (int plane, std::vector<float> fx1x2y1y2) {
     for(auto f:fx1x2y1y2) x1y1x2y2[plane].push_back(f);
     if (WireTimeAngle[plane]==-100)
     WireTimeAngle[plane] = atan2(x1y1x2y2[plane][3]-x1y1x2y2[plane][1],
     x1y1x2y2[plane][2]-x1y1x2y2[plane][0]);
     };
     
     
     
     // finders
     bool           IsWireTimeAlongTrack ( Int_t fplane, Int_t fwire , Float_t fPeakTime );
     
     
     
     // getters
     box              GetROI (int plane) {return roi[plane];};
     Int_t        GetCaloPDG (int plane) {return CalorimetryPDG[plane];};
     std::vector<Float_t> GetEdepYInfo (int step) {
     std::vector<Float_t> result = {residual_range_Y.at(step), dqdx_Y.at(step), dEdx_Y.at(step), Edep_Y.at(step)};
     return result;
     };
     std::vector<Float_t> GetTrackLengthVector (int plane) {
     switch (plane) {
     case 0:
     //                return residual_range_U;
     //                break;
     case 1:
     //                return residual_range_V;
     //                break;
     case 2:
     return residual_range_Y;
     break;
     default:
     return residual_range_Y;
     break;
     }};
     std::vector<Float_t> GetTrack_dEdxVector  (int plane) {
     switch (plane) {
     case 0:
     //                return dEdx_U;
     //                break;
     case 1:
     //                return dEdx_V;
     //                break;
     case 2:
     return dEdx_Y;
     break;
     default:
     return dEdx_Y;
     break;
     }};
     //    Int_t GetNSWtrigger () {return (int)swtrigger_name.size();};
     
     std::vector<float> GetX1Y1X2Y2forTrack( int plane = 0 ){
     if(!x1y1x2y2[plane].empty()) return x1y1x2y2[plane];
     else return std::vector<float> {0,0,0,0};
     }
     
     
     // operators
     inline bool operator==(const PandoraNuTrack & t) {
     return std::tie( run, subrun, event, track_id ) == std::tie(t.run, t.subrun, t.event, t.track_id);
     }
     
     
     
     // features
     bool        IsStartContained, IsEndContained, IsFullyContained;
     
     Int_t       CalorimetryPDG[3];
     Int_t       nhits       , is_flipped;
     // Int_t       NNeighborTracks;
     
     Float_t     startx  , starty , startz , endx , endy , endz;
     
     
     Float_t     momrange    , mommsllhd , momeavgrangellhd;
     Float_t     start_dqdx  , end_dqdx  , tot_dqdx , avg_dqdx;
     Float_t     dqdx_diff   , dqdx_ratio;
     Float_t     dQtotal;
     Float_t     pidpida     , pidchi    , cosmicscore   , coscontscore;
     Float_t     cftime      , cftimewidth   , cfzcenter , cfzwidth, cfycenter , cfywidth  , cftotalpe , cfdistance;
     Float_t     cfdistance_start;
     
     // charge deposition around start and end points
     Float_t     dqdx_around_start[3]                       , dqdx_around_end[3];
     Float_t     dqdx_around_start_total                    , dqdx_around_end_total;
     Float_t     dqdx_around_start_track_associated[3]      , dqdx_around_end_track_associated[3];
     Float_t     dqdx_around_start_track_associated_total   , dqdx_around_end_track_associated_total;
     Float_t     trajectory_slope[3], trajectory_intersect[3];
     // The trkpurtruth - purity variable is defined as the ratio of the energy of the particle that contributed most to this track in a given plane to the total energy coming from all particles that contribute to this track in that plane
     Float_t     purtruth_U  , purtruth_V    , purtruth_Y;
     Float_t     slope[3]    , intercept[3]  , WireTimeAngle[3];
     std::vector<float> x1y1x2y2[3];
     
     TString     TopBottDir  , ForBackDir    , LefRghtDir;
     
     // box         start_box[3], end_box[3];   // boxed around the start and end points of the track
     
     // dE/dx
     Int_t       NEdepYsteps;
     // std::vector <Float_t> track_dx_U, residual_range_U, dEdx_U , Edep_U, dqdx_U;
     // std::vector <Float_t> track_dx_V, residual_range_V, dEdx_V , Edep_V, dqdx_V;
     std::vector <Float_t> track_dx_Y, residual_range_Y, dEdx_Y , Edep_Y, dqdx_Y;
     
     
     
     
     
     //    // software trigger
     //    std::vector<std::string> swtrigger_name;       // the name of the trigger algorithm
     //    std::vector<bool>        swtrigger_triggered;  // true = event is triggered; false = event is not triggered based on the relative algorithm logic
     
     
     
     */
    
};
#endif
/** @} */ // end of doxygen group


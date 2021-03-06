#ifndef GENIEFILE_CXX
#define GENIEFILE_CXX

#include "GenieFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GenieFile::GenieFile( TString fPath
                     ,TString fRootFileName     // without the .root suffix (!)
                     ,TString fRootTreeName
                     ,int fdebug
                     ):
debug(fdebug),
Path( fPath ),
RootFileName( fRootFileName + ".root"),
RootTreeName( fRootTreeName ),
OutputCSVname( fRootFileName + ".csv")
{
    SetLArTools();
    SetInTree();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GenieFile::GenieFile( TString fPath
                     ,TString fRootFileName     // without the .root suffix (!)
                     ,TString fRootTreeName
                     ,int fdebug
                     ,TString fAccMapPath
                     ,TString fPmuThetaAccMapName,TString fPpThetaAccMapName,TString fQ2AccMapName):
debug(fdebug),
Path( fPath ),
RootFileName( fRootFileName + ".root"),
RootTreeName( fRootTreeName ),
OutputCSVname( fRootFileName + ".csv")
{
    SetLArTools();
    SetInTree();
    SetPmuThetaAcceptanceMaps(fAccMapPath,fPmuThetaAccMapName);
    SetPpThetaAcceptanceMaps(fAccMapPath,fPpThetaAccMapName);
    SetQ2AcceptanceMaps(fAccMapPath,fQ2AccMapName);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<double> GenieFile::Read1dArrayFromFile(TString filename){
    // read a 1d array from a file
    // which is of unknown size
    // and each line is a different number
    // e.g.
    // 1
    // 2.1
    // 3.05
    // 5.2
    // ....
    Debug(3,"GenieFile::Read1dArrayFromFile(%)",filename);
    ifstream fin;
    fin.open(filename);
    std::vector<double> numbers;
    double number;
    while(fin >> number) //  >> delim
        numbers.push_back(number);
    fin.close();
    if(debug>4){
        SHOWstdVector(numbers);
    }
    return numbers;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<std::vector<double>> GenieFile::Read2dArrayFromFile(TString filename){
    // read a 2d array from a file
    // which is of unknown size
    // and each line is a different number
    // e.g.
    // 1,3,2.2,10
    // 2.1,3.0,8.2,1
    // 3.05,1.1,29.2,3
    // ....
    Debug(3,"GenieFile::Read2dArrayFromFile(%)",filename);
    ifstream fin;
    fin.open(filename);
    std::vector<std::vector<double>> numbers;
    string line;
    while(!fin.eof()){
        getline ( fin, line, '\n' );
        stringstream sep(line);
        string numbers_in_line;
        numbers.push_back(std::vector<double>());
        while (getline(sep, numbers_in_line, ',')) {
            numbers.back().push_back(stod(numbers_in_line));
        }
    }
    fin.close();
    if(debug>5){
        for (auto row : numbers) {
            for (auto number : row) {
                cout << std::setprecision(5) << number << ' ';
            }
            cout << '\n';
        };
    }
    return numbers;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::SetPmuThetaAcceptanceMaps(TString fAccMapPath,TString fAccMapName){
    TString xbinsfilename   = fAccMapPath + fAccMapName + "_xbins.csv";
    TString ybinsfilename   = fAccMapPath + fAccMapName + "_ybins.csv";
    TString accfilename     = fAccMapPath + fAccMapName + "_acceptance.csv";
    TString errfilename     = fAccMapPath + fAccMapName + "_acc_err.csv";
    
    Pmu_xbins               = Read1dArrayFromFile(xbinsfilename);
    Pmu_theta_ybins         = Read1dArrayFromFile(ybinsfilename);
    Pmu_theta_acceptance    = Read2dArrayFromFile(accfilename);
    Pmu_theta_acc_err       = Read2dArrayFromFile(errfilename);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::SetPpThetaAcceptanceMaps(TString fAccMapPath,TString fAccMapName){
    TString xbinsfilename   = fAccMapPath + fAccMapName + "_xbins.csv";
    TString ybinsfilename   = fAccMapPath + fAccMapName + "_ybins.csv";
    TString accfilename     = fAccMapPath + fAccMapName + "_acceptance.csv";
    TString errfilename     = fAccMapPath + fAccMapName + "_acc_err.csv";
    
    Pp_xbins                = Read1dArrayFromFile(xbinsfilename);
    Pp_theta_ybins          = Read1dArrayFromFile(ybinsfilename);
    Pp_theta_acceptance     = Read2dArrayFromFile(accfilename);
    Pp_theta_acc_err        = Read2dArrayFromFile(errfilename);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::SetQ2AcceptanceMaps(TString fAccMapPath,TString fAccMapName){
    TString Q2binsfilename   = fAccMapPath + fAccMapName + "_bins.csv";
    TString accfilename     = fAccMapPath + fAccMapName + "_acceptance.csv";
    TString errfilename     = fAccMapPath + fAccMapName + "_acc_err.csv";
    
    Q2_bins           = Read1dArrayFromFile(Q2binsfilename);
    Q2_acceptance     = Read1dArrayFromFile(accfilename);
    Q2_acc_err        = Read1dArrayFromFile(errfilename);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::SetQ2_gen_rec_map(TString fMapPath,TString fMapName){
    TString Q2_gen_rec_bins_filename   = fMapPath + fMapName + "_bins.csv";
    TString Q2_gen_rec_map_filename    = fMapPath + fMapName + "_map.csv";
    
    Q2_gen_rec_bins    = Read1dArrayFromFile(Q2_gen_rec_bins_filename);
    Q2_gen_rec_map     = Read2dArrayFromFile(Q2_gen_rec_map_filename);
    SetQ2_rec_1d_Probabilities();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::SetQ2_rec_1d_Probabilities(){
    int fDebug = 3;
    Debug(fDebug,"GenieFile::SetQ2_rec_1d_Probabilities()");
    // use the Q2 gen/rec 2D map to create a set of 1D distributions
    // from which we generate "reco"-Q2 for each "gen"-Q2
    for ( size_t i_gen=0; i_gen < Q2_gen_rec_bins.size(); i_gen++ ) {
        // for each i_gen we build 1D distirbution from all i_rec
        std::vector<double> Q2_rec_vector = Q2_gen_rec_map.at(i_gen);
        // std::discrete_distribution seems to have a constructor that takes an iterator range.
        std::discrete_distribution<double> Q2_rec_distribution(Q2_rec_vector.begin(), Q2_rec_vector.end());
        Q2_rec_distributions.push_back( Q2_rec_distribution );
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::SetInTree(){
    InputGenieFile = new TFile( Path + RootFileName );
    GenieTree = (TTree*)InputGenieFile->Get(RootTreeName);
    SHOW( GenieTree->GetEntries() );
    
    
    GenieTree -> SetBranchAddress("neu"            , &neu);
    GenieTree -> SetBranchAddress("cc"             , &cc);
    GenieTree -> SetBranchAddress("qel"            , &qel);
    GenieTree -> SetBranchAddress("res"            , &res);
    GenieTree -> SetBranchAddress("dis"            , &dis);
    GenieTree -> SetBranchAddress("coh"            , &coh);
    GenieTree -> SetBranchAddress("mec"            , &mec);
    
    GenieTree -> SetBranchAddress("Q2"             , &Q2);
    GenieTree -> SetBranchAddress("W"              , &W);
    GenieTree -> SetBranchAddress("x"              , &x);
    // neutrino
    GenieTree -> SetBranchAddress("pxv"            , &pxv);
    GenieTree -> SetBranchAddress("pyv"            , &pyv);
    GenieTree -> SetBranchAddress("pzv"            , &pzv);
    // initial state hit nucleon
    GenieTree -> SetBranchAddress("pxn"            , &pxn);
    GenieTree -> SetBranchAddress("pyn"            , &pyn);
    GenieTree -> SetBranchAddress("pzn"            , &pzn);
    // lepton
    GenieTree -> SetBranchAddress("pxl"            , &pxl);
    GenieTree -> SetBranchAddress("pyl"            , &pyl);
    GenieTree -> SetBranchAddress("pzl"            , &pzl);
    
    GenieTree -> SetBranchAddress("nf"             , &nf);
    GenieTree -> SetBranchAddress("pxf"            , &pxf);
    GenieTree -> SetBranchAddress("pyf"            , &pyf);
    GenieTree -> SetBranchAddress("pzf"            , &pzf);
    GenieTree -> SetBranchAddress("pdgf"           , &pdgf);

    // kth particle in ‘primary’ hadronic system (in GeV).
    GenieTree -> SetBranchAddress("ni"             , &ni);
    GenieTree -> SetBranchAddress("pxi"            , &pxi);
    GenieTree -> SetBranchAddress("pyi"            , &pyi);
    GenieTree -> SetBranchAddress("pzi"            , &pzi);
    GenieTree -> SetBranchAddress("pdgi"           , &pdgi);
    
    GenieTree -> SetBranchAddress("vtxx"           , &vtxx);
    GenieTree -> SetBranchAddress("vtxy"           , &vtxy);
    GenieTree -> SetBranchAddress("vtxz"           , &vtxz);

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::Initialize(){
    
    for (int i=0; i<NMAX ; i++ ){
        pdgf[i] = -9999;
        pxf[i] = pyf[i] = pzf[i] = -9999;
    }
    IsCC_Np_200MeVc = IsCC_1p_200MeVc = IsCC_1p_200MeVc_0pi=false;
    Pmiss = muon  = proton = q = nu = proton_before_FSI = TLorentzVector();
    protons.clear();
    neutrons.clear();
    reco_l_mu = reco_l_p = 0;
    reco_Q2 = reco_Ev = 0;
    reco_Pnu = reco_Pp = reco_Pmu = reco_q = TLorentzVector();
    i_rec = i_gen = gen_Q2_gen_rec = rec_Q2_gen_rec = 0;
    muonTrajectory_end = muonTrack_end = TVector3();
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::ReadEvent (int fi_event){
    i_event = fi_event;
    Initialize();
    GenieTree->GetEntry(i_event);
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::HeaderCSV (){
    csv_file.open( Path + OutputCSVname );
    
    csv_file
    << "cc"  << "," << "qel" <<  "," << "res" <<  "," << "dis" <<  "," << "coh" <<  ","
    << "mec" << ","
    << "truth_Q2"  << "," << "W" <<  "," << "x" <<  ","
    << "pxn" << "," << "pyn" <<  "," << "pzn" <<  "," // Initial state hit nucleon
    << "pxv" << "," << "pyv" <<  "," << "pzv" <<  ","
    << "pxl" << "," << "pyl" <<  "," << "pzl" << ","
    << "nf"  << ",";
    
    for (int i=0; i<NMAX ; i++ ){
        csv_file
        << Form("pdgf_%d",i) <<  ","
        << Form("pxf_%d",i) <<  "," << Form("pyf_%d",i) <<  "," << Form("pzf_%d",i) << ",";
    }
    
    csv_file
    << "neu" << ","
    << "Pv_x" << ","
    << "Pv_y" << ","
    << "Pv_z" << ","
    << "Pv_theta" << ","
    << "truth_Ev" << ","
    ;
    
    // my variables for CC1p0pi
    csv_file
    << "IsCC_Np_200MeVc"    << ","
    << "IsCC_1p_200MeVc"    << ","
    << "IsCC_1p_200MeVc_0pi"<< ","
    << "truth_Pmu"        << "," << "truth_Pmu_theta"   << "," << "truth_Pmu_cos_theta"   << ","
    << "truth_Pmu_phi"   << ","
    << "truth_Pmu_x"      << "," << "truth_Pmu_y"       << "," << "truth_Pmu_z"   << ","
    << "truth_Pp"         << "," << "truth_Pp_theta"    << "," << "truth_Pp_cos_theta"    << ","
    << "truth_Pp_phi"    << ","
    << "truth_Pp_x"       << "," << "truth_Pp_y"        << "," << "truth_Pp_z"    << ","
    << "Pmiss"      << ","
    << "Pmiss_x"    << "," << "Pmiss_y"     << "," << "Pmiss_z" << ","
    << "alpha_p"    << ","
    << "alpha_q"    << ","
    << "alpha_miss"
    << "Pp_before_FSI"         << "," << "Pp_before_FSI_theta"    << ","
    << "Pp_before_FSI_x"       << "," << "Pp_before_FSI_y"        << "," << "Pp_before_FSI_z" << ",";
    
    // weight to the event based on MicroBooNE acceptance
    csv_file
    << "uBacc_truth_Q2"     << ","
    << "uBacc_truth_muon"   << ","  << "uBacc_truth_proton"     << ","   << "uBacc_truth_muon_proton" << ","
    << "uBacc_reco_Q2"      << ","
    << "uBacc_reco_muon"    << ","  << "uBacc_reco_proton"      << ","   << "uBacc_muon_reco_proton" << ",";
    
    
    
    // mimic the detector volume
    // and "ruin" the kinematical reconstruction of the muon momentum and Q2
    csv_file
    << "v_x"                << "," << "v_y"                 << "," << "v_z"                     << ","
    << "muonTrajectory_endx"<< "," << "muonTrajectory_endy" << "," << "muonTrajectory_endz"     << ","
    << "muonTrack_endx"     << "," << "muonTrack_endy"      << "," << "muonTrack_endz"          << ","
    << "reco_Emu"           << "," << "reco_Pmu"            << ","
    << "reco_Pmu_x"         << "," << "reco_Pmu_y"          << "," << "reco_Pmu_z"              << ","
    << "reco_Pmu_theta"     << "," << "reco_Pmu_phi"        << ","
    << "reco_Ep"            << "," << "reco_Pp"             << ","
    << "reco_Pp_x"          << "," << "reco_Pp_y"           << "," << "reco_Pp_z"               << ","
    << "reco_Pp_theta"      << "," << "reco_Pp_phi"         << ","
    << "reco_q"             << "," << "reco_omega"          << "," << "reco_Q2"                 << ","
    << "reco_l_mu"          << "," << "reco_l_p"            << ","
    << "reco_delta_phi"     << "," << "reco_Pt"             << ","
    << "gen_Q2_gen_rec"     << "," << "rec_Q2_gen_rec"      << ",";


    
    // start/end points, for FV cuts
    csv_file
    << "startx_muCandidate" << ","  << "starty_muCandidate" << ","  << "startz_muCandidate" << ","
    << "startx_pCandidate"  << ","  << "starty_pCandidate"  << ","  << "startz_pCandidate"  << ","
    << "endx_muCandidate"   << ","  << "endy_muCandidate"   << ","  << "endz_muCandidate"   << ","
    << "endx_pCandidate"    << ","  << "endy_pCandidate"    << ","  << "endz_pCandidate"    << ",";

    // kinematical cuts
    csv_file
    << "delta_phi"          << ","  << "Pt"                 << ","  << "theta_12";
    
    csv_file << endl;
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::StreamToCSV (){
    
    counter++;
    
    csv_file
    << cc  << "," << qel <<  "," << res <<  "," << dis <<  "," << coh <<  ","
    << mec << ","
    << Q2  << "," << W <<  "," << x <<  ","
    << pxn << "," << pyn <<  "," << pzn <<  "," // Initial state hit nucleon
    << pxv << "," << pyv <<  "," << pzv <<  ","
    << pxl << "," << pyl <<  "," << pzl << ","
    << nf  << ",";
    
    for (int i=0; i<NMAX ; i++ ){
        csv_file
        << pdgf[i] <<  ","
        << pxf[i] <<  "," << pyf[i] <<  "," << pzf[i] << ",";
    }
    
    csv_file
    << neu      << ","
    << nu.Px()  << ","
    << nu.Py()  << ","
    << nu.Pz()  << ","
    << nu.Theta()  << ","
    << nu.E()   << ",";
    
    // my variables for CC1p0pi
    csv_file
    << IsCC_Np_200MeVc  << ","
    << IsCC_1p_200MeVc  << ","
    << IsCC_1p_200MeVc_0pi  << ","
    << muon.P()     << "," << muon.Theta()      << "," << TMath::Cos(muon.Theta())  << ","
    << muon.Phi()   << ","
    << muon.Px()    << "," << muon.Py()         << "," << muon.Pz()     << ","
    << proton.P()   << "," << proton.Theta()    << "," << TMath::Cos(proton.Theta())  << ","
    << proton.Phi()   << ","
    << proton.Px()  << "," << proton.Py()       << "," << proton.Pz()   << ","
    << Pmiss.P()    << ","
    << Pmiss.Px()   << "," << Pmiss.Py()        << "," << Pmiss.Pz()    << ","
    << LightConeAlpha( proton ) << ","
    << LightConeAlpha( q )      << ","
    << LightConeAlpha( proton ) - LightConeAlpha( q )   
    << proton_before_FSI.P()   << "," << proton_before_FSI.Theta()    << ","
    << proton_before_FSI.Px()  << "," << proton_before_FSI.Py()       << "," << proton_before_FSI.Pz() << ",";
    
    // weight to the event based on MicroBooNE acceptance
    csv_file
    << uBacc_truth_Q2           << ","
    << uBacc_truth_muon         << "," << uBacc_truth_proton            << "," << uBacc_truth_muon * uBacc_truth_proton << ","
    << uBacc_reco_Q2            << ","
    << uBacc_reco_muon          << "," << uBacc_reco_proton            << "," << uBacc_reco_muon * uBacc_reco_proton << ",";
    
    
    
    // mimic the detector volume
    // and "ruin" the kinematical reconstruction of the muon momentum and Q2
    csv_file
    << vertex_position.X()      << "," << vertex_position.Y()       << "," << vertex_position.Z()       << ","
    << muonTrajectory_end.X()   << "," << muonTrajectory_end.Y()    << "," << muonTrajectory_end.Z()    << ","
    << muonTrack_end.X()        << "," << muonTrack_end.Y()         << "," << muonTrack_end.Z()         << ","
    << reco_Pmu.E()             << "," << reco_Pmu.P()              << ","
    << reco_Pmu.Px()            << "," << reco_Pmu.Py()             << "," << reco_Pmu.Pz()             << ","
    << reco_Pmu.Theta()         << "," << reco_Pmu.Phi()            << ","
    << reco_Pp.E()              << "," << reco_Pp.P()               << ","
    << reco_Pp.Px()             << "," << reco_Pp.Py()              << "," << reco_Pp.Pz()             << ","
    << reco_Pp.Theta()          << "," << reco_Pp.Phi()             << ","
    << reco_q.P()               << "," << reco_q.E()                << "," << reco_Q2                   << ","
    << reco_l_mu                << "," << reco_l_p                  << ","
    << reco_delta_phi           << "," << reco_Pt                   << ","
    << gen_Q2_gen_rec           << "," << rec_Q2_gen_rec            << ",";

    // start/end points, for FV cuts
    csv_file
    << vertex_position.x()  << ","  << vertex_position.y()  << ","  << vertex_position.z()  << ","
    << vertex_position.x()  << ","  << vertex_position.y()  << ","  << vertex_position.z()  << ","
    << muonTrack_end.x()    << ","  << muonTrack_end.y()    << ","  << muonTrack_end.z()    << ","
    << protonTrack_end.x()  << ","  << protonTrack_end.y()  << ","  << protonTrack_end.z()  << ",";
    
    // kinematical cuts
    csv_file
    << r2d*fabs(proton.Phi() - muon.Phi())                  << ","  << (proton.Vect() + muon.Vect()).Pt()  << "," <<     r2d*(proton.Vect().Angle( muon.Vect() ));

    
    csv_file << endl;
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::SetTopology (){
    nu.SetXYZM( pxv, pyv, pzv, 0 );
    if (debug>5) {
        SHOW4(cc,qel,nf,pdgf[0]);
        for (int i=0; i<nf; i++) {
            SHOW(pdgf[i]);
        }
    }
    
    
    
    int Nmu=0, Npi=0, Nel=0, Ngamma=0, Npi_0=0;
    int Np_200MeVc=0, Nn_200MeVc=0, Npi_70MeVc=0;
    
    muon.SetXYZM( pxl, pyl, pzl, 0.1056583745 );

    for( int i_f=0; i_f < nf ; i_f++ ){
        Debug(5,"pdgf[i_f]:%",pdgf[i_f]);
        switch (pdgf[i_f]) {
            case 2212:
                pmomentum.SetXYZM(pxf[i_f], pyf[i_f], pzf[i_f], 0.9383);
                protons.push_back( pmomentum );
                if (pmomentum.P() >= 0.2){
                    Np_200MeVc++;
                }
                if (protons.back().P() > proton.P()) {
                    proton = protons.back();
                }
                break;
                
            case 2112:
                pmomentum.SetXYZM(pxf[i_f], pyf[i_f], pzf[i_f], 0.9395);
                neutrons.push_back( pmomentum );
                if (pmomentum.P() >= 0.2){
                    Nn_200MeVc++;
                }
                break;
                
            case 211:
                Npi++;
                if (pmomentum.P() >= 0.07){
                    Npi_70MeVc++;
                }
                break;
            case -211:
                Npi++;
                if (pmomentum.P() >= 0.07){
                    Npi_70MeVc++;
                }
                break;
            case 111:
                Npi++;
                Npi_0++;
                break;

            case 22:
                Ngamma++;
                break;
                
            case 11:
            case -11:
                Nel++;
                break;
                
            default:
                break;
        }
    }
        
    
    
    // CCQE with FSI (the A-1 system can not exit the nucleus)
    Debug(2,"cc:% , qel:%, Nmu:%, Np_200MeVc:%, Npi:%, Nn_200MeVc:%, Nel:%, Ngamma:%"
          ,cc,qel,Nmu,Np_200MeVc,Npi,Nn_200MeVc,Nel,Ngamma);
    // IsCC_Np_200MeVc
    if ( cc==true
        && Np_200MeVc>=1 ){
        IsCC_Np_200MeVc = true;
        
        // IsCC_1p_200MeVc
        // an interaction with at least 1 muon and exactly 1 proton > 200 MeV/c and no charged pions with momentum greater than 70 MeV/c, and no neutral pions in the final state
        if ( Np_200MeVc==1 && Npi_70MeVc==0 && Npi_0==0){
            IsCC_1p_200MeVc = true;
            
            // IsCC_1p_200MeVc_0pi
            // an interaction with at least 1 muon and 1 proton > 200 MeV/c and no pions, and no photons or electrons outside the nucleus
            // which means that the final state includes
            // 1 protons with momentum >= 200 MeV/c
            // any number of neutrons
            // 0 pions
            // 0 electrons
            // 0 photons
            if ( Npi==0 && Nel==0 && Ngamma==0 ){
                IsCC_1p_200MeVc_0pi = true;
                q = nu - muon;
                Q2 = -q.Mag2();
                Pmiss = proton - q;
            }
        }
    }

    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int GenieFile::FindWhichBin( double x, std::vector<double> bins ){
    for (size_t i=0; i<bins.size()-1; i++) {
        if ( bins.at(i) < x && x < bins.at(i+1) ){
            Debug(5,"GenieFile::FindWhichBin( %, std::vector<double> ): bin %",x,i);
            return i;
        }
    }
    return -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::SetMicroBooNEWeights (){
    int fDebug=3;
    Debug(fDebug,"GenieFile::SetMicroBooNEWeights");
    Debug(fDebug,"p(mu): %, 𝜃(mu): %π,p(p): %, 𝜃(p): %π",muon.P(),muon.Theta()/3.1415,proton.P(),proton.Theta()/3.1415);
    
    // assign a weight to the event based on MicroBooNE acceptance
    // as a function of the muon momentum and scattering angle
    //
    // To account for the uncertainties in the acceptance,
    // weight generated from a Gaussian distribution fuction
    // with the mean being the acceptance from the overlay
    // and the sigma being the uncertainty in this acceptance
    //
    // multiple the weights by 1.0e3 is to avoid from very small weights...
    
    // MicroBooNE weights as a function of the truth kinematics
    uBacc_truth_Q2 = uBacc_truth_muon = uBacc_truth_proton = 0;
    
    // (1) find the xbin and ybin of the event
    int Pmu_bin = FindWhichBin( muon.P(), Pmu_xbins );
    int Pmu_theta_bin = FindWhichBin( muon.Theta(), Pmu_theta_ybins );
    Debug(fDebug,"Pmu_bin: %, Pmu_theta_bin: %",Pmu_bin,Pmu_theta_bin);
    // (2) apply the weight
    if ( Pmu_bin<0 || Pmu_theta_bin<0 ) { uBacc_truth_muon = 0; } else {
        double mean  = Pmu_theta_acceptance.at(Pmu_theta_bin).at(Pmu_bin);
        Debug(fDebug+1,"MicroBooNE Weight[Pmu-theta bin %][Pmu bin %]: %",Pmu_theta_bin,Pmu_bin,mean);
        double sigma = Pmu_theta_acc_err.at(Pmu_theta_bin).at(Pmu_bin);
        uBacc_truth_muon = rand.Gaus( mean , sigma ) * 1.e3;
    }

    
    
    int Pp_bin = FindWhichBin( proton.P(), Pp_xbins );
    int Pp_theta_bin = FindWhichBin( proton.Theta(), Pp_theta_ybins );
    Debug(fDebug,"Pp_bin: %, Pp_theta_bin: %",Pp_bin,Pp_theta_bin);
    if ( Pp_bin<0 || Pp_theta_bin<0 ) { uBacc_truth_proton = 0; } else {
        double mean  = Pp_theta_acceptance.at(Pp_theta_bin).at(Pp_bin);
        Debug(fDebug+1,"MicroBooNE Weight[Pp-theta bin %][Pp bin %]: %",Pp_theta_bin,Pp_bin,mean);
        double sigma = Pp_theta_acc_err.at(Pp_theta_bin).at(Pp_bin);
        uBacc_truth_proton = rand.Gaus( mean , sigma ) * 1.e3;
    }
    
    
    
    int Q2_bin = FindWhichBin( Q2 , Q2_bins );
    if ( Q2_bin<0 ) { uBacc_truth_Q2 = 0;    } else {
        double mean  = Q2_acceptance.at(Q2_bin);
        Debug(fDebug+1,"MicroBooNE Weight[Q2 bin %]: %",Q2_bin,mean);
        double sigma = Q2_acc_err.at(Q2_bin);
        uBacc_truth_Q2 = rand.Gaus( mean , sigma ) * 1.e3;
    }
    
    Debug(fDebug,"muon weight: %, proton weight: %",uBacc_truth_muon,uBacc_truth_proton);
    Debug(fDebug,"Q2 weight: %",uBacc_truth_Q2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double GenieFile::LinePlaneIntersectionDistance(TVector3 l,TVector3 l0,TVector3 p0,TVector3 n,double epsilon){
    // epsilon is a cutoff for the case of trajectory parallel to either one of the detector sides
    // [https://en.wikipedia.org/wiki/Line–plane_intersection]
    // where l is a vector in the direction of the line, l0 is a point on the line,
    // n is the normal to the plane and p0 is a point on the plan
    Debug(4, "GenieFile::LinePlaneIntersectionDistance(), n: (%,%,%),  l.Dot(n): %" ,n.x(),n.y(),n.z(),l.Dot(n));
    if (fabs(l.Dot(n)) < epsilon){
        return ( ((p0 - l0).Dot(n))/epsilon );
    }
    return ( ((p0 - l0).Dot(n))/l.Dot(n) );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::ProjectMuonTrajectory(){
    int fDebug=4;
    Debug(fDebug,"GenieFile::ProjectMuonTrajectory()");
    
    // the muon direction is given by its momentum...
    muonTrajectory_dir = muon.Vect().Unit();
    double Pmu_MeVc = 1000*muon.P();
    dmu_cm = Get_muonRangeFromMomentum( Pmu_MeVc );
    muonTrajectory_end = vertex_position + muonTrajectory_dir * dmu_cm;
    
    
    
    Debug(fDebug,"Pmu: % MeV/c, lmu: % cm",Pmu_MeVc,dmu_cm);
    Debug(fDebug,"muon trajectory: (%,%,%) -> (%,%,%)"
          ,vertex_position.X(),vertex_position.Y(),vertex_position.Z()
          ,muonTrajectory_end.X(),muonTrajectory_end.Y(),muonTrajectory_end.Z());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::CutMuonTrajectory(){
    int fDebug=3;
    
    muonTrack_end = muonTrajectory_end; // start with 'generated' end point
    
    double dx=1000,dy=1000,dz=1000,d_min=1000; // muon end = start + direction * d;
    while (
           (side_x_left-1 > muonTrack_end.X())      || (muonTrack_end.X() > side_x_right+1)
           ||
           (side_y_dw-1 > muonTrack_end.Y())        || (muonTrack_end.Y() > side_y_up+1)
           ||
           (side_z_downstream-1 > muonTrack_end.Z())|| (muonTrack_end.Z() > side_z_upstream+1)
           ) {
        
        // check from which side does the track exit, add 1 cm for ~ "thickness"
        if ( (muonTrajectory_dir.X()>0) && (side_x_right+1 < muonTrack_end.X()) ) {
            dx = LinePlaneIntersectionDistance( muonTrajectory_dir , vertex_position
                                        , TVector3( side_x_right , 0 , 0 ) , TVector3( 1 , 0 , 0 ));
            Debug(fDebug,"trajectory exited the detector from the right side, dx=%",dx);
        }
        else if ( (muonTrajectory_dir.X()<0) && (side_x_left-1 > muonTrack_end.X()) ){
            dx = LinePlaneIntersectionDistance( muonTrajectory_dir , vertex_position
                                        , TVector3( side_x_left , 0 , 0 ) , TVector3( 1 , 0 , 0 ));
            Debug(fDebug,"trajectory exited the detector from the left side, dx=%",dx);
        }
        
        if ( (muonTrajectory_dir.Y()>0) && (side_y_up+1 < muonTrack_end.Y()) ) {
            dy = LinePlaneIntersectionDistance( muonTrajectory_dir , vertex_position
                                        , TVector3( 0 , side_y_up , 0 ) , TVector3( 0 , 1 , 0 ));
            Debug(fDebug,"trajectory exited the detector from the top side, dy=%",dy);
        }
        else if ( (muonTrajectory_dir.Y()<0) && (side_y_dw-1 > muonTrack_end.Y()) ){
            dy = LinePlaneIntersectionDistance( muonTrajectory_dir , vertex_position
                                        , TVector3( 0 , side_y_dw , 0 ) , TVector3( 0 , 1 , 0 ));
            Debug(fDebug,"trajectory exited the detector from the bottom side, dy=%",dy);
        }
        
        if ( (muonTrajectory_dir.Z()>0) && (side_z_upstream+1 < muonTrack_end.Z()) ) {
            dz = LinePlaneIntersectionDistance( muonTrajectory_dir , vertex_position
                                        , TVector3( 0 , 0 , side_z_upstream ) , TVector3( 0  , 0 , 1 ));
            Debug(fDebug,"trajectory exited the detector from the top side, dz=%",dz);
        }
        else if ( (muonTrajectory_dir.Z()<0) && (side_z_downstream-1 > muonTrack_end.Z()) ){
            dz = LinePlaneIntersectionDistance( muonTrajectory_dir , vertex_position
                                        , TVector3( 0 , 0 , side_z_downstream ) , TVector3( 0 , 0 , 1 ));
            Debug(fDebug,"trajectory exited the detector from the bottom side, dz=%",dz);
        }
        
        if (fabs(dx)<fabs(dy) && fabs(dx)<fabs(dz)) {
            d_min = dx;
        }
        else if (fabs(dy)<fabs(dx) && fabs(dy)<fabs(dz)) {
            d_min = dy;
        }
        else{
            d_min = dz;
        }
        Debug(fDebug,"d_min: %",d_min);

        muonTrack_end = vertex_position + muonTrajectory_dir * (0.9999 * d_min);
        Debug(fDebug,"muon track direction: (%,%,%)",muonTrajectory_dir.x(),muonTrajectory_dir.y(),muonTrajectory_dir.z());
        
        Debug(fDebug,"changed muon track to: (%,%,%) -> (%,%,%)"
              ,vertex_position.X(),vertex_position.Y(),vertex_position.Z()
              ,muonTrack_end.X(),muonTrack_end.Y(),muonTrack_end.Z());
    }
    Debug(fDebug,"muon track: (%,%,%) -> (%,%,%)"
          ,vertex_position.X(),vertex_position.Y(),vertex_position.Z()
          ,muonTrack_end.X(),muonTrack_end.Y(),muonTrack_end.Z());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::ProjectProtonTrajectory(){
    int fDebug=4;
    Debug(fDebug,"GenieFile::ProjectProtonTrajectory()");
    
    // the muon direction is given by its momentum...
    protonTrajectory_dir = proton.Vect().Unit();
    double Pp_MeVc = 1000*proton.P();
    dp_cm = Get_protonRangeFromMomentum( Pp_MeVc );
    protonTrajectory_end = vertex_position + protonTrajectory_dir * dp_cm;
    
    Debug(fDebug,"Pp: % MeV/c, lp: % cm",Pp_MeVc,dp_cm);
    Debug(fDebug,"proton trajectory: (%,%,%) -> (%,%,%)"
          ,vertex_position.X(),vertex_position.Y(),vertex_position.Z()
          ,protonTrajectory_end.X(),protonTrajectory_end.Y(),protonTrajectory_end.Z());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::CutProtonTrajectory(){
    int fDebug=3;
    
    protonTrack_end = protonTrajectory_end; // start with 'generated' end point
    
    double dx=1000,dy=1000,dz=1000,d_min=1000; // proton end = start + direction * d;
    while (
           (side_x_left-1 > protonTrack_end.X())      || (protonTrack_end.X() > side_x_right+1)
           ||
           (side_y_dw-1 > protonTrack_end.Y())        || (protonTrack_end.Y() > side_y_up+1)
           ||
           (side_z_downstream-1 > protonTrack_end.Z())|| (protonTrack_end.Z() > side_z_upstream+1)
           ) {
        
        // check from which side does the track exit, add 1 cm for ~ "thickness"
        if ( (protonTrajectory_dir.X()>0) && (side_x_right+1 < protonTrack_end.X()) ) {
            dx = LinePlaneIntersectionDistance( protonTrajectory_dir , vertex_position
                                               , TVector3( side_x_right , 0 , 0 ) , TVector3( 1 , 0 , 0 ));
            Debug(fDebug,"trajectory exited the detector from the right side, dx=%",dx);
        }
        else if ( (protonTrajectory_dir.X()<0) && (side_x_left-1 > protonTrack_end.X()) ){
            dx = LinePlaneIntersectionDistance( protonTrajectory_dir , vertex_position
                                               , TVector3( side_x_left , 0 , 0 ) , TVector3( 1 , 0 , 0 ));
            Debug(fDebug,"trajectory exited the detector from the left side, dx=%",dx);
        }
        
        if ( (protonTrajectory_dir.Y()>0) && (side_y_up+1 < protonTrack_end.Y()) ) {
            dy = LinePlaneIntersectionDistance( protonTrajectory_dir , vertex_position
                                               , TVector3( 0 , side_y_up , 0 ) , TVector3( 0 , 1 , 0 ));
            Debug(fDebug,"trajectory exited the detector from the top side, dy=%",dy);
        }
        else if ( (protonTrajectory_dir.Y()<0) && (side_y_dw-1 > protonTrack_end.Y()) ){
            dy = LinePlaneIntersectionDistance( protonTrajectory_dir , vertex_position
                                               , TVector3( 0 , side_y_dw , 0 ) , TVector3( 0 , 1 , 0 ));
            Debug(fDebug,"trajectory exited the detector from the bottom side, dy=%",dy);
        }
        
        if ( (protonTrajectory_dir.Z()>0) && (side_z_upstream+1 < protonTrack_end.Z()) ) {
            dz = LinePlaneIntersectionDistance( protonTrajectory_dir , vertex_position
                                               , TVector3( 0 , 0 , side_z_upstream ) , TVector3( 0  , 0 , 1 ));
            Debug(fDebug,"trajectory exited the detector from the top side, dz=%",dz);
        }
        else if ( (protonTrajectory_dir.Z()<0) && (side_z_downstream-1 > protonTrack_end.Z()) ){
            dz = LinePlaneIntersectionDistance( protonTrajectory_dir , vertex_position
                                               , TVector3( 0 , 0 , side_z_downstream ) , TVector3( 0 , 0 , 1 ));
            Debug(fDebug,"trajectory exited the detector from the bottom side, dz=%",dz);
        }
        
        if (fabs(dx)<fabs(dy) && fabs(dx)<fabs(dz)) {
            d_min = dx;
        }
        else if (fabs(dy)<fabs(dx) && fabs(dy)<fabs(dz)) {
            d_min = dy;
        }
        else{
            d_min = dz;
        }
        Debug(fDebug,"d_min: %",d_min);
        protonTrack_end = vertex_position + protonTrajectory_dir * (0.9999 * d_min);
        Debug(fDebug,"proton track direction: (%,%,%)",protonTrajectory_dir.x(),protonTrajectory_dir.y(),protonTrajectory_dir.z());
        
        Debug(fDebug,"changed proton track to: (%,%,%) -> (%,%,%)"
              ,vertex_position.X(),vertex_position.Y(),vertex_position.Z()
              ,protonTrack_end.X(),protonTrack_end.Y(),protonTrack_end.Z());
    }
    Debug(fDebug,"proton track: (%,%,%) -> (%,%,%)"
          ,vertex_position.X(),vertex_position.Y(),vertex_position.Z()
          ,protonTrack_end.X(),protonTrack_end.Y(),protonTrack_end.Z());
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int GenieFile::SampleFromDistribution( std::vector<double> probabilities ){
    // sample from an array of probabilities
    while (1) { // eternal loop until we find a good index and return it....
        // (A) randomly sample an index from a uniform distribution
        int i = rand.Uniform( 0 , probabilities.size() );
        // (B) sample height from a uniform distribution funtion
        // and if its lower or equal the probability keep the index
        double height = rand.Uniform(0 , 1);
        if ( height <= probabilities.at(i) ) {
            return i;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::SetRecoKinematics(){
    reco_l_mu   = (muonTrack_end-vertex_position).Mag();
    reco_l_p    = (protonTrack_end-vertex_position).Mag();
    
    reco_Pmu = muon * (reco_l_mu/dmu_cm);
    reco_Pp  = proton * (reco_l_p/dp_cm);

    reco_Ev  = reco_Pmu.E() + (reco_Pp.E() - reco_Pp.Mag()) + 0.040 ; // Eµ + Tp + Sn + T(A-1)
    reco_Pnu = TLorentzVector( 0 , 0 , reco_Ev , reco_Ev );

    reco_q   = reco_Pnu - reco_Pmu;
    reco_Q2  = - reco_q.Mag2();
    reco_delta_phi = r2d*fabs(reco_Pmu.Phi() - reco_Pp.Phi());
    reco_Pt  = (reco_Pmu + reco_Pp).Pt();
    
    
    // reconstructed Q2 from the overlay map
    i_gen = FindWhichBin( Q2 , Q2_gen_rec_bins );
    Debug(4,"Q2: %, i_gen: %",Q2,i_gen);
    if ( Q2 < Q2_gen_rec_bins.at(0) ) { // under bin
        Debug(4,"Q2 < Q2_gen_rec_bins.at(0)!!");
        gen_Q2_gen_rec = -1;
        rec_Q2_gen_rec = -1;
    }
    else if ( Q2_gen_rec_bins.at( Q2_gen_rec_bins.size() - 1 ) < Q2 ) { // over bin
        Debug(4,"Q2_gen_rec_bins.at( Q2_gen_rec_bins.size() - 1 ) < Q2!!");
        gen_Q2_gen_rec = 9999;
        rec_Q2_gen_rec = 9999;
    }
    else if (i_gen == -1){
        Debug(4,"i_gen == -1!!");
        gen_Q2_gen_rec = -9999;
        rec_Q2_gen_rec = -9999;
    }
    else{
        gen_Q2_gen_rec = 0.5*(Q2_gen_rec_bins.at(i_gen) + Q2_gen_rec_bins.at(i_gen+1));
        i_rec = Q2_rec_distributions.at(i_gen)(generator);
        rec_Q2_gen_rec = 0.5*(Q2_gen_rec_bins.at(i_rec) + Q2_gen_rec_bins.at(i_rec+1));
    }
    Debug(4,"gen_Q2_gen_rec: %, i_rec=%, rec_Q2_gen_rec: %",gen_Q2_gen_rec,i_rec,rec_Q2_gen_rec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::MimicDetectorVolume(){
    
    // mimic the effect of the detector volume on the "reconstructed" kinematics
    // by projecting the muon track length from range-based momentum in LAr
    // and cutting it short if it touches the detector edges
    
    // vtxx is given between -1.282 and 1.282 m
    // vtxy is given between -1.165 and 1.165 m
    // vtxz is given between -5.184 and 5.184 m
    // we want to vertex position in cm
    SetVertexPosition( 100*(vtxx + 1.282) , 100*vtxy , 100*(vtxz + 5.184));
    
    ProjectMuonTrajectory();
    //    CutMuonTrajectory();
    ProjectProtonTrajectory();
    //    CutProtonTrajectory();    
    // SetRecoKinematics();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::EndJob(){
    csv_file.close();
    cout << "wrote " << counter << " events to csv file:" << endl << Path+OutputCSVname << endl;
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::Print(){
    SHOW(Q2);
    EndEventBlock();
    return true;
}

#endif

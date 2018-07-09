#ifndef GENIEFILE_CXX
#define GENIEFILE_CXX

#include "GenieFile.h"

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
    Debug(0,"GenieFile::Read1dArrayFromFile(%)",filename);
    ifstream fin;
    fin.open(filename);
    std::vector<double> numbers;
    double number;
    while(fin >> number) //  >> delim
        numbers.push_back(number);
    fin.close();
    if(debug>0){
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
    Debug(0,"GenieFile::Read2dArrayFromFile(%)",filename);
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
    if(debug>0){
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
//    dx = dy = dz = 0;
    reco_Q2 = reco_Ev = 0;
    reco_Pnu = reco_Pp = reco_Pmu = reco_q = TLorentzVector();
    
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
    << "Q2"  << "," << "W" <<  "," << "x" <<  ","
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
    << "Ev" << ","
    ;
    
    // my variables for CC1p0pi
    csv_file
    << "IsCC_Np_200MeVc"    << ","
    << "IsCC_1p_200MeVc"    << ","
    << "IsCC_1p_200MeVc_0pi"<< ","
    << "Pmu"        << "," << "Pmu_theta"   << ","
    << "Pmu_x"      << "," << "Pmu_y"       << "," << "Pmu_z"   << ","
    << "Pp"         << "," << "Pp_theta"    << ","
    << "Pp_x"       << "," << "Pp_y"        << "," << "Pp_z"    << ","
    << "Pmiss"      << ","
    << "Pmiss_x"    << "," << "Pmiss_y"     << "," << "Pmiss_z" << ","
    << "alpha_p"    << ","
    << "alpha_q"    << ","
    << "alpha_miss"
    << "Pp_before_FSI"         << "," << "Pp_before_FSI_theta"    << ","
    << "Pp_before_FSI_x"       << "," << "Pp_before_FSI_y"        << "," << "Pp_before_FSI_z" << ",";
    
    // weight to the event based on MicroBooNE acceptance
    csv_file
    << "MicroBooNEWeight_Q2" << ","
    << "MicroBooNEWeight_Pmu_theta" << ","
    << "MicroBooNEWeight_Pp_theta"  << ","
    << "MicroBooNEWeight_Pmu_theta_Pp_theta" << ",";
    
    
    
    // mimic the detector volume
    // and "ruin" the kinematical reconstruction of the muon momentum and Q2
    csv_file
    << "v_x"                << "," << "v_y"                 << "," << "v_z"                     << ","
    << "muonTrajectory_endx"<< "," << "muonTrajectory_endy" << "," << "muonTrajectory_endz"     << ","
    << "muonTrack_endx"     << "," << "muonTrack_endy"      << "," << "muonTrack_endz"          << ","
    << "reco_Emu"           << "," << "reco_Pmu"            << ","
    << "reco_Pmu_x"         << "," << "reco_Pmu_y"          << "," << "reco_Pmu_z"              << ","
    << "reco_q"             << "," << "reco_omega"          << "," << "reco_Q2";

    
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
    << muon.P()     << "," << muon.Theta()      << ","
    << muon.Px()    << "," << muon.Py()         << "," << muon.Pz()     << ","
    << proton.P()   << "," << proton.Theta()    << ","
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
    << MicroBooNEWeight_Q2                  << ","
    << MicroBooNEWeight_Pmu_theta           << ","
    << MicroBooNEWeight_Pp_theta            << ","
    << (MicroBooNEWeight_Pmu_theta * MicroBooNEWeight_Pp_theta ) << ",";
    
    
    
    // mimic the detector volume
    // and "ruin" the kinematical reconstruction of the muon momentum and Q2
    csv_file
    << vertex_position.X()      << "," << vertex_position.Y()       << "," << vertex_position.Z()       << ","
    << muonTrajectory_end.X()   << "," << muonTrajectory_end.Y()    << "," << muonTrajectory_end.Z()    << ","
    << muonTrack_end.X()        << "," << muonTrack_end.Y()         << "," << muonTrack_end.Z()         << ","
    << reco_Pmu.E()             << "," << reco_Pmu.P()              << ","
    << reco_Pmu.Px()            << "," << reco_Pmu.Py()             << "," << reco_Pmu.Pz()             << ","
    << reco_q.P()               << "," << reco_q.E()                << "," << reco_Q2;

    
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
    
    
    
    int Nmu=0, Npi=0, Nel=0, Ngamma=0;
    int Np_200MeVc=0, Nn_200MeVc=0;
    
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
            case -211:
            case 111:
                Npi++;
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
        //        && ( Nmu>=1 ) // if cc==0 then we have a muon!
        && Np_200MeVc>=1 ){
        IsCC_Np_200MeVc = true;
        
        // IsCC_1p_200MeVc
        if ( Np_200MeVc==1 ){
            IsCC_1p_200MeVc = true;
            
            // IsCC_1p_200MeVc_0pi
            // which means that the final state includes
            // 1 protons with momentum >= 200 MeV/c
            // any number of neutrons
            // 0 pions
            // 0 electrons
            // 0 photons
            if ( Npi==0 && Nel==0 && Ngamma==0 ){
                IsCC_1p_200MeVc_0pi = true;
                q = nu - muon;
                Pmiss = proton - q;
//                // Delete this!
//                CC1p0pi = true;
            }
        }
    }

    
//    
//    if ((cc==true)
//        && ( Np_200MeVc==1 )
//        //        && ( Nmu>=1 ) // if cc==0 then we have a muon!
//        && ( Npi==0 && Nel==0 && Ngamma==0 )
//        && Nn_200MeVc==0
//        //        && ((nf == 1) && (pdgf[0]==2212) )
//        )
//    {
//        CC1p0pi = true;
//        q = nu - muon;
//        Pmiss = proton - q;
//        // proton_before_FSI.SetXYZM( pxi[0],pyi[0],pzi[0], 0.9382720813 ); // correct this!
//    }
    
//    // CCQE with no FSI
//    else if ((cc==true)
//             && (nf == 2) && (pdgf[0]==1000180390) && (pdgf[1]==2212) ){
//        CC1p0pi = true;
//        muon.SetXYZM( pxl, pyl, pzl, 0.1056583745 );
//        q = nu - muon;
//        proton.SetXYZM( pxf[1], pyf[1], pzf[1], 0.9382720813 );
//        Pmiss = proton - q;
//        proton_before_FSI.SetXYZM( pxi[1],pyi[1],pzi[1], 0.9382720813 );
//    }
    
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
void GenieFile::SetMicroBooNEWeight (){
    int fDebug=3;
    Debug(fDebug,"GenieFile::SetMicroBooNEWeight");
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
    
    MicroBooNEWeight_Q2 = MicroBooNEWeight_Pmu_theta = MicroBooNEWeight_Pp_theta = 0;
    
    // (1) find the xbin and ybin of the event
    int Pmu_bin = FindWhichBin( muon.P(), Pmu_xbins );
    int Pmu_theta_bin = FindWhichBin( muon.Theta(), Pmu_theta_ybins );
    Debug(fDebug,"Pmu_bin: %, Pmu_theta_bin: %",Pmu_bin,Pmu_theta_bin);
    // (2) apply the weight
    if ( Pmu_bin<0 || Pmu_theta_bin<0 ) {
        MicroBooNEWeight_Pmu_theta = 0;
    } else {
        double mean  = Pmu_theta_acceptance.at(Pmu_theta_bin).at(Pmu_bin);
        Debug(fDebug+1,"MicroBooNE Weight[Pmu-theta bin %][Pmu bin %]: %",Pmu_theta_bin,Pmu_bin,mean);
        double sigma = Pmu_theta_acc_err.at(Pmu_theta_bin).at(Pmu_bin);
        MicroBooNEWeight_Pmu_theta = rand.Gaus( mean , sigma ) * 1.e3;
    }

    
    
    // (1) find the xbin and ybin of the event
    int Pp_bin = FindWhichBin( proton.P(), Pp_xbins );
    int Pp_theta_bin = FindWhichBin( proton.Theta(), Pp_theta_ybins );
    Debug(fDebug,"Pp_bin: %, Pp_theta_bin: %",Pp_bin,Pp_theta_bin);
    
    // (2) apply the weight
    if ( Pp_bin<0 || Pp_theta_bin<0 ) {
        MicroBooNEWeight_Pp_theta = 0;
    } else {
        double mean  = Pp_theta_acceptance.at(Pp_theta_bin).at(Pp_bin);
        Debug(fDebug+1,"MicroBooNE Weight[Pp-theta bin %][Pp bin %]: %",Pp_theta_bin,Pp_bin,mean);
        double sigma = Pp_theta_acc_err.at(Pp_theta_bin).at(Pp_bin);
        MicroBooNEWeight_Pp_theta = rand.Gaus( mean , sigma ) * 1.e3;
    }
    
    
    
    // (1) find the Q2-bin
    int Q2_bin = FindWhichBin( Q2 , Q2_bins );
    Debug(fDebug,"Q2_bin: %",Q2_bin);
    
    // (2) apply the weight
    if ( Q2_bin<0 ) {
        MicroBooNEWeight_Q2 = 0;
    } else {
        double mean  = Q2_acceptance.at(Q2_bin);
        Debug(fDebug+1,"MicroBooNE Weight[Q2 bin %]: %",Q2_bin,mean);
        double sigma = Q2_acc_err.at(Q2_bin);
        MicroBooNEWeight_Q2 = rand.Gaus( mean , sigma ) * 1.e3;
    }
    

    
    Debug(fDebug,"muon weight: %, proton weight: %",MicroBooNEWeight_Pmu_theta,MicroBooNEWeight_Pp_theta);
    Debug(fDebug,"Q2 weight: %",MicroBooNEWeight_Q2);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::ProjectMuonTrajectory(){
    int fDebug=0;
    Debug(fDebug,"GenieFile::ProjectMuonTrajectory()");
    
    // the muon direction is given by its momentum...
    muonTrajectory_dir = muon.Vect().Unit();
    double Pmu_MeVc = 1000*muon.P();
    double dmu_cm = Get_muonRangeFromMomentum( Pmu_MeVc );
    muonTrajectory_end = vertex_position + muonTrajectory_dir * dmu_cm;
    
    
    
    Debug(fDebug,"Pmu: % MeV/c, lmu: % cm",Pmu_MeVc,lmu_cm);
    Debug(fDebug,"muon trajectory: (%,%,%) -> (%,%,%)"
          ,vertex_position.X(),vertex_position.Y(),vertex_position.Z()
          ,muonTrajectory_end.X(),muonTrajectory_end.Y(),muonTrajectory_end.Z());
}

double d_LinePlaneIntersection(TVector3 l,TVector3 l0,TVector3 p0,TVector3 n,float epsilon){
    // epsilon is a cutoff for the case of trajectory parallel to either one of the detector sides
    // [https://en.wikipedia.org/wiki/Line–plane_intersection]
    // where l is a vector in the direction of the line, l0 is a point on the line,
    // n is the normal to the plane and p0 is a point on the plan
    return ((p0 - l0).Dot(n))/(std::max(l.Dot(n),epsilon));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::MuonIntersectionWithSides(){
    // find the distance from the vertex position
    // in which the track meets the upper side of the detector (with normal in direction y
    // [https://en.wikipedia.org/wiki/Line–plane_intersection]
    // where l is a vector in the direction of the line, l0 is a point on the line,
    // n is the normal to the plane and p0 is a point on the plan
    TVector3 p0;
    // left and right sides
    if (muonTrajectory_dir.X() > 0){
        p0 = TVector3( side_x_right , 0 , 0 );
    }
    else if (muonTrajectory_dir.X() <= 0){
        p0 = TVector3( side_x_left , 0 , 0 );
    }
    dx = d_LinePlaneIntersection(muonTrajectory_dir,vertex_position,p0,TVector3( 1 , 0 , 0 ));
    
    // upper and lower sides
    if (muonTrajectory_dir.Y() > 0){
        p0 = TVector3( 0 , side_y_up , 0 );
    }
    else if (muonTrajectory_dir.Y() <= 0){
        p0 = TVector3( 0 , side_y_dw , 0 );
    }
    dy = d_LinePlaneIntersection(muonTrajectory_dir,vertex_position,p0,TVector3( 0 , 1 , 0 ));
    
    // up and downstream sides
    if (muonTrajectory_dir.Z() > 0){
        p0 = TVector3( side_z_upstream , 0 , 0 );
    }
    else if (muonTrajectory_dir.Z() <= 0){
        p0 = TVector3( side_z_downstream , 0 , 0 );
    }
    dz = d_LinePlaneIntersection(muonTrajectory_dir,vertex_position,p0,TVector3( 0 , 0 , 1 ));
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::CutMuonTrajectory(){
    int fDebug=0;
    
    if (muonTrajectory_dir.X()>0) {
        muonTrack_end_x = std::min( side_x_right , muonTrajectory_end.X() );
    }
    else {
        muonTrack_end_x = std::max( side_x_left , muonTrajectory_end.X() );
    }

    if (muonTrajectory_dir.Y()>0) {
        muonTrack_end_y = std::min( side_y_up , muonTrajectory_end.Y() );
    }
    else {
        muonTrack_end_y = std::max( side_y_dw , muonTrajectory_end.Y() );
    }

    if (muonTrajectory_dir.Z()>0) {
        muonTrack_end_z = std::min( side_z_upstream , muonTrajectory_end.Z() );
    }
    else {
        muonTrack_end_z = std::max( side_z_downstream , muonTrajectory_end.Z() );
    }
    
    muonTrack_end = TVector3( muonTrack_end_x , muonTrack_end_y , muonTrack_end_z );
    
    Debug(fDebug,"muon track: (%,%,%) -> (%,%,%)"
          ,vertex_position.X(),vertex_position.Y(),vertex_position.Z()
          ,muonTrack_end.X(),muonTrack_end.Y(),muonTrack_end.Z());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GenieFile::SetRecoKinematics(){
    reco_Pp  = proton;
    reco_Pmu = muon * (muonTrack_end.Mag()/muonTrajectory_end.Mag());

    reco_Ev  = reco_Pmu.E() + (reco_Pp.E() - reco_Pp.Mag()) + 0.040 ; // Eµ + Tp + Sn + T(A-1)
    reco_Pnu = TLorentzVector( 0 , 0 , reco_Ev , reco_Ev );

    reco_q   = reco_Pnu - reco_Pmu;
    reco_Q2  = - reco_q.Mag2();
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
    //    MuonIntersectionWithSides();
    CutMuonTrajectory();
    SetRecoKinematics();
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

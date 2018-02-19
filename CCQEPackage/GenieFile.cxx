#ifndef GENIEFILE_CXX
#define GENIEFILE_CXX

#include "GenieFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GenieFile::GenieFile( TString fPath
                     ,TString fRootFileName     // without the .root suffix (!)
                     ,TString fRootTreeName
                     ,int fdebug  ):
Path( fPath ),
RootFileName( fRootFileName + ".root"),
RootTreeName( fRootTreeName ),
OutputCSVname( fRootFileName + ".csv"),
debug(fdebug)
{
    SetInTree();
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
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::Initialize(){
    
    for (int i=0; i<NMAX ; i++ ){
        pdgf[i] = -9999;
        pxf[i] = pyf[i] = pzf[i] = -9999;
    }
    CC1p0pi = false;
    Pmiss = muon  = proton = q = nu = proton_before_FSI = TLorentzVector();
    protons.clear();
    neutrons.clear();
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
    << "CC1p0pi"    << ","
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
    << "Pp_before_FSI_x"       << "," << "Pp_before_FSI_y"        << "," << "Pp_before_FSI_z"
    << endl;
    
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
    << CC1p0pi      << ","
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
    << proton_before_FSI.Px()  << "," << proton_before_FSI.Py()       << "," << proton_before_FSI.Pz()
    << endl;
    

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::SetTopology (){
    nu.SetXYZM( pxv, pyv, pzv, 0 );
    if (debug>2) {
        SHOW4(cc,qel,nf,pdgf[0]);
        for (int i=0; i<nf; i++) {
            SHOW(pdgf[i]);
        }
    }
    
    
    
    int Nmu=0, Npi=0, Nel=0, Ngamma=0;
    int Np_200MeVc=0, Nn_200MeVc=0;
    
    muon.SetXYZM( pxl, pyl, pzl, 0.1056583745 );

    for( int i_f=0; i_f < nf ; i_f++ ){
        Debug(3,"pdgf[i_f]:%",pdgf[i_f]);
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
    if ((cc==true)
        && ( Np_200MeVc==1 )
        //        && ( Nmu>=1 ) // if cc==0 then we have a muon!
        && ( Npi==0 && Nel==0 && Ngamma==0 )
        && Nn_200MeVc==0
        //        && ((nf == 1) && (pdgf[0]==2212) )
        )
    {
        CC1p0pi = true;
        q = nu - muon;
        Pmiss = proton - q;
        // proton_before_FSI.SetXYZM( pxi[0],pyi[0],pzi[0], 0.9382720813 ); // correct this!
    }
     // CCQE with no FSI
    else if ((cc==true)
             && (nf == 2) && (pdgf[0]==1000180390) && (pdgf[1]==2212) ){
        CC1p0pi = true;
        muon.SetXYZM( pxl, pyl, pzl, 0.1056583745 );
        q = nu - muon;
        proton.SetXYZM( pxf[1], pyf[1], pzf[1], 0.9382720813 );
        Pmiss = proton - q;
        proton_before_FSI.SetXYZM( pxi[1],pyi[1],pzi[1], 0.9382720813 );
    }
    
    return true;
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
    return true;
}

#endif

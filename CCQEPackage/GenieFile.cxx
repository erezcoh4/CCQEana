#ifndef GENIEFILE_CXX
#define GENIEFILE_CXX

#include "GenieFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GenieFile::GenieFile( TString fPath
                     ,TString fRootFileName     // without the .root suffix (!)
                     ,TString fRootTreeName  ):
Path( fPath ),
RootFileName( fRootFileName + ".root"),
RootTreeName( fRootTreeName ),
OutputCSVname( fRootFileName + ".csv")
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
    Pmiss = muon  = proton = q = nu = TLorentzVector();
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
    << "Q2"  <<  "," << "W" <<  "," << "x" <<  ","
    << "pxn" <<  "," << "pyn" <<  "," << "pzn" <<  "," // Initial state hit nucleon
    << "pxv" <<  "," << "pyv" <<  "," << "pzv" <<  ","
    << "pxl" <<  "," << "pyl" <<  "," << "pzl" << ","
    << "nf"  <<  ",";
    
    for (int i=0; i<NMAX ; i++ ){
        csv_file
        << Form("pdgf_%d",i) <<  ","
        << Form("pxf_%d",i) <<  "," << Form("pyf_%d",i) <<  "," << Form("pzf_%d",i) << ",";
    }
    
    csv_file
    << "neu" << ","
    << "v_Px" << ","
    << "v_Py" << ","
    << "v_Pz" << ","
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
    << endl;
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::StreamToCSV (){
    
    counter++;
    
    csv_file
    << cc  << "," << qel <<  "," << res <<  "," << dis <<  "," << coh <<  ","
    << Q2  <<  "," << W <<  "," << x <<  ","
    << pxn <<  "," << pyn <<  "," << pzn <<  "," // Initial state hit nucleon
    << pxv <<  "," << pyv <<  "," << pzv <<  ","
    << pxl <<  "," << pyl <<  "," << pzl << ","
    << nf  <<  ",";
    
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
    << endl;
    

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool GenieFile::SetTopology (){
    if ((cc==true) && (qel==true) && (nf == 1) && pdgf[0]==2212 ){
        CC1p0pi = true;
        nu.SetXYZM( pxv, pyv, pzv, 0 );
        muon.SetXYZM( pxl, pyl, pzl, 0.1056583745 );
        q = nu - muon;
        proton.SetXYZM( pxf[0], pyf[0], pzf[0], 0.9382720813 );
        Pmiss = proton - q;
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

#ifndef ANALYSEEVENTS_CXX
#define ANALYSEEVENTS_CXX

#include "AnalyseEvents.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AnalyseEvents::AnalyseEvents (TTree * tree):
InTree( tree )
{
    InitInputTree();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnalyseEvents::InitInputTree(){
    
    InTree -> SetBranchAddress("Ntracks"        , &Ntracks);
    InTree -> SetBranchAddress("Nhits"          , &Nhits);
    
    Nentries = InTree -> GetEntries();
    if(debug>1) cout << "AnalyzeVertex input-tree ready (" << InTree -> GetName() <<"), " <<  Nentries << " entries" << endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnalyseEvents::InitEvent(){
    
    Debug ( 2 , "AnalyzeVertex::InitEvent");
    
    tracks.clear();
    hits.clear();
    
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnalyseEvents::GetEntry (int entry){
    
    InitEvent();
    
    std::vector <PandoraNuTrack> * ftracks = 0;
    InTree -> SetBranchAddress("tracks" , &ftracks);

    std::vector <hit> * fhits = 0;
    InTree -> SetBranchAddress("hits" , &fhits);
    

    std::vector <pairVertex> * fvertices = 0;
    InTree -> SetBranchAddress("vertices" , &fvertices);

    InTree -> GetEntry(entry);
    
    hits = *fhits;
    tracks = *ftracks;
    vertices = *fvertices;

    delete ftracks;
    delete fhits;
    delete fvertices;
}
#endif

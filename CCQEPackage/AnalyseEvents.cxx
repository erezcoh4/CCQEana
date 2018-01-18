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
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnalyseEvents::InitEvent(){
    
    tracks.clear();
    hits.clear();
    vertices.clear();
    
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnalyseEvents::GetEntry (int entry){
    
    InitEvent();
    
    std::vector <PandoraNuTrack> * ftracks = 0;
    InTree -> SetBranchAddress("tracks" , &ftracks);
    
    std::vector <PandoraNuTrack> * fcosmic_tracks = 0;
    InTree -> SetBranchAddress("cosmic_tracks" , &fcosmic_tracks);

    std::vector <hit> * fhits = 0;
    InTree -> SetBranchAddress("hits" , &fhits);
    

    std::vector <pairVertex> * fvertices = 0;
    InTree -> SetBranchAddress("vertices" , &fvertices);

    std::vector <pairVertex> * fcosmic_vertices = 0;
    InTree -> SetBranchAddress("cosmic_vertices" , &fcosmic_vertices);
    
    InTree -> GetEntry(entry);
    
    hits = *fhits;
    tracks = *ftracks;
    cosmic_tracks = *fcosmic_tracks;
    vertices = *fvertices;
    cosmic_vertices = *fcosmic_vertices;

    delete ftracks;
    delete fcosmic_tracks;
    delete fhits;
    delete fvertices;
    delete fcosmic_vertices;
}
#endif

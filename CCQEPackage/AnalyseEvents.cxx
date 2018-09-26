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
    // InTree -> SetBranchAddress("Nhits"          , &Nhits);
    
    Nentries = InTree -> GetEntries();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnalyseEvents::InitEvent(){
    
    tracks.clear();
    hits.clear();
    vertices.clear();
    tripleVertices.clear();
    mcparticles.clear();
    
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnalyseEvents::GetEntry (int entry){
    
    InitEvent();
    
    //    std::vector <PandoraNuTrack> * fmcparticles = 0;
    //    InTree -> SetBranchAddress("truth_trajectories" , &fmcparticles);
    
    std::vector <PandoraNuTrack> * ftracks = 0;
    InTree -> SetBranchAddress("tracks" , &ftracks);
    
    //    std::vector <PandoraNuTrack> * fcosmic_tracks = 0;
    //    InTree -> SetBranchAddress("cosmic_tracks" , &fcosmic_tracks);
    
    std::vector <hit> * fhits = 0;
    InTree -> SetBranchAddress("hits" , &fhits);
    
    
    
    std::vector <pairVertex> * fvertices = 0;
    InTree -> SetBranchAddress("vertices" , &fvertices);
    
    //    std::vector <pairVertex> * fcosmic_vertices = 0;
    //    InTree -> SetBranchAddress("cosmic_vertices" , &fcosmic_vertices);
    
    InTree -> GetEntry(entry);
    
    hits = *fhits;
    //    mcparticles = *fmcparticles;
    tracks = *ftracks;
    //    cosmic_tracks = *fcosmic_tracks;
    vertices = *fvertices;
    //    cosmic_vertices = *fcosmic_vertices;
    
    delete ftracks;
    //    delete fmcparticles;
    //    delete fcosmic_tracks;
    delete fhits;
    delete fvertices;
    //    delete fcosmic_vertices;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnalyseEvents::GetEntryOnlyTracks (int entry){
    
    InitEvent();
    cout << "InitEvent()" << endl;
    std::vector <PandoraNuTrack> * ftracks = 0;
    InTree -> SetBranchAddress("tracks" , &ftracks);
    cout << "SetBranchAddress(tracks , &ftracks)" << endl;
    
    std::vector <tripleVertex> * ftripleVertices = 0;
    InTree -> SetBranchAddress("vertices" , &ftripleVertices);
    cout << "SetBranchAddress(vertices , &ftripleVertices)" << endl;
    
    InTree -> GetEntry(entry);
    cout << "Got Entry("<<entry<<")" << endl;
    
    tracks = *ftracks;
    cout << "tracks = *ftracks;" << endl;
    
    tripleVertices = *ftripleVertices;
    cout << "tripleVertices = *ftripleVertices;" << endl;

    delete ftracks;
    cout << "delete ftracks;" << endl;
    delete ftripleVertices;
    cout << "delete ftripleVertices;" << endl;
}

#endif

/**
 * \file AnalyseEvents.h
 *
 * \ingroup CCQEPackage
 *
 * \brief Class def header for a class AnalyseEvents
 *
 * @author erezcohen
 */

/** \addtogroup CCQEPackage
 
 @{*/
#ifndef ANALYSEEVENTS_H
#define ANALYSEEVENTS_H

#include <iostream>
//#include "../../mySoftware/MySoftwarePackage/myIncludes.h"
//#include "../../AnalysisTreesInformation/AnaTreesPackage/PandoraNuTrack.h"
//#include "../../AnalysisTreesInformation/AnaTreesPackage/hit.h"
#include "TTree.h"
#include "Rtypes.h"
#include "PandoraNuTrack.h"
#include "hit.h"
#include "box.h"
#include "GENIEinteraction.h"
#include "pairVertex.h"

/**
 \class AnalyseEvents
 User defined class AnalyseEvents ... these comments are used to generate
 doxygen documentation!
 */
class AnalyseEvents{
    
public:
    
    /// Default constructor
    AnalyseEvents(){}
    ~AnalyseEvents(){}
    AnalyseEvents (TTree * tree);
    
    
    // SETters
    void                  SetInTree (TTree * tree)       {InTree = tree;};
    
    // GETters
    void                           GetEntry ( int );
    std::vector<hit>                GetHits ()  const {return hits;};
    std::vector<PandoraNuTrack>   GetTracks ()  const {return tracks;};
    std::vector<pairVertex>     GetVertices ()  const {return vertices;};
    
    
    // INITializers
    void                  InitEvent ();
    void              InitInputTree ();
    
    

    
    
private:
    
    Int_t   Nentries, Nhits , Ntracks;
    
    TTree   * InTree;
    
    std::vector<PandoraNuTrack> tracks;
    
    std::vector<hit>            hits;

    std::vector<pairVertex>     vertices;

};

#endif
/** @} */ // end of doxygen group


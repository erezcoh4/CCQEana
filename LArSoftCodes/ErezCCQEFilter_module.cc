////////////////////////////////////////////////////////////////////////
// Class:       ErezCCQEFilter
// Plugin Type: filter (art v2_05_01)
// File:        ErezCCQEFilter_module.cc
//
// Generated at Tue Jun 19 10:46:33 2018 by Erez Cohen using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

namespace ub {
    class ErezCCQEFilter;
}


class ub::ErezCCQEFilter : public art::EDFilter {
public:
    explicit ErezCCQEFilter(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    ErezCCQEFilter(ErezCCQEFilter const &) = delete;
    ErezCCQEFilter(ErezCCQEFilter &&) = delete;
    ErezCCQEFilter & operator = (ErezCCQEFilter const &) = delete;
    ErezCCQEFilter & operator = (ErezCCQEFilter &&) = delete;
    
    // Required functions.
    bool filter(art::Event & e) override;
    
    // Selected optional functions.
    void beginJob() override;
    
private:
    
    // Declare member data here.
    
};


ub::ErezCCQEFilter::ErezCCQEFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
    // Call appropriate produces<>() functions here.
}

bool ub::ErezCCQEFilter::filter(art::Event & evt)
{
    int run = evt.run();
    //    int subrun = evt.subRun();
    int event = evt.id().event();
    std::cout << "run: "<< run << ", event: "<< event << std::endl;
    
    if (run==7004 && (event == 21860 || event==21859|| event==21857)) return true;

    return false;
}

void ub::ErezCCQEFilter::beginJob()
{
    // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ub::ErezCCQEFilter)

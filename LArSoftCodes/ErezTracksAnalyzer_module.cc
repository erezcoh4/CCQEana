////////////////////////////////////////////////////////////////////////
// Class:       ErezTracksAnalyzer
// Plugin Type: analyzer (art v2_07_03)
// File:        ErezTracksAnalyzer_module.cc
//
// Generated at Fri Oct 13 07:54:46 2017 by Erez Cohen using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace ub {
  class ErezTracksAnalyzer;
}


class ub::ErezTracksAnalyzer : public art::EDAnalyzer {
public:
  explicit ErezTracksAnalyzer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ErezTracksAnalyzer(ErezTracksAnalyzer const &) = delete;
  ErezTracksAnalyzer(ErezTracksAnalyzer &&) = delete;
  ErezTracksAnalyzer & operator = (ErezTracksAnalyzer const &) = delete;
  ErezTracksAnalyzer & operator = (ErezTracksAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.

};


ub::ErezTracksAnalyzer::ErezTracksAnalyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void ub::ErezTracksAnalyzer::analyze(art::Event const & e)
{
  // Implementation of required member function here.
}

void ub::ErezTracksAnalyzer::beginJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ub::ErezTracksAnalyzer)

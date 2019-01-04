#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "PU14/EventMixer.hh"
#include "PU14/CmdLine.hh"
#include "PU14/PU14.hh"

#include "include/rhoEstimator.hh"

#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"

#include "include/Angularity.hh"

using namespace std;
using namespace fastjet;

// ASO local: ./runRhoEstimator -hard /Users/albasotoontoso/Work/Jet_substraction/JetToyHI/samples/PythiaEventsTune14PtHat120_0.pu14 -pileup /Users/albasotoontoso/Work/Jet_substraction/JetToyHI/samples/ThermalEventsMult7000PtAv1.20_0.pu14 -nev 10

int main (int argc, char ** argv) {


  auto start_time = std::chrono::steady_clock::now();

  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  std::cout << "will run on " << nEvent << " events" << std::endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  TFile *fout = new TFile("JetToyHIResultRhoEstimator_test.root","RECREATE");
  fout->cd();
  treeWriter trw("jetTree");

  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  fastjet::JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 3.0;
  fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( mixer.next_event() && iev < nEvent )
  {
    // increment event number
    iev++;

    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    std::vector<fastjet::PseudoJet> particlesMerged = mixer.particles();

    std::vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());
    // cluster hard event only
    std::vector<fastjet::PseudoJet> particlesBkg, particlesSig;
    SelectorIsHard().sift(particlesMerged, particlesSig, particlesBkg); // this sifts the full event into two vectors of PseudoJet, one for the hard event, one for the underlying event



    // run the clustering, extract the bkg jets
    fastjet::JetDefinition jet_def_bkgd(fastjet::antikt_algorithm, 0.4);
    fastjet::AreaDefinition area_def_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
    fastjet::ClusterSequenceArea csBkg(particlesBkg, jet_def_bkgd, area_def_bkgd);
    fastjet::Selector bkg_selector = fastjet::SelectorAbsRapMax(jetRapMax);
    jetCollection jetCollectionBkg(sorted_by_pt(bkg_selector(csBkg.inclusive_jets())));

     vector<double> rho_True;
     rho_True.reserve(jetCollectionBkg.getJet().size());
     vector<double> nConstituents_True;
     nConstituents_True.reserve(jetCollectionBkg.getJet().size());

    for(PseudoJet jet : jetCollectionBkg.getJet()) {
       rho_True.push_back(jet.pt()/jet.area());
       std::vector<fastjet::PseudoJet> particles, ghosts;
       fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
      nConstituents_True.push_back(particles.size());
    }
    jetCollectionBkg.addVector("rho_Truth", rho_True);
    jetCollectionBkg.addVector("nconstituents_Truth", nConstituents_True);

    //---------------------------------------------------------------------------
    //   background Estimation with Rho Estimator
    //---------------------------------------------------------------------------

    rhoEstimator rhoComputation(R,0.001,ghostRapMax,jetRapMax);
    rhoComputation.setInputParticles(particlesMerged);
    vector<double> rho_Estimate(rhoComputation.doEstimation());

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------

    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported


    trw.addCollection("rhoTruth",           rho_True);
    trw.addCollection("nconstituentsTruth",  nConstituents_True);
    trw.addCollection("rhoEstimated",      rho_Estimate);


    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();
  cout << iev << endl;
  fout->cd();

 TTree *trOut = trw.getTree();
 trOut->Write();
 trOut->Scan();

 fout->Write();
 fout->Close();

  std::cout << "Check JetToyHIResultRhoEstimator_test.root for results" << std::endl;

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}

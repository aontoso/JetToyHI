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

#include "include/jetCollection.hh"
#include "include/csSubtractor.hh"
#include "include/csSubtractorFullEvent.hh"
#include "include/skSubtractor.hh"
#include "include/softDropGroomer.hh"
#include "include/softDropCounter.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/randomCones.hh"
#include "include/Angularity.hh"
#include "include/jewelMatcher.hh"
#include "include/gridSubtractor.hh"

using namespace std;
using namespace fastjet;

// Pb+Pb@5.02 TeV w/o Recoils

// ./runSDGenVarious -hard  /Users/albasotoontoso/Work/Jet_substraction/JetToyHI/samples/PbPb5p02TeV_dijet_cent05_woRecoils_50pthatInf_default_0.pu14 -nev 1

// p+p@5.02 TeV
// ./runSDGenVarious -hard  /Users/albasotoontoso/Work/Jet_substraction/JetToyHI/samples/pp5p02TeV_dijet_50pthatInf_1.pu14 -nev 1

// Pb+Pb@5.02 TeV w Recoils

// ./runSDGenVarious -hard  /Users/albasotoontoso/Work/Jet_substraction/JetToyHI/samples/PbPb5p02TeV_dijet_cent05_wRecoils_50pthatInf_default_0.pu14 -nev 1


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

  Angularity width(1.,1.,R);
  Angularity pTD(0.,2.,R);

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
    SelectorIsHard().sift(particlesMerged, particlesSig, particlesBkg); //  this sifts the full event into two vectors of PseudoJet, one for the hard event, one for the underlying event
    std::vector<fastjet::PseudoJet> scatteringCenters;

    // Find the scatteringCenters (i.e. recoil partons)

    for(fastjet::PseudoJet p : particlesMerged) {
         if (p.user_info<PU14>().vertex_number() == -1){ scatteringCenters.push_back(p);
         }
      }
    //  cout << scatteringCenters.size() << endl;
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    //--------------------------------------------------------------------------
    //-
    // -------------------------------------------------------------------------
    //- GridSubtraction
    //--------------------------------------------------------------------------

    gridSubtractor jewelSub(-1., -M_PI, 1., M_PI,0.4,0.05);
  //  std::vector<fastjet::PseudoJet> prueba = jewelSub.doGridSub1(jetCollectionSig,scatteringCenters);
    std::vector<fastjet::PseudoJet> subtracted_particles = jewelSub.doGridSub1(jetCollectionSig,scatteringCenters);
    //cout << subtracted_particles.size() << endl;
    fastjet::ClusterSequenceArea csSub(subtracted_particles, jet_def, area_def);
    jetCollection jetCollectionSub(sorted_by_pt(jet_selector(csSub.inclusive_jets(10.))));


    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------

    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
  //  softDropGroomer sdgSigBeta00Z01(0.1, 0, R);
    //jetCollection jetCollectionSigSD(sdgSigBeta00Z01.doGrooming(jetCollectionSub));
  //  jetCollectionSigSD.addVector("zgSigSDBeta00Z01",    sdgSigBeta00Z01.getZgs());
    //jetCollectionSigSD.addVector("ndropSigSDBeta00Z01", sdgSigBeta00Z01.getNDroppedSubjets());
  //  jetCollectionSigSD.addVector("dr12SigSDBeta00Z01",  sdgSigBeta00Z01.getDR12());

    //--------------------------------------------------------------------------
    // Count the number of soft drop splittings
    //--------------------------------------------------------------------------

    softDropCounter sdgCounter(0.1,0.,R,0.1);
    sdgCounter.run(jetCollectionSub);
    jetCollectionSub.addVector("nsd", sdgCounter.calculateNSD(0.,0.));
    jetCollectionSub.addVector("tf", sdgCounter.calculateTf());
    jetCollectionSub.addVector("deltaR", sdgCounter.calculateDeltaR());
    jetCollectionSub.addVector("zg", sdgCounter.calculateZg());
    jetCollectionSub.addVector("kt", sdgCounter.calculatekT());
  //  std::vector<std::vector<double>> Result = sdgCounter.calculateTf();
  //  if (nEvent==1)cout << Result.at(0).size() << endl;

//jetCollectionSigSD.addDoubleCollection("bla",Result);
    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------

    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trw.addCollection("eventWeight",   eventWeight);
    trw.addCollection("subJet",        jetCollectionSub, true);
  //  trw.addCollection("sdJet",      jetCollectionSigSD, true);
    //trw.addDoubleVectorCollection("tf", Result);
//

    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile("JetToyHIResultSD_JEWELwRecoil_PbPb_Zcut01Rcut01_test.root","RECREATE");
  trOut->Write();
//  trOut->Scan();
  fout->Write();
  fout->Close();
  //trOut->GetEntry(0);
    std::cout << "Check JetToyHIResultSD_JEWELwRecoil_PbPb_Zcut01Rcut01_test.root for results" << std::endl;

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}

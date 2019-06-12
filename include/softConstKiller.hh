#ifndef softConstKiller_h
#define softConstKiller_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <map>
#include <random>
#include <numeric>
#include "TH1.h"
#include "TROOT.h"


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "skEstimator.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"

#include "../PU14/PU14.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs a jet-by-jet constituent subtraction with SoftKiller
// background estimator
// Author: A. Soto-Ontoso, M. Verweij, Y. Mehtar-Tani
//---------------------------------------------------------------

class softConstKiller {

private :
  double jetRParam_;
  double rKTParam_;
  double ghostArea_;
  double ghostRapMax_;
  double jetRapMax_;


  std::vector<fastjet::PseudoJet> fjInputs_;
  std::vector<fastjet::PseudoJet> fjJetInputs_;

  std::vector<fastjet::PseudoJet> fjJetParticles_;
  std::vector<std::vector<fastjet::PseudoJet>> Hard;
  std::vector<std::vector<fastjet::PseudoJet>> Soft;

  std::random_device rd_;

  fastjet::ClusterSequenceArea *csJetSub;

public :
   softConstKiller(double rJet = 0.4,
                        double rKT = 0.19,
                        double ghostArea = 0.001,
                        double ghostRapMax = 3.0,
                        double jetRapMax = 3.0) :
    jetRParam_(rJet),
    rKTParam_(rKT),
    ghostArea_(ghostArea),
    ghostRapMax_(ghostRapMax),
    jetRapMax_(jetRapMax)
    {

    }

  void setGhostArea(double a) { ghostArea_ = a; }

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }
  void setInputJets(std::vector<fastjet::PseudoJet> v)      { fjJetInputs_ = v; }

  void clearMemory() {

  if(csJetSub) delete csJetSub;
}

std::vector<fastjet::PseudoJet> doSubtraction() {

  Hard.clear();
  Soft.clear();

fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, ghostArea_);
fastjet::GhostedAreaSpec ghost_spec_sub(ghostRapMax_, 1, 1.);
std::vector<fastjet::PseudoJet> jets = fjJetInputs_;

// do the clustering with ghosts and get the jets
//----------------------------------------------------------
fastjet::JetDefinition jet_def(antikt_algorithm, jetRParam_);
fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec);

fastjet::ClusterSequenceArea cs(fjInputs_, jet_def, area_def);
fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax_);
jets = fastjet::sorted_by_pt(jet_selector(cs.inclusive_jets()));

//cout << "SoftConstKiller: " << jets.at(0).pt()<< endl;

fastjet::JetDefinition jet_defSub(antikt_algorithm, 20.);
fastjet::AreaDefinition area_def_sub = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec_sub);
fastjet::ClusterSequenceArea csSub(fjInputs_, jet_defSub, area_def_sub);

// SoftKiller estimate for the Background pT
// -----------------------------------------------------------

skEstimator rhoComputation(jetRParam_,rKTParam_,ghostArea_,ghostRapMax_,jetRapMax_);
rhoComputation.setInputParticles(fjInputs_);
vector<double> rho_Estimate(rhoComputation.doEstimation());


  // Jet-by-jet subtraction
  //----------------------------------------------------------
  std::vector<fastjet::PseudoJet> csjets;
  csjets.reserve(jets.size());
  int ijet = -1;

  for(fastjet::PseudoJet& jet : jets) {
      ijet++;
      if(jet.is_pure_ghost()) continue;
    // Constituent Subtractor JetDefinition
    //------------------------------------------------------------

     std::vector<fastjet::PseudoJet> particles_two, ghosts_two;
     fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts_two, particles_two);
     if(particles_two.size()<1 || jet.pt()<1.) continue;


      contrib::ConstituentSubtractor subtractor_(rho_Estimate[ijet]/jet.area(), rho_Estimate[ijet]/jet.area());
      subtractor_.set_distance_type(contrib::ConstituentSubtractor::deltaR);
      subtractor_.set_max_distance(-1.); //free parameter for the maximal allowed distance between particle i and ghost k
      subtractor_.set_alpha(1.); // free parameter for the distance measure (the exponent of particle pt). Note that in older versions of the package alpha was multiplied by two but in newer versions this is not the case anymore

    fastjet::PseudoJet subtracted_jet = subtractor_(jet);
    std::vector<fastjet::PseudoJet> particles, ghosts;
    fastjet::SelectorIsPureGhost().sift(subtracted_jet.constituents(), ghosts, particles);

    std::vector<fastjet::PseudoJet> A, B;
    SelectorIsHard().sift(jet.constituents(), A, B);

     std::vector<fastjet::PseudoJet> hard, soft;
     for(fastjet::PseudoJet p : subtracted_jet.constituents())
     {
        double BestA = -1, BestB = -1;
        for(fastjet::PseudoJet x : A)
        {
           double DR2 = p.squared_distance(x);
           if(BestA < 0 || BestA > DR2)
              BestA = DR2;
        }
        for(fastjet::PseudoJet x : B)
        {
           double DR2 = p.squared_distance(x);
           if(BestB < 0 || BestB > DR2)
              BestB = DR2;
        }
        if(BestA < BestB)
           hard.push_back(p);
        else
           soft.push_back(p);
     }

     if(particles.size() > 0)
     {
        std::vector<fastjet::PseudoJet> combinedparticles;
        for(fastjet::PseudoJet p : particles)
           combinedparticles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E()));

          for(fastjet::PseudoJet p : jet.constituents())
           if(p.E() < 1e-5)
              combinedparticles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E()));

        csjets.push_back(fastjet::PseudoJet(join(combinedparticles)));
        Hard.push_back(hard);
        Soft.push_back(soft);
     }
     //cout << csjets.at(0).pt() << endl;
  }

//  std::cout << "\n n subtracted jets: " << csjets.size() << std::endl;

  return csjets;

} // end of doSubtraction


};

#endif
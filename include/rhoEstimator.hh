#ifndef rhoEstimator_h
#define rhoEstimator_h

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
#include "TGraph.h"
#include "TF1.h"


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "../PU14/PU14.hh"

#include "Angularity.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the jet-by-jet shared layer subtraction
// Author: M. Verweij, Y. Mehtar-Tani
//---------------------------------------------------------------

class rhoEstimator{

private :
  double jetRParam_;
  double ghostArea_;
  double ghostRapMax_;
  double jetRapMax_;
  double rho_;
  double rhoSigma_;

  std::vector<fastjet::PseudoJet> fjInputs_;
  std::vector<fastjet::PseudoJet> fjJetInputs_;

//  std::random_device rd_;

  //fastjet::ClusterSequenceArea *csJetSub;

public :
  rhoEstimator(double rJet = 0.4,
                        double ghostArea = 0.001,
                        double ghostRapMax = 3.0,
                        double jetRapMax = 3.0
                      ):
    jetRParam_(rJet),
    ghostArea_(ghostArea),
    ghostRapMax_(ghostRapMax),
    jetRapMax_(jetRapMax)

{

}
  void setGhostArea(double a) { ghostArea_ = a; }

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }
  void setInputJets(std::vector<fastjet::PseudoJet> v)      { fjJetInputs_ = v; }

  double getRho()  const { return rho_; }
  double getRhoSigma() const { return rhoSigma_; }

  // By calling this function we can get the value of the "real" pt given the linear assumption
   double correlation_function(double slope, int n_constituents){
   return slope*n_constituents;
   };

  std::vector<double> doEstimation() {

    fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, ghostArea_);
    fastjet::GhostedAreaSpec ghost_spec_sub(ghostRapMax_, 1, 1.);
    std::vector<fastjet::PseudoJet> jets = fjJetInputs_;

    // Create what we need for the background estimation
    //----------------------------------------------------------
    fastjet::JetDefinition jet_def_bkgd(fastjet::kt_algorithm, jetRParam_);
    fastjet::AreaDefinition area_def_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
    fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4) * (!fastjet::SelectorNHardest(2));

    fastjet::ClusterSequenceArea csKt(fjInputs_, jet_def_bkgd, area_def_bkgd);
    std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

    fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    bkgd_estimator.set_particles(fjInputs_);
    bkgd_estimator.set_jets(bkgd_jets);

    double rhoMedian_ = 0;
    double rhoSigma_ = 0;

    // Compute the rho median for each event and its standard deviation

    rhoMedian_ = bkgd_estimator.rho();
    rhoSigma_ = bkgd_estimator.sigma();

    // Determination of the slope ("temperature") of the correlation line
    //-----------------------------------------------------------

    double temperature = 1.2;

  //  for(fastjet::PseudoJet& jet : bkgd_jets){


  //  }

    // Do the anti-kt clustering
    //------------------------------------------------------------

    fastjet::JetDefinition jet_def(antikt_algorithm, jetRParam_);
    fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec);

    fastjet::ClusterSequenceArea cs(fjInputs_, jet_def, area_def);
    fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax_);
    jets = fastjet::sorted_by_pt(jet_selector(cs.inclusive_jets()));

    std::vector<double> rho_estimate; // one for each patch;
    rho_estimate.reserve(jets.size());

    int ijet = -1;

    for(fastjet::PseudoJet& jet : jets) {
      ++ijet;
      if(jet.is_pure_ghost()) continue;

      //Get ghosts and true particles
      //----------------------------------------------------------
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
      if(particles.size()<1 || jet.pt()<1.) continue;

      // Order particles from soft to high pT

      std::vector<size_t> idx(particles.size());
      iota(idx.begin(), idx.end(), 0);

      std::sort(idx.begin(), idx.end(),
                [&particles](size_t i1, size_t i2) {return particles[i1].pt() < particles[i2].pt();});
      std::vector<fastjet::PseudoJet> particles_pTOrdered;
      for (int i = 0; i<particles.size(); i++){
         particles_pTOrdered.push_back(particles[idx[i]]);

      }

      // Set user_index of all particles to position particles vector

       for(int i = 0; i<(int)particles_pTOrdered.size(); ++i) {
         particles_pTOrdered[i].set_user_index(i);
       }

      std::vector<fastjet::PseudoJet> bkgEstimate;  //list of particles until we reach the correlation line
      std::vector<int> avail_part(particles_pTOrdered.size()); // list of particles still available
      double pTcurrent = 0; // to make sure it goes into the while loop
      int nparticles = -1;

      double pTcorrelation = 1;

      std::fill(avail_part.begin(),avail_part.end(),1);

      while(pTcurrent<pTcorrelation && std::accumulate(avail_part.begin(),avail_part.end(),0)>0){

      nparticles++;
      pTcorrelation = correlation_function(temperature, nparticles+1);
      fastjet::PseudoJet partSel = particles_pTOrdered[nparticles];
      bkgEstimate.push_back(partSel);
      pTcurrent+= partSel.pt();
      avail_part[nparticles] = 0;
    //  if (ijet ==0) cout <<pTcurrent << endl;
      }
      rho_estimate.push_back(pTcurrent);
    //  cout << "Estimate: " << pTcurrent << "Rho median: " << rhoMedian_*jet.area()<< endl;
    } // jets_loop

    return rho_estimate;


};

};

#endif

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
#include "TCanvas.h"


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "skSubtractor.hh"

#include "../PU14/PU14.hh"

#include "Angularity.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the jet-by-jet rho Estimator
// Author: M. Verweij, Y. Mehtar-Tani, Alba Soto-Ontoso
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
   double correlation_line(double slope, int n_constituents){
   return slope*(n_constituents+1); //+1 because the arrays start in zero
   };

  std::vector<double> doEstimation() {

  //  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, ghostArea_);
    std::vector<fastjet::PseudoJet> jets = fjJetInputs_;

    // Create what we need for the background estimation
    //---------------------------------------------------------

   fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, 0.001);

   fastjet::JetDefinition jet_estimate_bkgd(fastjet::kt_algorithm, 0.4);
   fastjet::AreaDefinition area_estimate_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
   fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4)*(!fastjet::SelectorNHardest(2));

   fastjet::ClusterSequenceArea csKt(fjInputs_, jet_estimate_bkgd, area_estimate_bkgd);
   std::vector<fastjet::PseudoJet> bkg_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

    fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_estimate_bkgd, area_estimate_bkgd);

    bkgd_estimator.set_particles(fjInputs_);
    bkgd_estimator.set_jets(bkg_jets);

  //  fastjet::GridMedianBackgroundEstimator grid_estimator(3., 0.7);
  //  grid_estimator.set_particles(fjInputs_);
  //  grid_estimator.set_jets(bkg_jets);

  //  cout << grid_estimator.rho() << endl;

    double rhoMedian_ = 0;
    double rhoCut_ = 0;
    double rhoGammaCut_ = 0;
    double lambda = 1;
    // Compute the rho median for each event and its standard deviation

     rhoMedian_ = bkgd_estimator.rho();
//    rhoSigma_ = bkgd_estimator.sigma();
//cout << rhoMedian_ << endl;

    // Determine the ptCut a la Soft Killer
    //-------------------------------------------------------------

    std::vector<double> ptConst_bkgd;
    std::vector<double> gammaConst_bkgd;

    for(fastjet::PseudoJet& jet : bkg_jets) {
        double maxpt_bkgd = 0;
        double maxgamma_bkgd = 0;
        if(jet.is_pure_ghost()) continue;
        std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
        fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
      //  cout << ghosts.size() << endl;
     // Find the maximum-pt of patch in the background
       for(fastjet::PseudoJet p : particles) {
            double deltaR = jet.delta_R(p);
            if (deltaR == 0) deltaR=1e-8;
            double momentum = p.pt();
            if (momentum > maxpt_bkgd) maxpt_bkgd = momentum;
            double gamma = momentum/pow(deltaR,lambda);
            if (gamma > maxgamma_bkgd) maxgamma_bkgd = gamma;
          }
       ptConst_bkgd.push_back(maxpt_bkgd);
       gammaConst_bkgd.push_back(maxgamma_bkgd);
      //    cout << maxpt_bkgd << endl;
       } // bkgd_jets loop

       std::nth_element(ptConst_bkgd.begin(), ptConst_bkgd.begin() + ptConst_bkgd.size()/2, ptConst_bkgd.end());
       double pTsoftKiller = ptConst_bkgd[ptConst_bkgd.size()/2];
      // cout << pTsoftKiller << endl;
        double alpha = 1.; // alpha is between (0,1). Controls how much we want to be like soft Killer. For alpha = 0 we recover the area-median result. For alpha = 1, soft Killer.
       double pTModifiedsoftKiller = alpha * pTsoftKiller;

       std::nth_element(gammaConst_bkgd.begin(), gammaConst_bkgd.begin() + gammaConst_bkgd.size()/2, gammaConst_bkgd.end());

       double gammaCut = alpha*gammaConst_bkgd[gammaConst_bkgd.size()/2];
    //   double gammaCut = 2.4;
      //  cout << gammaCut << endl;
    //   double pTModifiedsoftKiller = 2.4;
       // Determine the rho(cut): rho computed after pT < pTcut have been removed
       //-------------------------------------------------------------

    //   skSubtractor skSub(0.4, 3.0);
    //   skSub.setInputParticles(fjInputs_);
    //   std::vector<fastjet::PseudoJet> skEvent = skSub.doSubtraction();
  //     std::vector<double> skPtThreshold;
  //     skPtThreshold.push_back(skSub.getPtThreshold()); //SoftKiller pT ///threshold

      // cout << skPtThreshold.at(0) << endl;

  //     double pTModifiedsoftKiller = skPtThreshold.at(0);

       std::vector<double> rhoMedian_bkgd;
       std::vector<double> rhoCut_bkgd;
       std::vector<double> rhoGammaCut_bkgd;

         for(fastjet::PseudoJet& jet : bkg_jets) {
             double rho_const = 0;
             double rho_cut = 0;
             double rho_gamma_cut = 0;
             std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
             fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
          // Find the maximum-pt of patch in the background
            for(fastjet::PseudoJet p : particles) {
                 double momentum = p.pt();
                 rho_const+=momentum;

                 if (momentum > pTModifiedsoftKiller){
                 rho_cut+=momentum;
                 }
                 double deltaR = jet.delta_R(p);
                 if (deltaR == 0) deltaR=1e-8;
                 double gamma = momentum/pow(deltaR,lambda);

                 if (gamma > gammaCut){
                  rho_gamma_cut+=momentum;
                 }
               }
              rhoMedian_bkgd.push_back(rho_const/jet.area());
              rhoCut_bkgd.push_back(rho_cut/jet.area());
              rhoGammaCut_bkgd.push_back(rho_gamma_cut/jet.area());
            //  cout <<"hh:" << jet.area() << endl;
            } // bkgd_jets loop

            std::nth_element(rhoMedian_bkgd.begin(), rhoMedian_bkgd.begin() + rhoMedian_bkgd.size()/2, rhoMedian_bkgd.end());

            std::nth_element(rhoCut_bkgd.begin(), rhoCut_bkgd.begin() + rhoCut_bkgd.size()/2, rhoCut_bkgd.end());

            std::nth_element(rhoGammaCut_bkgd.begin(), rhoGammaCut_bkgd.begin() + rhoGammaCut_bkgd.size()/2, rhoGammaCut_bkgd.end());

           rhoMedian_ = rhoMedian_bkgd[rhoMedian_bkgd.size()/2];
            //   cout << rhoCut_bkgd.size() << endl;
           rhoCut_ = rhoCut_bkgd[rhoCut_bkgd.size()/2];

           rhoGammaCut_ = rhoGammaCut_bkgd[rhoGammaCut_bkgd.size()/2];


    //   cout << pTsoftKiller << endl;

    // Determination of the slope ("temperature") of the correlation line
    //----------------------------------------------------------

     double temperature = 1.2;

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

      // Raw values
      std::vector<double> particles_deltaR;
      std::vector<double> particles_gamma;
      double truePatchPt = 0;
      int truePatchN = 0;
      double pTsignal = 0;
    //  double trueRho = 0;
      for(fastjet::PseudoJet p : particles) {
      // Obtain the true rho for each patch by using only the background particles
       truePatchPt+=p.pt();
       truePatchN++;
       double radial = jet.delta_R(p);
       if (radial == 0) radial = 1e-8;
       particles_deltaR.push_back(radial);
       }

      // Order particles from soft to high pT

      std::vector<size_t> idx(particles.size());
      iota(idx.begin(), idx.end(), 0);

      std::sort(idx.begin(), idx.end(),
                [&particles](size_t i1, size_t i2) {return particles[i1].pt() < particles[i2].pt();});
      std::vector<fastjet::PseudoJet> particles_pTOrdered;
      for (int i = 0; i<particles.size(); i++){
         particles_pTOrdered.push_back(particles[idx[i]]);
      }

      // Order particles from low to high values of Gamma

      std::vector<size_t> idx2(particles.size());
      iota(idx2.begin(), idx2.end(), 0);

      std::sort(idx2.begin(), idx2.end(),
                [&particles, &particles_deltaR](size_t i1, size_t i2, double lambda = 0.5) {
                return particles[i1].pt()/pow(particles_deltaR[i1],lambda) < particles[i2].pt()/pow(particles_deltaR[i2],lambda);});

      std::vector<fastjet::PseudoJet> particles_gammaOrdered;
      for (int i = 0; i<particles.size(); i++){
         particles_gammaOrdered.push_back(particles[idx2[i]]);
      }

      // Set user_index of all particles to position particles vector

       for(int i = 0; i<(int)particles_pTOrdered.size(); ++i) {
         particles_pTOrdered[i].set_user_index(i);
       }

       for(int i = 0; i<(int)particles_gammaOrdered.size(); ++i) {
         particles_gammaOrdered[i].set_user_index(i);
       }

       for(int i = 0; i<(int)particles_gammaOrdered.size(); ++i) {
         double gamma_value = particles_gammaOrdered[i].pt()/pow(jet.delta_R(particles_gammaOrdered[i]),lambda);
         particles_gamma.push_back(gamma_value);
       }

    //   cout << particles_gammaOrdered[0] << endl;

      std::vector<fastjet::PseudoJet> bkgEstimate;  //list of particles until we reach the correlation line
      std::vector<int> avail_part(particles_pTOrdered.size()); // list of particles still available

      double pTcurrent = 0; // to make sure it goes into the while loop
      int nparticles = -1;

      std::fill(avail_part.begin(),avail_part.end(),1);

    //  while(std::accumulate(avail_part.begin(),avail_part.end(),0)>0 && particles_pTOrdered[nparticles+1].pt()<pTModifiedsoftKiller){
      while(std::accumulate(avail_part.begin(),avail_part.end(),0)>0 && particles_gamma[nparticles+1]<gammaCut){

      nparticles++;
      fastjet::PseudoJet partSel = particles_gammaOrdered[nparticles];
      bkgEstimate.push_back(partSel);
      pTcurrent+= partSel.pt();
      avail_part[nparticles] = 0;
      //if (partSel.user_info<PU14>().vertex_number() == 1){
      //  pTsignal+=partSel.pt();
    //  if (ijet==0) cout << partSel.pt()*pow(jet.delta_R(partSel),2) << endl;
    //
      }
  //      cout << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;

      //Trial to correct for the signal contramination below the pt cut
      //------------------------------------------------------------------

    //    double pT_cut_corrected = nparticles*temperature*(1-(1+2*pTModifiedsoftKiller/temperature+2*pow(pTModifiedsoftKiller,2)/pow(temperature,2))*exp(-2*pTModifiedsoftKiller/temperature));
      //  if (pTcurrent > pT_cut_corrected) pTcurrent = pT_cut_corrected;
        double pT_estimate = pTcurrent+rhoGammaCut_*jet.area();
       if (pT_estimate > truePatchPt) pT_estimate = truePatchPt;

      rho_estimate.push_back(pT_estimate);

  }

    return rho_estimate;


};

};

#endif

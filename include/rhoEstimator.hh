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
  // fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4);

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
    // Compute the rho median for each event and its standard deviation

     rhoMedian_ = bkgd_estimator.rho();
//    rhoSigma_ = bkgd_estimator.sigma();
//cout << rhoMedian_ << endl;

    // Determine the ptCut a la Soft Killer
    //-------------------------------------------------------------

    std::vector<double> ptConst_bkgd;


    for(fastjet::PseudoJet& jet : bkg_jets) {
        double maxpt_bkgd = 0;
        if(jet.is_pure_ghost()) continue;
        std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
        fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
      //  cout << ghosts.size() << endl;
     // Find the maximum-pt of patch in the background
       for(fastjet::PseudoJet p : particles) {
            double momentum = p.pt();
            if (momentum > maxpt_bkgd) maxpt_bkgd = momentum;
          }
       ptConst_bkgd.push_back(maxpt_bkgd);
       } // bkgd_jets loop
       std::nth_element(ptConst_bkgd.begin(), ptConst_bkgd.begin() + ptConst_bkgd.size()/2, ptConst_bkgd.end());
       double pTsoftKiller = ptConst_bkgd[ptConst_bkgd.size()/2];
      // cout << pTsoftKiller << endl;
      //  double alpha = 0.; // alpha is between (0,1). Controls how much we want to be like soft Killer. For alpha = 0 we recover the area-median result. For alpha = 1, soft Killer.
      // double pTModifiedsoftKiller = alpha * pTsoftKiller;

       double pTModifiedsoftKiller = 2.4;
       double deltaR_min = 0.1;
      // double deltaR_max = 0.4;
       // Determine the rho(cut): rho computed after pT < pTcut have been removed
       //-------------------------------------------------------------

       std::vector<double> rhoMedian_bkgd;
       std::vector<double> rhoCut_bkgd;

         for(fastjet::PseudoJet& jet : bkg_jets) {
             double rho_const = 0;
             double rho_cut = 0;

             std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
             fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
          // Find the maximum-pt of patch in the background
            for(fastjet::PseudoJet p : particles) {
                 double momentum = p.pt();
                 rho_const+=momentum;
                 double deltaR = jet.delta_R(p);
              //   if (deltaR > 0.4) cout << "YES: " << deltaR << endl;
                 if (momentum > pTModifiedsoftKiller || deltaR < deltaR_min){
                 rho_cut+=momentum;
                 }

               }
              rhoMedian_bkgd.push_back(rho_const/jet.area());
              rhoCut_bkgd.push_back(rho_cut/jet.area());
            //  cout <<"hh:" << jet.area() << endl;
            } // bkgd_jets loop

            std::nth_element(rhoMedian_bkgd.begin(), rhoMedian_bkgd.begin() + rhoMedian_bkgd.size()/2, rhoMedian_bkgd.end());

            std::nth_element(rhoCut_bkgd.begin(), rhoCut_bkgd.begin() + rhoCut_bkgd.size()/2, rhoCut_bkgd.end());

           rhoMedian_ = rhoMedian_bkgd[rhoMedian_bkgd.size()/2];
            //   cout << rhoCut_bkgd.size() << endl;
           rhoCut_ = rhoCut_bkgd[rhoCut_bkgd.size()/2];

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
      double truePatchPt = 0;
      double truePatchPx = 0;
      double truePatchPy = 0;
      int truePatchN = 0;
    //  double trueRho = 0;
      for(fastjet::PseudoJet p : particles) {
      // Obtain the true rho for each patch by using only the background particles
       truePatchPx+=p.px();
       truePatchPy+=p.py();
      // truePatchPt+=p.pt();
       truePatchN++;
       }
       truePatchPt = sqrt(pow(truePatchPx,2)+pow(truePatchPy,2));
    //if (ijet==0){cout << "Vectorial: "<< sqrt(pow(truePatchPx,2)+pow(truePatchPy,2)) << endl;
    //   cout << "Scalar: " << truePatchPt << endl;
    //   cout << "True:" << jet.pt() << endl;
     //}

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
      double pTsignal = 0;
      double pXcurrent = 0;
      double pYcurrent=0;
      double pXsignal = 0;
      double pYsignal = 0;

      std::fill(avail_part.begin(),avail_part.end(),1);

    //  while(std::accumulate(avail_part.begin(),avail_part.end(),0)>0 && particles_pTOrdered[nparticles+1].pt()<pTModifiedsoftKiller){
      while(std::accumulate(avail_part.begin(),avail_part.end(),0)>0 && particles_pTOrdered[nparticles+1].pt()<pTModifiedsoftKiller && jet.delta_R(particles_pTOrdered[nparticles+1])>deltaR_min){

      nparticles++;
      fastjet::PseudoJet partSel = particles_pTOrdered[nparticles];
      bkgEstimate.push_back(partSel);
    //  pTcurrent+= partSel.pt();
      pXcurrent+= partSel.px();
      pYcurrent+= partSel.py();
      avail_part[nparticles] = 0;
      if (partSel.user_info<PU14>().vertex_number() == 0){
        pXsignal+=partSel.px();
        pYsignal+=partSel.py();
       }
      }

      pTcurrent = sqrt(pow(pXcurrent,2)+pow(pYcurrent,2));
      pTsignal = sqrt(pow(pXsignal,2)+pow(pYsignal,2));
    //  if (ijet==0) cout << pTcurrent << endl;

  //    double pt_below = pTcurrent;

    //  while(std::accumulate(avail_part.begin(),avail_part.end(),0)>0 && pTcurrent<pt_below+rhoCut_*jet.area()){

  //    nparticles++;
  //    fastjet::PseudoJet partSel = particles_pTOrdered[nparticles];
  //    bkgEstimate.push_back(partSel);
  //    pTcurrent+= partSel.pt();
  //    avail_part[nparticles] = 0;
  //    }

    //Compute the shortest distance to the correlation line
    //----------------------------------------------------------

  //    double minus_slope_perp = 1/temperature;
  //    double intercept_perp = pTcurrent + (nparticles+1)*minus_slope_perp;
  //    int n_shortest = int(intercept_perp/(temperature+minus_slope_perp));
  //    double pTshortest = correlation_line(temperature, n_shortest);
  //    pTcurrent = pTshortest;
  //    }

      //double pT_estimate = pTcurrent;


    //  if (ijet==0) cout << pT_estimate << endl;

  //      cout << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;

      //Trial to correct for the signal contramination below the pt cut
      //------------------------------------------------------------------

    //    double pT_cut_corrected = nparticles*temperature*(1-(1+2*pTModifiedsoftKiller/temperature+2*pow(pTModifiedsoftKiller,2)/pow(temperature,2))*exp(-2*pTModifiedsoftKiller/temperature));
      //  if (pTcurrent > pT_cut_corrected) pTcurrent = pT_cut_corrected;
       double pT_estimate = pTcurrent+rhoCut_*jet.area();
       if (pT_estimate > truePatchPt) pT_estimate = truePatchPt;

      rho_estimate.push_back(pT_estimate);

  }

    return rho_estimate;


};

};

#endif

#ifndef simpleSubtractor_h
#define simpleSubtractor_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <map>
#include <random>
#include <numeric>
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TRandom3.h"


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "skEstimator.hh"
//#include "randomCones.hh"

#include "../PU14/PU14.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs a simplified version of the
// jet-by-jet shared layer subtraction
// Author: A. Soto-Ontoso, M. Verweij, Y. Mehtar-Tani
//---------------------------------------------------------------

class simpleSubtractor {

private :
  double jetRParam_;
  double ghostArea_;
  double ghostRapMax_;
  double jetRapMax_;
  double rho_;
  double rhoSigma_;

  std::vector<fastjet::PseudoJet> fjInputs_;
  std::vector<fastjet::PseudoJet> fjJetInputs_;

  std::vector<fastjet::PseudoJet> fjJetParticles_;
  //std::vector<std::vector<fastjet::PseudoJet>> fSigParticles;

  std::random_device rd_;

  fastjet::ClusterSequenceArea *csJetSub;

public :
   simpleSubtractor(double rJet = 0.4,
                        double ghostArea = 0.001,
                        double ghostRapMax = 3.0,
                        double jetRapMax = 3.0) :
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

  void clearMemory() {

  if(csJetSub) delete csJetSub;
}

std::vector<fastjet::PseudoJet> doSubtraction() {

fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, 0.001);
fastjet::GhostedAreaSpec ghost_spec_sub(ghostRapMax_, 1, 1.);
std::vector<fastjet::PseudoJet> jets = fjJetInputs_;
// do the clustering with ghosts and get the jets
//----------------------------------------------------------
fastjet::JetDefinition jet_def(antikt_algorithm, jetRParam_);
fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec);

fastjet::ClusterSequenceArea cs(fjInputs_, jet_def, area_def);
fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax_);
jets = fastjet::sorted_by_pt(jet_selector(cs.inclusive_jets()));

/// create what we need for the background estimation
//----------------------------------------------------------

fastjet::JetDefinition jet_estimate_bkgd(fastjet::kt_algorithm, 0.4);
fastjet::AreaDefinition area_estimate_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4)*(!fastjet::SelectorNHardest(2));

fastjet::ClusterSequenceArea csKt(fjInputs_, jet_estimate_bkgd, area_estimate_bkgd);
std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

   double total_pt_bkgd = 0;
   double maxpt_bkgd = 6;
   int nbins_bkg_pt = 5;
   double total_area_bkgd = 0;

    TH1D *h = (TH1D *)gROOT->FindObject("p_{T} const");
    delete h;
    TH1D *h_pt_constituents = new TH1D("p_{T} const", "p_{T} const", nbins_bkg_pt, 0.,maxpt_bkgd);

 std::vector<double> rhoMedian_bkgd;
 double rho_Median = 0;
for(fastjet::PseudoJet& jet : bkgd_jets) {
  double rho_const = 0;
  std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    for(fastjet::PseudoJet p : particles) {
        h_pt_constituents->Fill(p.pt(),p.pt());
         total_pt_bkgd+=p.pt();
         if(p.pt()<maxpt_bkgd)rho_const+=p.pt();
     }
     total_area_bkgd+= jet.area();
     rhoMedian_bkgd.push_back(rho_const/jet.area());
  } // bkgd jets loop

  std::nth_element(rhoMedian_bkgd.begin(), rhoMedian_bkgd.begin() + rhoMedian_bkgd.size()/2, rhoMedian_bkgd.end());

  rho_Median=rhoMedian_bkgd[rhoMedian_bkgd.size()/2];

  //cout << rho_Median << endl;

  h_pt_constituents->Scale(1./(total_pt_bkgd*total_area_bkgd)); //normalize


  double pt_binWidth = h_pt_constituents->GetBinWidth(0); // will use it later
  int nbins = h_pt_constituents->GetNbinsX();
  double ptmax = h_pt_constituents->GetBinCenter(nbins)+pt_binWidth/2;

  // Jet-by-jet subtraction
  //----------------------------------------------------------
  std::mt19937 rndSeed(rd_()); //rnd number generator seed
  std::vector<fastjet::PseudoJet> subtracted_jets;
  subtracted_jets.reserve(jets.size());
  int ijet = -1;
  int sampling = 0;
  int median = 0;
  int reduced = 0;
  for(fastjet::PseudoJet& jet : jets) {
    ++ijet;
    if(jet.is_pure_ghost()) continue;

    std::vector<fastjet::PseudoJet> particles, ghosts;
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    if(particles.size()<1 || jet.pt()<1.) continue;

    //set user_index of all particles to position particles vector

    for(int i = 0; i<(int)particles.size(); ++i) {
      particles[i].set_user_index(i);
    }

    // Create a new list of particles whose pT is below the bkg max Pt
  std::vector<fastjet::PseudoJet> particlesReduced;
 // Create a new list of particles whose pT is above the bkg max Pt and, therefore, considered to be signal
  std::vector<fastjet::PseudoJet> signalParticles;

    // Histogram to save the pT-information of each jet using the area

     TH1D *h_jet = (TH1D *)gROOT->FindObject("p_{T} const jet");
     delete h_jet;
     TH1D *h_pt_constituents_jet = new TH1D("p_{T} const jet", "p_{T} const jet", nbins_bkg_pt, 0.,maxpt_bkgd);

     double pT_reduced = 0; // pT stored by the particles below the bkg max Pt
     double patch_pT = 0; // total pT of the patch

     for(fastjet::PseudoJet p : particles) {
       double momentum = p.pt();
       patch_pT+=momentum;
       if (momentum<=maxpt_bkgd){
         h_pt_constituents_jet->Fill(momentum,momentum);
         particlesReduced.push_back(p);
         pT_reduced+=momentum;
       }
       else signalParticles.push_back(p);
    }


  for(int i = 0; i<(int)particlesReduced.size(); ++i) {
    particlesReduced[i].set_user_index(i);
  }

   double trueRho = 0;
    for(fastjet::PseudoJet p : particles) {
    if (abs(p.user_info<PU14>().vertex_number()) == 1) trueRho+=p.pt();
    }

    double pTbkg_estimate = rho_Median*jet.area();
  //  if(pTbkg_estimate > )  pTbkg_estimate = ;
  //  cout << "True: " << trueRho << "Area Median: "<< pTbkg_estimate <<  endl;


    //fjJetParticles_.clear();
    std::uniform_real_distribution<double> particleSelect(0., particlesReduced.size());

    // Depending on the values of the background estimate, the total pT and the pT below the bkg max Pt we have 4 options

    if(pTbkg_estimate == 0){
      for(fastjet::PseudoJet p : particles){
       fjJetParticles_.push_back(p);
     }
    }
    else if (pTbkg_estimate >= patch_pT)
        {
        //  fjJetParticles_.clear();
          median++;

        }

    else if (pTbkg_estimate >= pT_reduced)
       {
      for(fastjet::PseudoJet p : signalParticles){
        fjJetParticles_.push_back(p);
      };
      reduced++;
       }

     else {
      sampling++;

      h_pt_constituents->Scale(pTbkg_estimate*jet.area());
      double maxPtCurrent = 0.;
      std::vector<fastjet::PseudoJet> notBkgParticles = signalParticles;

  //Next step: figure out the probability of each particle to be background /only with momentum up to now
//----------------------------------------------------------
std::vector<double> prob_idx(particlesReduced.size(),0);
std::vector<fastjet::PseudoJet> BkgParticles(particlesReduced.size());
// One loop for the momentum
//------------------------------------------------------------------------
  for(int ic = 0; ic<(int)particlesReduced.size(); ++ic) {
    double candidate_pt = particlesReduced[ic].pt();
   int candidate_ptbin = int(candidate_pt*nbins/ptmax)+1;
   double candidate_pt_mean = h_pt_constituents->GetBinContent(candidate_ptbin);
   if(candidate_pt_mean==0) candidate_pt_mean = (h_pt_constituents->GetBinContent(candidate_ptbin+1)+h_pt_constituents->GetBinContent(candidate_ptbin-1))/2;

   double candidate_pt_error = h_pt_constituents->GetBinError(candidate_ptbin);
   double candidate_pt_prob = candidate_pt_mean + candidate_pt_error;

   double upper_limit_mean = h_pt_constituents_jet->GetBinContent(candidate_ptbin);
   double upper_limit_error = h_pt_constituents_jet->GetBinError(candidate_ptbin);
   double upper_limit = upper_limit_mean + upper_limit_error;

   if (candidate_pt_prob >= upper_limit){ //When the signal+background is below the background
    prob_idx[particlesReduced[ic].user_index()]=1;
   }
   else{
     std::uniform_real_distribution<double> distPt(0., upper_limit);

     double random_prob = distPt(rndSeed);
     double prob_bkg = candidate_pt_prob/random_prob;
     if (prob_bkg>1) prob_bkg =1;
     prob_idx[particlesReduced[ic].user_index()]=prob_bkg;
    }
    prob_idx[particlesReduced[ic].user_index()] = abs(prob_idx[particlesReduced[ic].user_index()]);
  }

  //sort according to their probability values
  //----------------------------------------------------------
  // initialize original index locations
  std::vector<size_t> ish(prob_idx.size());
  iota(ish.begin(), ish.end(), 0);
  std::sort(ish.begin(), ish.end(),
            [&prob_idx](size_t i1, size_t i2) {return prob_idx[i1] > prob_idx[i2];});
// For Fisher method one has to revert the order
  //create final UE object
  //----------------------------------------------------------

        for(auto userIndex : ish) {
          fastjet::PseudoJet part = particlesReduced[userIndex];
          if(maxPtCurrent<=pTbkg_estimate) { //assign as bkgd particle
              maxPtCurrent+=part.pt();
              BkgParticles.push_back(part);
              }
             else {
               notBkgParticles.push_back(part);
            //   if(ijet == 0) cout << "eo" << endl;
              }
            }
          //   if(ijet==0) cout << maxPtCurrent << pTbkg_estimate << endl;
             for (fastjet::PseudoJet p : notBkgParticles){
             fjJetParticles_.push_back(p);
             }

             if(maxPtCurrent>pTbkg_estimate && BkgParticles.size() > 1) {
              double maxPtPrev = 0;
              std::vector<fastjet::PseudoJet> initCondParticles;
               for(int ic = 0; ic<BkgParticles.size()-1; ++ic) {
                 initCondParticles.push_back(BkgParticles[ic]);
                 maxPtPrev+=initCondParticles.at(ic).pt();
               }
          //     if(ijet ==0) cout << "Current: " << maxPtCurrent << "Prev: " << maxPtPrev << endl;
               double distance_one = sqrt((pTbkg_estimate-maxPtCurrent)*(pTbkg_estimate-maxPtCurrent));
               double distance_two = sqrt((pTbkg_estimate-maxPtPrev)*(pTbkg_estimate-maxPtPrev));

               if (distance_one > distance_two) {
             //    if(ijet==0) cout<<"yes"<< maxPtPrev<<endl;
                 fjJetParticles_.push_back(BkgParticles[BkgParticles.size()-1]);
              //   BkgParticles.pop_back();


               }
             }
} // end of the else sampling loop


 h_pt_constituents->Scale(1/(pTbkg_estimate*jet.area()));
//   fSigParticles.push_back(fjJetParticles_);
  } // end of jets-loop
  // cout << fjJetParticles_.size() << endl;
   //fastjet::JetDefinition jet_defSub(antikt_algorithm, 0.4);
   //fastjet::AreaDefinition area_def_sub = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec_sub);
   //fastjet::ClusterSequenceArea csSub(fjJetParticles_, jet_defSub, area_def_sub);
   //subtracted_jets = fastjet::sorted_by_pt(selector(csSub.inclusive_jets()));;

//  std::cout << "\n n signal particles " << fjJetParticles_.size() << std::endl;
  return fjJetParticles_;

} // end of doSubtraction


};

#endif




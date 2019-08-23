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
  std::vector<fastjet::PseudoJet> constituents1_;
  std::vector<fastjet::PseudoJet> constituents2_;
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

fastjet::JetDefinition jet_defSub(antikt_algorithm, 20.);
fastjet::AreaDefinition area_def_sub = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec_sub);
fastjet::ClusterSequenceArea csSub(fjInputs_, jet_defSub, area_def_sub);

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


for(fastjet::PseudoJet& jet : bkgd_jets) {

  std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    for(fastjet::PseudoJet p : particles) {
        h_pt_constituents->Fill(p.pt(),p.pt());
         total_pt_bkgd+=p.pt();
     }
     total_area_bkgd+= jet.area();
  }


  h_pt_constituents->Scale(1./((double)total_pt_bkgd)); //normalize


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
  constituents1_.reserve(fjInputs_.size());
  constituents2_.reserve(fjInputs_.size());
  for(fastjet::PseudoJet& jet : jets) {
    ++ijet;
    if(jet.is_pure_ghost()) continue;
    fjJetParticles_.clear();
    std::vector<fastjet::PseudoJet> particles, ghosts;
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    if(particles.size()<1 || jet.pt()<1.) continue;

    fastjet::JetDefinition jet_def_ca(fastjet::cambridge_algorithm,999.);
    fastjet::ClusterSequence cs(particles, jet_def);
    std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());
    fastjet::PseudoJet transformedJet = tempJets[0];
    fastjet::PseudoJet j1, j2;

    if(transformedJet.has_parents(j1, j2) == true)
    {
    //  std::vector<fastjet::PseudoJet> particles, ghosts;
    //  fastjet::SelectorIsPureGhost().sift(j1.constituents(), ghosts, particles);
       constituents1_=j1.constituents();
       constituents2_=j2.constituents();
    }

    // Now I have to repeat the same steps for the two subjets

    // First subjet

    //set user_index of all particles to position particles vector
  //  cout << constituents1_[0].[0].pt() << endl;

    for(int i = 0; i<(int)constituents1_.size(); ++i) {
      constituents1_[i].set_user_index(i);
    }

    // Create a new list of particles whose pT is below the bkg max Pt
  std::vector<fastjet::PseudoJet> particlesReduced_one;
 // Create a new list of particles whose pT is above the bkg max Pt and, therefore, considered to be signal
  std::vector<fastjet::PseudoJet> signalParticles_one;

    // Histogram to save the pT-information of each jet using the area

     TH1D *h_jet_one = (TH1D *)gROOT->FindObject("p_{T} const jet one");
     delete h_jet_one;
     TH1D *h_pt_constituents_jet_one = new TH1D("p_{T} const jet one", "p_{T} const jet one", nbins_bkg_pt, 0.,maxpt_bkgd);

     double pT_reduced_one = 0; // pT stored by the particles below the bkg max Pt
     double patch_pT_one = 0; // total pT of the patch

     for(fastjet::PseudoJet p : constituents1_) {
       double momentum = p.pt();
       patch_pT_one+=momentum;
       if (momentum<=maxpt_bkgd){
         h_pt_constituents_jet_one->Fill(momentum,momentum);
         particlesReduced_one.push_back(p);
         pT_reduced_one+=momentum;
       }
       else signalParticles_one.push_back(p);
    }


  for(int i = 0; i<(int)particlesReduced_one.size(); ++i) {
    particlesReduced_one[i].set_user_index(i);
  }

   double trueRho_one = 0;
    for(fastjet::PseudoJet p : constituents1_) {
    if (abs(p.user_info<PU14>().vertex_number()) == 1) trueRho_one+=p.pt();
   }

    double pTbkg_estimate_one = trueRho_one;


    //fjJetParticles_.clear();
    std::uniform_real_distribution<double> particleSelect(0., particlesReduced_one.size());

    // Depending on the values of the background estimate, the total pT and the pT below the bkg max Pt we have 4 options

    if(pTbkg_estimate_one == 0){
      for(fastjet::PseudoJet p : constituents1_){
       fjJetParticles_.push_back(p);
     }
    }
    else if (pTbkg_estimate_one >= patch_pT_one)
        {
        //  fjJetParticles_.clear();
          median++;

        }

    else if (pTbkg_estimate_one >= pT_reduced_one)
       {
      for(fastjet::PseudoJet p : signalParticles_one){
        fjJetParticles_.push_back(p);
      };
      reduced++;
       }

     else {
      sampling++;

      h_pt_constituents->Scale(pTbkg_estimate_one);
      double maxPtCurrent_one = 0.;
      std::vector<fastjet::PseudoJet> notBkgParticles_one = signalParticles_one;

  //Next step: figure out the probability of each particle to be background /only with momentum up to now
//----------------------------------------------------------
std::vector<double> prob_idx_one(particlesReduced_one.size(),0);
std::vector<fastjet::PseudoJet> BkgParticles_one(particlesReduced_one.size());
// One loop for the momentum
//------------------------------------------------------------------------
  for(int ic = 0; ic<(int)particlesReduced_one.size(); ++ic) {
    double candidate_pt = particlesReduced_one[ic].pt();
   int candidate_ptbin = int(candidate_pt*nbins/ptmax)+1;
   double candidate_pt_mean = h_pt_constituents->GetBinContent(candidate_ptbin);
   if(candidate_pt_mean==0) candidate_pt_mean = (h_pt_constituents->GetBinContent(candidate_ptbin+1)+h_pt_constituents->GetBinContent(candidate_ptbin-1))/2;

   double candidate_pt_error = h_pt_constituents->GetBinError(candidate_ptbin);
   double candidate_pt_prob = candidate_pt_mean + candidate_pt_error;

   double upper_limit_mean = h_pt_constituents_jet_one->GetBinContent(candidate_ptbin);
   double upper_limit_error = h_pt_constituents_jet_one->GetBinError(candidate_ptbin);
   double upper_limit = upper_limit_mean + upper_limit_error;

   if (candidate_pt_prob >= upper_limit){ //When the signal+background is below the background
    prob_idx_one[particlesReduced_one[ic].user_index()]=1;
   }
   else{
     std::uniform_real_distribution<double> distPt(0., upper_limit);

     double random_prob = distPt(rndSeed);
     double prob_bkg = candidate_pt_prob/random_prob;
     if (prob_bkg>1) prob_bkg =1;
     prob_idx_one[particlesReduced_one[ic].user_index()]=prob_bkg;
    }
    prob_idx_one[particlesReduced_one[ic].user_index()] = abs(prob_idx_one[particlesReduced_one[ic].user_index()]);
  }

  //sort according to their probability values
  //----------------------------------------------------------
  // initialize original index locations
  std::vector<size_t> ish_one(prob_idx_one.size());
  iota(ish_one.begin(), ish_one.end(), 0);
  std::sort(ish_one.begin(), ish_one.end(),
            [&prob_idx_one](size_t i1, size_t i2) {return prob_idx_one[i1] > prob_idx_one[i2];});
// For Fisher method one has to revert the order
  //create final UE object
  //----------------------------------------------------------

        for(auto userIndex : ish_one) {
          fastjet::PseudoJet part = particlesReduced_one[userIndex];
          if(maxPtCurrent_one<=pTbkg_estimate_one) { //assign as bkgd particle
              maxPtCurrent_one+=part.pt();
              BkgParticles_one.push_back(part);
              }
             else {
               notBkgParticles_one.push_back(part);
            //   if(ijet == 0) cout << "eo" << endl;
              }
            }
          //   if(ijet==0) cout << maxPtCurrent << pTbkg_estimate << endl;
             for (fastjet::PseudoJet p : notBkgParticles_one){
             fjJetParticles_.push_back(p);
             }

             if(maxPtCurrent_one>pTbkg_estimate_one && BkgParticles_one.size() > 1) {
              double maxPtPrev = 0;
              std::vector<fastjet::PseudoJet> initCondParticles;
               for(int ic = 0; ic<BkgParticles_one.size()-1; ++ic) {
                 initCondParticles.push_back(BkgParticles_one[ic]);
                 maxPtPrev+=initCondParticles.at(ic).pt();
               }
          //     if(ijet ==0) cout << "Current: " << maxPtCurrent << "Prev: " << maxPtPrev << endl;
               double distance_one = sqrt((pTbkg_estimate_one-maxPtCurrent_one)*(pTbkg_estimate_one-maxPtCurrent_one));
               double distance_two = sqrt((pTbkg_estimate_one-maxPtPrev)*(pTbkg_estimate_one-maxPtPrev));

               if (distance_one > distance_two) {
             //    if(ijet==0) cout<<"yes"<< maxPtPrev<<endl;
                 fjJetParticles_.push_back(BkgParticles_one[BkgParticles_one.size()-1]);
              //   BkgParticles.pop_back();


               }
             }
} // end of the else sampling loop


 h_pt_constituents->Scale(1/pTbkg_estimate_one);

// Second subjet

for(int i = 0; i<(int)constituents2_.size(); ++i) {
  constituents2_[i].set_user_index(i);
}

// Create a new list of particles whose pT is below the bkg max Pt
std::vector<fastjet::PseudoJet> particlesReduced_two;
// Create a new list of particles whose pT is above the bkg max Pt and, therefore, considered to be signal
std::vector<fastjet::PseudoJet> signalParticles_two;

// Histogram to save the pT-information of each jet using the area

 TH1D *h_jet_two = (TH1D *)gROOT->FindObject("p_{T} const jet two");
 delete h_jet_two;
 TH1D *h_pt_constituents_jet_two = new TH1D("p_{T} const jet two", "p_{T} const jet two", nbins_bkg_pt, 0.,maxpt_bkgd);

 double pT_reduced_two = 0; // pT stored by the particles below the bkg max Pt
 double patch_pT_two = 0; // total pT of the patch

 for(fastjet::PseudoJet p : constituents2_) {
   double momentum = p.pt();
   patch_pT_two+=momentum;
   if (momentum<=maxpt_bkgd){
     h_pt_constituents_jet_two->Fill(momentum,momentum);
     particlesReduced_two.push_back(p);
     pT_reduced_two+=momentum;
   }
   else signalParticles_two.push_back(p);
}


for(int i = 0; i<(int)particlesReduced_two.size(); ++i) {
particlesReduced_two[i].set_user_index(i);
}

double trueRho_two = 0;
for(fastjet::PseudoJet p : constituents2_) {
if (abs(p.user_info<PU14>().vertex_number()) == 1) trueRho_two+=p.pt();
}

double pTbkg_estimate_two = trueRho_two;


//fjJetParticles_.clear();
std::uniform_real_distribution<double> particleSelect_two(0., particlesReduced_two.size());

// Depending on the values of the background estimate, the total pT and the pT below the bkg max Pt we have 4 options

if(pTbkg_estimate_two == 0){
  for(fastjet::PseudoJet p : constituents2_){
   fjJetParticles_.push_back(p);
 }
}
else if (pTbkg_estimate_two >= patch_pT_two)
    {
    //  fjJetParticles_.clear();
      median++;

    }

else if (pTbkg_estimate_two >= pT_reduced_two)
   {
  for(fastjet::PseudoJet p : signalParticles_two){
    fjJetParticles_.push_back(p);
  };
  reduced++;
   }

 else {
  sampling++;

  h_pt_constituents->Scale(pTbkg_estimate_two);
  double maxPtCurrent_two = 0.;
  std::vector<fastjet::PseudoJet> notBkgParticles_two = signalParticles_two;

//Next step: figure out the probability of each particle to be background /only with momentum up to now
//----------------------------------------------------------
std::vector<double> prob_idx_two(particlesReduced_two.size(),0);
std::vector<fastjet::PseudoJet> BkgParticles_two(particlesReduced_two.size());
// One loop for the momentum
//------------------------------------------------------------------------
for(int ic = 0; ic<(int)particlesReduced_two.size(); ++ic) {
double candidate_pt = particlesReduced_two[ic].pt();
int candidate_ptbin = int(candidate_pt*nbins/ptmax)+1;
double candidate_pt_mean = h_pt_constituents->GetBinContent(candidate_ptbin);
if(candidate_pt_mean==0) candidate_pt_mean = (h_pt_constituents->GetBinContent(candidate_ptbin+1)+h_pt_constituents->GetBinContent(candidate_ptbin-1))/2;

double candidate_pt_error = h_pt_constituents->GetBinError(candidate_ptbin);
double candidate_pt_prob = candidate_pt_mean + candidate_pt_error;

double upper_limit_mean = h_pt_constituents_jet_two->GetBinContent(candidate_ptbin);
double upper_limit_error = h_pt_constituents_jet_two->GetBinError(candidate_ptbin);
double upper_limit = upper_limit_mean + upper_limit_error;

if (candidate_pt_prob >= upper_limit){ //When the signal+background is below the background
prob_idx_two[particlesReduced_two[ic].user_index()]=1;
}
else{
 std::uniform_real_distribution<double> distPt(0., upper_limit);

 double random_prob = distPt(rndSeed);
 double prob_bkg = candidate_pt_prob/random_prob;
 if (prob_bkg>1) prob_bkg =1;
 prob_idx_two[particlesReduced_two[ic].user_index()]=prob_bkg;
}
prob_idx_two[particlesReduced_two[ic].user_index()] = abs(prob_idx_two[particlesReduced_two[ic].user_index()]);
}

//sort according to their probability values
//----------------------------------------------------------
// initialize original index locations
std::vector<size_t> ish_two(prob_idx_two.size());
iota(ish_two.begin(), ish_two.end(), 0);
std::sort(ish_two.begin(), ish_two.end(),
        [&prob_idx_two](size_t i1, size_t i2) {return prob_idx_two[i1] > prob_idx_two[i2];});
// For Fisher method one has to revert the order
//create final UE object
//----------------------------------------------------------

    for(auto userIndex : ish_two) {
      fastjet::PseudoJet part = particlesReduced_two[userIndex];
      if(maxPtCurrent_two<=pTbkg_estimate_two) { //assign as bkgd particle
          maxPtCurrent_two+=part.pt();
          BkgParticles_two.push_back(part);
          }
         else {
           notBkgParticles_two.push_back(part);
        //   if(ijet == 0) cout << "eo" << endl;
          }
        }
      //   if(ijet==0) cout << maxPtCurrent << pTbkg_estimate << endl;
         for (fastjet::PseudoJet p : notBkgParticles_two){
         fjJetParticles_.push_back(p);
         }

         if(maxPtCurrent_two>pTbkg_estimate_two && BkgParticles_two.size() > 1) {
          double maxPtPrev = 0;
          std::vector<fastjet::PseudoJet> initCondParticles;
           for(int ic = 0; ic<BkgParticles_two.size()-1; ++ic) {
             initCondParticles.push_back(BkgParticles_two[ic]);
             maxPtPrev+=initCondParticles.at(ic).pt();
           }
      //     if(ijet ==0) cout << "Current: " << maxPtCurrent << "Prev: " << maxPtPrev << endl;
           double distance_one = sqrt((pTbkg_estimate_two-maxPtCurrent_two)*(pTbkg_estimate_two-maxPtCurrent_two));
           double distance_two = sqrt((pTbkg_estimate_two-maxPtPrev)*(pTbkg_estimate_two-maxPtPrev));

           if (distance_one > distance_two) {
         //    if(ijet==0) cout<<"yes"<< maxPtPrev<<endl;
             fjJetParticles_.push_back(BkgParticles_two[BkgParticles_two.size()-1]);
          //   BkgParticles.pop_back();


           }
         }
} // end of the else sampling loop


h_pt_constituents->Scale(1/pTbkg_estimate_two);

if(fjJetParticles_.size()>0) {
   csJetSub = new fastjet::ClusterSequenceArea(fjJetParticles_, jet_defSub,    area_def_sub);
     std::vector<fastjet::PseudoJet> jetSub = fastjet::sorted_by_pt(csJetSub->inclusive_jets());
     if(jetSub[0].pt()>0.) subtracted_jets.push_back(jetSub[0]);
     if(subtracted_jets.size()>0 && subtracted_jets.size()<2)
     csJetSub->delete_self_when_unused();
   }
//   fSigParticles.push_back(fjJetParticles_);
  } // end of jets-loop
//}
//}

  return subtracted_jets;

} // end of doSubtraction


};

#endif

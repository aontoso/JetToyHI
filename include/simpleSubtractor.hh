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

// create what we need for the background estimation with Random cones
//----------------------------------------------------------------------


  //TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 1200, 800);
  //h_bkg_constituents->Draw("LEGO");
  //c1->SaveAs("Prueba_3.C");

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


    double total_pt_bkgd = 0;
    double maxdeltar_bkgd = 0.5;
    double maxpt_bkgd = 6;
    int nbins_bkg_deltar = 5;
    int nbins_bkg_pt = 10;
    TH2D *h_bkg = (TH2D *)gROOT->FindObject("bkg const");
    delete h_bkg;
    TH2D *h_bkg_constituents = new TH2D("bkg const", "bkg const", nbins_bkg_deltar, 0.,maxdeltar_bkgd, nbins_bkg_pt, 0.,maxpt_bkgd);

    double cone_rapidity = -3+jetRParam_;
    double etaMax_ = cone_rapidity;
    double minPhi = 0.;
    double maxPhi = 2.*TMath::Pi();
    int nCones_= 5;
    TRandom3*rnd_ = new TRandom3(0);

    for(unsigned int i = 0; i<nCones_; ++i) {

      //pick random position for random cone
      double etaRC = rnd_->Rndm() * (etaMax_ - -1.*etaMax_) + -1.*etaMax_;
      double phiRC = rnd_->Rndm() * (maxPhi - minPhi) + minPhi;
      for(fastjet::PseudoJet part : fjInputs_) {
        if(part.pt()<maxpt_bkgd){
      //  double dr = deltaR(part.phi(),phiRC,part.eta(),etaRC);
        double dPhi = part.phi() - phiRC;
        double dEta =part.eta() - etaRC;
        dPhi = TVector2::Phi_mpi_pi(dPhi);
        double dr = TMath::Sqrt(dPhi * dPhi + dEta * dEta);
        if(dr<jetRParam_) {
          h_bkg_constituents->Fill(dr,part.pt(),part.pt());
          total_pt_bkgd+=part.pt();
      }
        }
      }
    }

      h_bkg_constituents->Scale(1./((double)total_pt_bkgd)); //normalize


    // Create a new list of particles whose pT is below the bkg max Pt
     std::vector<fastjet::PseudoJet> particlesReduced;
    // Create a new list of particles whose pT is above the bkg max Pt and, therefore, considered to be signal
     std::vector<fastjet::PseudoJet> signalParticles;

    // Histogram to save the pT-information of each jet using the area

     TH2D *h_jet = (TH2D *)gROOT->FindObject("const jet");
     delete h_jet;
     TH2D *h_constituents_jet = new TH2D("const jet", "const jet", nbins_bkg_deltar, 0.,maxdeltar_bkgd,nbins_bkg_pt,0.,maxpt_bkgd);

     double pT_reduced = 0; // pT stored by the particles below the bkg max Pt
     double patch_pT = 0; // total pT of the patch


     for(fastjet::PseudoJet p : particles) {
       double momentum = p.pt();
       double deltar = jet.delta_R(p);
       patch_pT+=momentum;
       if (momentum<=maxpt_bkgd && deltar<maxdeltar_bkgd){
         particlesReduced.push_back(p);
         pT_reduced+=momentum;
         h_constituents_jet->Fill(deltar, momentum, momentum);
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

        double pTbkg_estimate = trueRho;

    fjJetParticles_.clear();
    std::uniform_real_distribution<double> particleSelect(0., particlesReduced.size());

    // Depending on the values of the background estimate, the total pT and the pT below the bkg max Pt we have 4 options

    if(pTbkg_estimate == 0){
      fjJetParticles_ = particles;
    }
    else if (pTbkg_estimate >= patch_pT)
        {
          fjJetParticles_.clear();
          median++;

        }

    else if (pTbkg_estimate >= pT_reduced)
       {
      fjJetParticles_ = signalParticles;
      reduced++;
       }

     else {
      sampling++;

      h_bkg_constituents->Scale(pTbkg_estimate);
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
   double candidate_deltar = jet.delta_R(particlesReduced[ic]);
   int candidate_bin = h_bkg_constituents->FindFixBin(candidate_deltar,candidate_pt);
   double candidate_mean = h_bkg_constituents->GetBinContent(candidate_bin);
   if(candidate_mean==0) candidate_mean = (h_bkg_constituents->GetBinContent(candidate_bin+1)+h_bkg_constituents->GetBinContent(candidate_bin-1))/2;

   double candidate_error = h_bkg_constituents->GetBinError(candidate_bin);
   double candidate_prob = candidate_mean + candidate_error;

   double upper_limit_mean = h_constituents_jet->GetBinContent(candidate_bin);
   double upper_limit_error = h_constituents_jet->GetBinError(candidate_bin);
   double upper_limit = upper_limit_mean + upper_limit_error;

   if (candidate_prob >= upper_limit){ //When the signal+background is below the background
    prob_idx[particlesReduced[ic].user_index()]=1;
   }
   else{
     std::uniform_real_distribution<double> distPt(0., upper_limit);

     double random_prob = distPt(rndSeed);
     double prob_bkg = candidate_prob/random_prob;
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
             }
            }
            fjJetParticles_=notBkgParticles;

} // end of the else sampling loop
 h_bkg_constituents->Scale(1/pTbkg_estimate);

    if(fjJetParticles_.size()>0) {
    csJetSub = new fastjet::ClusterSequenceArea(fjJetParticles_, jet_defSub,    area_def_sub);
      std::vector<fastjet::PseudoJet> jetSub = fastjet::sorted_by_pt(csJetSub->inclusive_jets());
      if(jetSub[0].pt()>0.) subtracted_jets.push_back(jetSub[0]);
      if(subtracted_jets.size()>0 && subtracted_jets.size()<2)
      csJetSub->delete_self_when_unused();
    //  if(ijet==0) cout << jetSub[0].pt() << endl;
    }
  } // end of jets-loop

//  std::cout << "\n n subtracted jets: " << subtracted_jets.size() << "Number of jets: " << counter << "Sampling applied: " << sampling << "Median larger than pT" << median << "Reduced part too small: " << reduced << std::endl;
  return subtracted_jets;

} // end of doSubtraction


};

#endif

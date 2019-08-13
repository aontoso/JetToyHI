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

// create what we need for the background estimation
//----------------------------------------------------------

//fastjet::JetDefinition jet_estimate_bkgd(fastjet::kt_algorithm, 0.4);
//fastjet::AreaDefinition area_estimate_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
//fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4)*(!fastjet::SelectorNHardest(2));

//fastjet::ClusterSequenceArea csKt(fjInputs_, jet_estimate_bkgd, area_estimate_bkgd);
//std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

double maxpt_bkgd = 0;
//int nparticles_bkg = 0;
double total_pt_bkgd = 0;

//for(fastjet::PseudoJet& jet : bkgd_jets) {
//  std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
//    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

 // Find the maximum-pt of the background
//   for(fastjet::PseudoJet p : particles) {
//        double momentum = p.pt();
//        if (momentum > maxpt_bkgd) maxpt_bkgd = momentum;
//        nparticles_bkg++;
//     }
//   } // bkgd_jets loop

  //   maxpt_bkgd = 6;

  //  int nbins_bkg = 10;
  //  TH1D *h = (TH1D *)gROOT->FindObject("p_{T} const");
  //  delete h;
  //  TH1D *h_pt_constituents = new TH1D("p_{T} const", "p_{T} const", nbins_bkg, 0.,maxpt_bkgd);


//for(fastjet::PseudoJet& jet : bkgd_jets) {

//  std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
//    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
//    for(fastjet::PseudoJet p : particles) {
//        h_pt_constituents->Fill(p.pt(),p.pt());
//         total_pt_bkgd+=p.pt();
//     }
//  }


  //h_pt_constituents->Scale(1./((double)total_pt_bkgd)); //normalize


//  double pt_binWidth = h_pt_constituents->GetBinWidth(0); // will use it later
//  int nbins = h_pt_constituents->GetNbinsX();
//  double ptmax = h_pt_constituents->GetBinCenter(nbins)+pt_binWidth/2;

double maxdeltar_bkgd = 0.5;
int nbins_bkg_deltar = 5;
TH1D *h_cone = (TH1D *)gROOT->FindObject("deltaR const");
delete h_cone;
TH1D *h_deltar_constituents = new TH1D("deltaR const", "deltaR const", nbins_bkg_deltar, 0.,maxdeltar_bkgd);

double binWidth_bkgd_deltaR = h_deltar_constituents->GetBinWidth(0);
double cone_rapidity = -3+jetRParam_;
double etaMax_ = cone_rapidity;
double minPhi = 0.;
double maxPhi = 2.*TMath::Pi();
int nCones_= 10;
TRandom3*rnd_ = new TRandom3(0);

maxpt_bkgd = 6;

int nbins_bkg = 10;
TH1D *h = (TH1D *)gROOT->FindObject("p_{T} const");
delete h;
TH1D *h_pt_constituents = new TH1D("p_{T} const", "p_{T} const", nbins_bkg, 0.,maxpt_bkgd);

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
      //std::cout << "phi: " << part.phi() << " " << phiRC << " eta: " << part.eta() << " " << etaRC << std::endl;
      h_deltar_constituents->Fill(dr);
      h_pt_constituents->Fill(part.pt(),part.pt());
      total_pt_bkgd+=part.pt();
  }
    }
  }
}

h_deltar_constituents->Scale(1./((double)h_deltar_constituents->GetEntries())); //normalize
double deltarmax = h_deltar_constituents->GetBinCenter(nbins_bkg_deltar)+binWidth_bkgd_deltaR/2;

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


  //  randomCones rigidcone(5, jetRParam_, cone_rapidity,maxpt_bkgd);
  //  rigidcone.setInputParticles(fjInputs_);
  //  std::vector<std::vector<double>> cones;
  //  cones = rigidcone.run();

// Random cones constituents pT-spectrum
//----------------------------------------------------------

//double maxdeltar_bkgd = 0;
//int nparticles_bkg = 0;

//for(int i = 0; i<cones.size(); i++) {
// Find the maximum-deltaR of the background
//for(int j = 0; j<cones[i].size(); j++) {
//     double deltar = cones[i][j];
//     if (deltar > maxdeltar_bkgd) maxdeltar_bkgd = deltar;
//     nparticles_bkg++;
//  }
 //} // random_cone loop
//
 //int nbins_bkg_deltar = 5;
 //TH1D *h_cone = (TH1D *)gROOT->FindObject("deltaR const");
 //delete h_cone;
 //TH1D *h_deltar_constituents = new TH1D("deltaR const", "deltaR const", //nbins_bkg_deltar, 0.,maxdeltar_bkgd);

//double binWidth_bkgd_deltaR = h_deltar_constituents->GetBinWidth(0);
//for(int i = 0; i<cones.size(); i++) {
//for(int j = 0; j<cones[i].size(); j++) {
//    double deltar = cones[i][j];
//    h_deltar_constituents->Fill(deltar);
 //}
//}

//      cout << maxpt_bkgd << endl;

    // Create a new list of particles whose pT is below the bkg max Pt
     std::vector<fastjet::PseudoJet> particlesReduced;
    // Create a new list of particles whose pT is above the bkg max Pt and, therefore, considered to be signal
     std::vector<fastjet::PseudoJet> signalParticles;

    // Histogram to save the pT-information of each jet using the area

     TH1D *h_jet = (TH1D *)gROOT->FindObject("p_{T} const jet");
     delete h_jet;
     TH1D *h_pt_constituents_jet = new TH1D("p_{T} const jet", "p_{T} const jet", nbins, 0.,maxpt_bkgd);

     double pT_reduced = 0; // pT stored by the particles below the bkg max Pt
     double patch_pT = 0; // total pT of the patch

     TH1D *h_deltar = (TH1D *)gROOT->FindObject("deltaR const jet");
     delete h_deltar;
     TH1D *h_deltar_constituents_jet = new TH1D("deltaR const jet", "deltaR const jet", nbins_bkg_deltar, 0.,maxdeltar_bkgd);


     for(fastjet::PseudoJet p : particles) {
       double momentum = p.pt();
       double deltar = jet.delta_R(p);
       patch_pT+=momentum;
       if (momentum<=maxpt_bkgd){
         h_pt_constituents_jet->Fill(momentum,momentum);
         particlesReduced.push_back(p);
         pT_reduced+=momentum;
         if (deltar < maxdeltar_bkgd) h_deltar_constituents_jet->Fill(deltar);
       }
       else signalParticles.push_back(p);
    }

        h_deltar_constituents_jet->Scale(1./((double)h_deltar_constituents_jet->GetEntries()));

  for(int i = 0; i<(int)particlesReduced.size(); ++i) {
    particlesReduced[i].set_user_index(i);
  }

   double trueRho = 0;
    for(fastjet::PseudoJet p : particles) {
    if (abs(p.user_info<PU14>().vertex_number()) == 1) trueRho+=p.pt();
   }

//  double pTbkg_estimate = rho_*jet.area();
        double pTbkg_estimate = trueRho;
//     double rho_ = trueRho/jet.area();

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
    //  h_pt_constituents_jet->Scale(pT_reduced);
    // cout << "Before: " << h_pt_constituents->GetBinContent(1) << endl;
      h_pt_constituents->Scale(pTbkg_estimate);
      double maxPtCurrent = 0.;
      std::vector<fastjet::PseudoJet> notBkgParticles = signalParticles;
      //     cout << "After: " << h_pt_constituents->GetBinContent(1) << endl;
      // cout << pTbkg_estimate << endl;
//      std::random_shuffle(particlesReduced.begin(), particlesReduced.end()); // randomize the vector

  //Next step: figure out the probability of each particle to be background /only with momentum up to now
//----------------------------------------------------------
std::vector<double> prob_idx(particlesReduced.size(),0);
std::vector<fastjet::PseudoJet> BkgParticles(particlesReduced.size());
//Double_t pbkg_momentum[particlesReduced.size()];
//Double_t pbkg_theta[particlesReduced.size()];
// One loop for the momentum
//------------------------------------------------------------------------
  for(int ic = 0; ic<(int)particlesReduced.size(); ++ic) {
//    if (ijet==0)cout << particles[indices[ic]].user_index() << endl;
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
     //cout << prob_bkg << endl;
     prob_idx[particlesReduced[ic].user_index()]=prob_bkg;
    }
    prob_idx[particlesReduced[ic].user_index()] = abs(prob_idx[particlesReduced[ic].user_index()]);
  //  pbkg_momentum[ic]=prob_idx[particlesReduced[ic].user_index()];
    //if(ijet==0) cout << pbkg_momentum[ic] << endl;
  }

  // One loop for the angle
  //------------------------------------------------------------------------

  for(int ic = 0; ic<(int)particlesReduced.size(); ++ic) {
//    if (ijet==0)cout << particles[indices[ic]].user_index() << endl;
   double candidate_deltar = jet.delta_R(particlesReduced[ic]);
   int candidate_deltarbin = int(candidate_deltar*nbins_bkg_deltar/deltarmax)+1;
   double candidate_deltar_mean = h_deltar_constituents->GetBinContent(candidate_deltarbin);
   if(candidate_deltar_mean==0) candidate_deltar_mean = (h_deltar_constituents->GetBinContent(candidate_deltarbin+1)+h_deltar_constituents->GetBinContent(candidate_deltarbin-1))/2;

   double candidate_deltar_error = h_deltar_constituents->GetBinError(candidate_deltarbin);
   double candidate_deltar_prob = candidate_deltar_mean + candidate_deltar_error;

   double upper_limit_mean = h_deltar_constituents_jet->GetBinContent(candidate_deltarbin);
   double upper_limit_error = h_deltar_constituents_jet->GetBinError(candidate_deltarbin);
   double upper_limit = upper_limit_mean + upper_limit_error;

   if (candidate_deltar_prob >= upper_limit){ //When the signal+background is below the background
    prob_idx[particlesReduced[ic].user_index()]=1;
  //  pbkg_theta[ic]=1;
   }
   else{
     std::uniform_real_distribution<double> distPt(0., upper_limit);

     double random_prob = distPt(rndSeed);
     double prob_bkg = candidate_deltar_prob/random_prob;
     if (prob_bkg>1) prob_bkg = 1;
//     if(ijet ==0 )cout << prob_bkg << endl;
     prob_idx[particlesReduced[ic].user_index()]*=prob_bkg;
  //   pbkg_theta[ic]=abs(prob_bkg);
    }
  //  if(ijet ==0 )cout << pbkg_theta[ic] << endl;
     prob_idx[particlesReduced[ic].user_index()] = abs(prob_idx[particlesReduced[ic].user_index()]);
  }

    if (ijet==0) {

      TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 1200, 800);
      h_deltar_constituents->Draw("HIST");
//           TGraph *correlation_ = new TGraph(particlesReduced.size(),pbkg_theta, pbkg_momentum);
//           correlation_->SetMarkerSize(1.4);
//           correlation_->SetMarkerStyle(20);
//           correlation_->Draw("AP");
           c1->SaveAs("Prueba_3.C");
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
 h_pt_constituents->Scale(1/pTbkg_estimate);

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

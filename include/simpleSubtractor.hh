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


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

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

fastjet::JetDefinition jet_estimate_bkgd(fastjet::kt_algorithm, 0.4);
fastjet::AreaDefinition area_estimate_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4)*(!fastjet::SelectorNHardest(2));
// fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4);

fastjet::ClusterSequenceArea csKt(fjInputs_, jet_estimate_bkgd, area_estimate_bkgd);
std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

 fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_estimate_bkgd, area_estimate_bkgd);

 bkgd_estimator.set_particles(fjInputs_);
 bkgd_estimator.set_jets(bkgd_jets);


rho_ = bkgd_estimator.rho();
rhoSigma_ = bkgd_estimator.sigma();


if(rho_ < 0)   rho_ = 0;

// Background constituents pT-spectrum
//----------------------------------------------------------

double maxpt_bkgd = 0;
int nparticles_bkg = 0;

for(fastjet::PseudoJet& jet : bkgd_jets) {
  std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

 // Find the maximum-pt of the background
   for(fastjet::PseudoJet p : particles) {
        double momentum = p.pt();
        if (momentum > maxpt_bkgd) maxpt_bkgd = momentum;
        nparticles_bkg++;
     }
   } // bkgd_jets loop

    int nbins_bkg = 50;
    TH1D *h = (TH1D *)gROOT->FindObject("p_{T} const");
    delete h;
    TH1D *h_pt_constituents = new TH1D("p_{T} const", "p_{T} const", nbins_bkg, 0.,maxpt_bkgd);


for(fastjet::PseudoJet& jet : bkgd_jets) {

  std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    for(fastjet::PseudoJet p : particles) {
        h_pt_constituents->Fill(p.pt());
     }
  }

  int nentries = h_pt_constituents->GetEntries();
  h_pt_constituents->Sumw2();
  h_pt_constituents->Scale(1./(double)nentries); // normalize


  double pt_binWidth = h_pt_constituents->GetBinWidth(0); // will use it later
  int nbins = h_pt_constituents->GetNbinsX();
  double ptmax = h_pt_constituents->GetBinCenter(nbins)+pt_binWidth/2;

  // Jet-by-jet subtraction
  //----------------------------------------------------------
  std::mt19937 rndSeed(rd_()); //rnd number generator seed
  std::vector<fastjet::PseudoJet> subtracted_jets;
  subtracted_jets.reserve(jets.size());
  int ijet = -1;
  int counter = jets.size();
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

    // Histogram to save the pT-information of each jet

     TH1D *h_jet = (TH1D *)gROOT->FindObject("p_{T} const jet");
     delete h_jet;
     TH1D *h_pt_constituents_jet = new TH1D("p_{T} const jet", "p_{T} const jet", nbins, 0.,maxpt_bkgd);

     double pT_reduced = 0; // pT stored by the particles below the bkg max Pt
     double patch_pT = jet.pt(); // total pT of the patch
     int nconst =0;

     for(fastjet::PseudoJet p : particles) {
       double momentum = p.pt();
       if (momentum<=maxpt_bkgd){
       h_pt_constituents_jet->Fill(momentum);
       particlesReduced.push_back(p);
       pT_reduced+=momentum;
       }
       else signalParticles.push_back(p);
       if (p.user_info<PU14>().vertex_number() == 0) nconst++;
     }
  //    if (ijet==0) cout << "Particles above the cut: " << signalParticles.size() << endl;

  //    if (ijet==0) cout << "Number of signal particles " << nconst << endl;

     h_pt_constituents_jet->Sumw2();
     int nentries_jet = h_pt_constituents_jet->GetEntries();
     h_pt_constituents_jet->Scale(1./(double)nentries_jet);

    double pTbkg_estimate = rho_*jet.area();
    //  cout << "Estimate: " << pTbkg_estimate << "pT_reduced: " << pT_reduced << endl;

    fjJetParticles_.clear();
    std::uniform_real_distribution<double> particleSelect(0., particlesReduced.size());

    // Depending on the values of the background estimate, the total pT and the pT below the bkg max Pt we have 3 options

    if (pTbkg_estimate > patch_pT)
        {
          fjJetParticles_.clear();
          median++;

        }

    else if (pTbkg_estimate > pT_reduced)
       {
      fjJetParticles_ = signalParticles;
      reduced++;
       }

     else {
      sampling++;
      std::random_shuffle(particlesReduced.begin(), particlesReduced.end()); // randomize the vector
      std::vector<fastjet::PseudoJet> notBkgParticles = signalParticles;
    //  cout << notBkgParticles.size() << endl;
      double maxPtCurrent = 0.;
      std::vector<int> avail_part(particlesReduced.size());
      std::fill(avail_part.begin(),avail_part.end(),1);
     //if(ijet==0) cout << particlesReduced.size() << endl;
      while(maxPtCurrent<pTbkg_estimate &&
      std::accumulate(avail_part.begin(),avail_part.end(),0)>0){

      int candidate = particleSelect(rndSeed);
      int ipSelected = -1;
    //  for (int i=0; i<particlesReduced.size(); i++){
    //    cout << i << endl;
        if(avail_part.at(candidate)!=0){

          double candidate_pt = particlesReduced[candidate].pt();
          int candidate_ptbin = int(candidate_pt*nbins/ptmax)+1;
          double candidate_pt_mean = h_pt_constituents->GetBinContent(candidate_ptbin);
          if(candidate_pt_mean==0) candidate_pt_mean = (h_pt_constituents->GetBinContent(candidate_ptbin+1)+h_pt_constituents->GetBinContent(candidate_ptbin-1))/2;

          double candidate_pt_error = h_pt_constituents->GetBinError(candidate_ptbin);
          double candidate_pt_prob = candidate_pt_mean + candidate_pt_error;

          double upper_limit_mean = h_pt_constituents_jet->GetBinContent(candidate_ptbin);
          double upper_limit_error = h_pt_constituents_jet->GetBinError(candidate_ptbin);
          double upper_limit = upper_limit_mean + upper_limit_error;


          if (candidate_pt_prob >= upper_limit){ //When the signal+background is below the background
           //  avail_part.at(ipSel) = 0;
             ipSelected = candidate;
             avail_part.at(ipSelected) = 0;
    //  if (ijet==0)    std::cout << " Upward Added new particle with pt = " << particlesReduced[i].pt() <<  " iparticle: " << candidate << " to init condition. total pt now " << maxPtCurrent << "/" << pTbkg_estimate  << " Particles left: " << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;
          //   break;

          }

         else{
           std::uniform_real_distribution<double> distPt(0., upper_limit);

           double random_prob = distPt(rndSeed);
          // cout << random_prob << endl;
          if (candidate_pt_prob <= upper_limit && candidate_pt_prob >= random_prob){
             ipSelected = candidate;
             avail_part.at(ipSelected) = 0;
      //      if (ijet==0)    std::cout << " Added new particle with pt = " << particlesReduced[i].pt() <<  " iparticle: " << candidate << " to init condition. total pt now " << maxPtCurrent << "/" << pTbkg_estimate  << " Particles left: " << std::accumulate(avail_part.begin(),avail_part.end(),0) << "Random probability: " << random_prob << "Candidate pT: " << candidate_pt_prob << endl;
      //    break;

         } // end if
       else continue; //if the condition is not fulfilled
       } // end else

      } // end if availability
       if(std::accumulate(avail_part.begin(),avail_part.end(),0)<5) {
         ipSelected = candidate;
         avail_part.at(ipSelected) = 0;}
    // } // end for
    //  if(ijet==0) cout << ipSelected << endl;
      if (ipSelected!=-1){
    //  ipSelected = candidate;
      fastjet::PseudoJet partSel = particlesReduced[ipSelected];
      maxPtCurrent+=partSel.pt();
     }
    } // end while

//    if (ijet==0)cout << "Signal particles after sampling: " << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;

        for (int i=0; i<particlesReduced.size(); i++){
         if(avail_part.at(i)!=0){
           fastjet::PseudoJet partSel = particlesReduced[i];
           notBkgParticles.push_back(partSel);
         }
        }

       fjJetParticles_ = notBkgParticles;

     }

  //  if (ijet==0)cout << "Estimate " << fjJetParticles_.size() << endl;

    if(fjJetParticles_.size()>0) {
    csJetSub = new fastjet::ClusterSequenceArea(fjJetParticles_, jet_defSub,    area_def_sub);
      std::vector<fastjet::PseudoJet> jetSub = fastjet::sorted_by_pt(csJetSub->inclusive_jets());
      if(jetSub[0].pt()>0.) subtracted_jets.push_back(jetSub[0]);
      if(subtracted_jets.size()>0 && subtracted_jets.size()<2)
      csJetSub->delete_self_when_unused();
    }
  } // end of jets-loop

  std::cout << "\n n subtracted jets: " << subtracted_jets.size() << "Number of jets: " << counter << "Sampling applied: " << sampling << "Median larger than pT" << median << "Reduced part too small: " << reduced << std::endl;
  return subtracted_jets;

} // end of doSubtraction


};

#endif

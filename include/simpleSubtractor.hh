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
#include "skEstimator.hh"

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

fastjet::ClusterSequenceArea csKt(fjInputs_, jet_estimate_bkgd, area_estimate_bkgd);
std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

// fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_estimate_bkgd, area_estimate_bkgd);

// bkgd_estimator.set_particles(fjInputs_);
// bkgd_estimator.set_jets(bkgd_jets);


//rho_ = bkgd_estimator.rho();
//rhoSigma_ = bkgd_estimator.sigma();


//if(rho_ < 0)   rho_ = 0;

// SoftKiller estimate for the Background pT
// -----------------------------------------------------------
//double rKTParam_ = 0.19;
//skEstimator rhoComputation(jetRParam_,rKTParam_,ghostArea_,ghostRapMax_,jetRapMax_);
//rhoComputation.setInputParticles(fjInputs_);
//std::vector<double> rho_Estimate(rhoComputation.doEstimation());

// Background constituents pT-spectrum
//----------------------------------------------------------

double maxpt_bkgd = 0;
int nparticles_bkg = 0;
double total_pt_bkgd = 0;

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

    int nbins_bkg = 25;
    TH1D *h = (TH1D *)gROOT->FindObject("p_{T} const");
    delete h;
    TH1D *h_pt_constituents = new TH1D("p_{T} const", "p_{T} const", nbins_bkg, 0.,maxpt_bkgd);

   double binWidth_bkgd = h_pt_constituents->GetBinWidth(0);

for(fastjet::PseudoJet& jet : bkgd_jets) {

  std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    for(fastjet::PseudoJet p : particles) {
        h_pt_constituents->Fill(p.pt(),p.pt());
         total_pt_bkgd+=p.pt();
     }
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

    // Histogram to save the pT-information of each jet using the area

     TH1D *h_jet = (TH1D *)gROOT->FindObject("p_{T} const jet");
     delete h_jet;
     TH1D *h_pt_constituents_jet = new TH1D("p_{T} const jet", "p_{T} const jet", nbins, 0.,maxpt_bkgd);

     double pT_reduced = 0; // pT stored by the particles below the bkg max Pt
     double patch_pT = 0; // total pT of the patch
     int nconst = 0;

     for(fastjet::PseudoJet p : particles) {
       double momentum = p.pt();
       patch_pT+=momentum;
     }


     for(fastjet::PseudoJet p : particles) {
       double momentum = p.pt();
       if (momentum<=maxpt_bkgd){
       h_pt_constituents_jet->Fill(momentum,momentum);
       particlesReduced.push_back(p);
       pT_reduced+=momentum;
       }
       else signalParticles.push_back(p);
     }

    // int nentries_jet = h_pt_constituents_jet->GetEntries();
    // double binwidth_pT = h_pt_constituents_jet->GetBinWidth(0);
    // h_pt_constituents_jet->Scale(1/((double)binwidth_pT));

  //    double rho_ = trueRho/jet.area();
  //  double pTbkg_estimate = rho_Estimate[ijet];

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
      //     cout << "After: " << h_pt_constituents->GetBinContent(1) << endl;
      // cout << pTbkg_estimate << endl;
      std::random_shuffle(particlesReduced.begin(), particlesReduced.end()); // randomize the vector
      std::vector<fastjet::PseudoJet> notBkgParticles = signalParticles;
    //  cout << notBkgParticles.size() << endl;
      double maxPtCurrent = 0.;
      std::vector<int> avail_part(particlesReduced.size());
      std::fill(avail_part.begin(),avail_part.end(),1);
      std::vector<int> selected(particlesReduced.size());
      std::vector<fastjet::PseudoJet> BkgParticles(particlesReduced.size());
     //if(ijet==0) cout << particlesReduced.size() << endl;
  //   cout << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;
      int trials = 0;
      while(maxPtCurrent<=pTbkg_estimate &&
      std::accumulate(avail_part.begin(),avail_part.end(),0)>1){
      trials++;
    //  if (ijet==123) cout << "eo" << endl;
      int candidate = particleSelect(rndSeed);
      int ipSelected = -1;
    //   if(ijet==0) cout << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;
    //  if(std::accumulate(avail_part.begin(),avail_part.end(),0)==1) cout << "Jet" <<  ijet << "Only one" << candidate << endl;
      //  if (ijet==4) cout << candidate << endl;
    //  for (int i=0; i<particlesReduced.size(); i++){
    //    cout << i << endl;
        if(avail_part.at(candidate)!=0){

          double candidate_pt = particlesReduced[candidate].pt();
          int candidate_ptbin = int(candidate_pt*nbins/ptmax)+1;
          double candidate_pt_mean = h_pt_constituents->GetBinContent(candidate_ptbin);
          if(candidate_pt_mean==0) candidate_pt_mean = (h_pt_constituents->GetBinContent(candidate_ptbin+1)+h_pt_constituents->GetBinContent(candidate_ptbin-1))/2;
      //    if(candidate_pt_mean==0) candidate_pt_mean = (h_pt_constituents->GetBinContent(candidate_ptbin+2)+h_pt_constituents->GetBinContent(candidate_ptbin-2))/2;

          double candidate_pt_error = h_pt_constituents->GetBinError(candidate_ptbin);
          double candidate_pt_prob = candidate_pt_mean + candidate_pt_error;

          double upper_limit_mean = h_pt_constituents_jet->GetBinContent(candidate_ptbin);
          double upper_limit_error = h_pt_constituents_jet->GetBinError(candidate_ptbin);
          double upper_limit = upper_limit_mean + upper_limit_error;

        //  cout << upper_limit_mean << endl;
      //    if (candidate_pt_mean == 0){
      //       cout << "eo" << endl;
          //   cout << pTbkg_estimate << endl;
      //     }
          //  cout << "yes" << endl;
          //  cout << candidate_pt << endl;
          //  cout << candidate_ptbin << endl;
          //}

          if (candidate_pt_prob >= upper_limit){ //When the signal+background is below the background

             ipSelected = candidate;
             avail_part.at(ipSelected) = 0;

        //     if(ijet==4) cout << "Particle accepted: " << candidate << endl;

          }

         else{
           std::uniform_real_distribution<double> distPt(0., upper_limit);

           double random_prob = distPt(rndSeed);
      //      if(ijet==0)  cout << "Bkg: " << candidate_pt_prob << "Bkg+Sig" << random_prob << endl;
          if (candidate_pt_prob >= random_prob){
             ipSelected = candidate;
             avail_part.at(ipSelected) = 0;
          //   if(ijet==4) cout << "Particle accepted: " << candidate << endl;

         } // end if
       else continue; //if the condition is not fulfilled
       } // end else

      } // end if availability

      if (ipSelected!=-1){

      fastjet::PseudoJet partSel = particlesReduced[ipSelected];
      maxPtCurrent+=partSel.pt();
      BkgParticles.push_back(partSel);
      selected.push_back(ipSelected);
     }
//   if(trials>10000){cout << "Jet: "<< ijet << "Number of particles available: " << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;}
//  if(ijet==4)    cout << "Jet: "<< ijet << "Number of particles available: " << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;
   //}
    } // end while
    h_pt_constituents->Scale(1/pTbkg_estimate);

// Complete the list

       if (maxPtCurrent<pTbkg_estimate && std::accumulate(avail_part.begin(),avail_part.end(),0)>0){
         while(maxPtCurrent<pTbkg_estimate && std::accumulate(avail_part.begin(),avail_part.end(),0)>0){
         for (int i=0; i<particlesReduced.size(); i++){
          if(avail_part.at(i)!=0){
            int ipSelected = i;
            avail_part.at(ipSelected) = 0;
            fastjet::PseudoJet partSel = particlesReduced[ipSelected];
            maxPtCurrent+=partSel.pt();
            BkgParticles.push_back(partSel);
            selected.push_back(ipSelected);
            break;
            }
         }
       }
     }

     if(maxPtCurrent>pTbkg_estimate && selected.size()>1) {
    double maxPtPrev = 0;
    int size = selected.size()-1;
    int last_particle = selected.at(size);
//  if(ijet==1) cout << last_particle << endl;
   int size_bkg = BkgParticles.size()-1;
   fastjet::PseudoJet partSel = BkgParticles[size_bkg];
 // cout << size << endl;
   maxPtPrev = maxPtCurrent - partSel.pt();
//  if(ijet==1)  cout << maxPtPrev << endl;
    double distance_one = sqrt((pTbkg_estimate-maxPtCurrent)*(pTbkg_estimate-maxPtCurrent));
   double distance_two = sqrt((pTbkg_estimate-maxPtPrev)*(pTbkg_estimate-maxPtPrev));

  if (distance_one > distance_two) {
     avail_part.at(last_particle) = 1; // convert it into signal particle
     maxPtCurrent=maxPtPrev;
   }
 }


  //   cout << "Truth: " << trueRho << "JetPt" << patch_pT << "Subtracted: " << maxPtCurrent << endl;

//cout << ijet << endl;
        for (int i=0; i<particlesReduced.size(); i++){
         if(avail_part.at(i)!=0){
           fastjet::PseudoJet partSel = particlesReduced[i];
           notBkgParticles.push_back(partSel);
         }
        }

       fjJetParticles_ = notBkgParticles;
//       if(ijet==0) cout << fjJetParticles_.size() << endl;

     }

    if(fjJetParticles_.size()>0) {
    csJetSub = new fastjet::ClusterSequenceArea(fjJetParticles_, jet_defSub,    area_def_sub);
      std::vector<fastjet::PseudoJet> jetSub = fastjet::sorted_by_pt(csJetSub->inclusive_jets());
      if(jetSub[0].pt()>0.) subtracted_jets.push_back(jetSub[0]);
      if(subtracted_jets.size()>0 && subtracted_jets.size()<2)
      csJetSub->delete_self_when_unused();
    }
  } // end of jets-loop

//  std::cout << "\n n subtracted jets: " << subtracted_jets.size() << "Number of jets: " << counter << "Sampling applied: " << sampling << "Median larger than pT" << median << "Reduced part too small: " << reduced << std::endl;
  return subtracted_jets;

} // end of doSubtraction


};

#endif

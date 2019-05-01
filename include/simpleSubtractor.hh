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

fastjet::JetDefinition jet_defSub(antikt_algorithm, 25.);
fastjet::AreaDefinition area_def_sub = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec_sub);
fastjet::ClusterSequenceArea csSub(fjInputs_, jet_defSub, area_def_sub);

// create what we need for the background estimation
//----------------------------------------------------------

fastjet::JetDefinition jet_estimate_bkgd(fastjet::kt_algorithm, 0.4);
fastjet::AreaDefinition area_estimate_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4)*(!fastjet::SelectorNHardest(2));

fastjet::ClusterSequenceArea csKt(fjInputs_, jet_estimate_bkgd, area_estimate_bkgd);
std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

 fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_estimate_bkgd, area_estimate_bkgd);

 bkgd_estimator.set_particles(fjInputs_);
 bkgd_estimator.set_jets(bkgd_jets);


rho_ = bkgd_estimator.rho();
rhoSigma_ = bkgd_estimator.sigma();


if(rho_ < 0)   rho_ = 0;

// SoftKiller estimate for the Background pT
// -----------------------------------------------------------
//double rKTParam_ = 0.19;
//skEstimator rhoComputation(jetRParam_,rKTParam_,ghostArea_,ghostRapMax_,jetRapMax_);
//rhoComputation.setInputParticles(fjInputs_);
//vector<double> rho_Estimate(rhoComputation.doEstimation());

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

    int nbins_bkg = 10;
    TH1D *h = (TH1D *)gROOT->FindObject("background");
    delete h;
    TH1D *h_pt_constituents_bkg = new TH1D("background", "background", nbins_bkg, 0.,maxpt_bkgd);
    double bin = h_pt_constituents_bkg->GetBinWidth(0);


for(fastjet::PseudoJet& jet : bkgd_jets) {

  std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
    fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
    for(fastjet::PseudoJet p : particles) {
        h_pt_constituents_bkg->Fill(p.pt(),p.pt()/bin);
        total_pt_bkgd+=p.pt();
     }
  //   total_pt_bkgd+=jet.pt();
  }

  h_pt_constituents_bkg->Scale(1./((double)total_pt_bkgd)); //normalize


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


     double pT_reduced = 0; // pT stored by the particles below the bkg max Pt
     int nconst = 0;
     double patch_pT=0;

     for(fastjet::PseudoJet p : particles) {
       double momentum = p.pt();
       patch_pT+=momentum;
       if (momentum<=maxpt_bkgd){
       particlesReduced.push_back(p);
       pT_reduced+=momentum;
       }
       else signalParticles.push_back(p);
     }

   double trueRho = 0;
    for(fastjet::PseudoJet p : particles) {
    if (abs(p.user_info<PU14>().vertex_number()) == 1) trueRho+=p.pt();
   }

  double pTbkg_estimate = trueRho;
//     double rho_ = trueRho/jet.area();

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
    //  h_pt_constituents_jet->Scale(pT_reduced);
  //    double area = h_pt_constituents_bkg->Integral(0,100,"width");
      h_pt_constituents_bkg->Scale(pTbkg_estimate);
    //  if(ijet==0) cout << h_pt_constituents_bkg->Integral(0,100,"width") << endl;
      std::random_shuffle(particlesReduced.begin(), particlesReduced.end()); // randomize the vector
      std::vector<fastjet::PseudoJet> notBkgParticles = signalParticles;
      std::vector<fastjet::PseudoJet> BkgParticles;
      BkgParticles.reserve(particlesReduced.size());
    //  if(ijet==0)   cout << "Number of particles inicial: " << notBkgParticles.size() << endl;
    //  cout << notBkgParticles.size() << endl;

    // Empty histogram to save the pT-information of each jet

     TH1D *h_jet = (TH1D *)gROOT->FindObject("estimate");
     delete h_jet;
     TH1D *h_pt_constituents_estimate = new TH1D("estimate", "estimate", nbins_bkg, 0.,maxpt_bkgd);
      double binWidth=h_pt_constituents_estimate->GetBinWidth(0);
      double maxPtCurrent = 0.;
      std::vector<int> avail_part(particlesReduced.size());
      std::fill(avail_part.begin(),avail_part.end(),1);
      std::vector<int> selected(particlesReduced.size());
      while(std::accumulate(avail_part.begin(),avail_part.end(),0)>0 && maxPtCurrent < pTbkg_estimate){
      // Choose a random particle from the list
      int candidate = particleSelect(rndSeed);
      int ipSelected = -1;
      // Find the bin to which it corresponds
      if (avail_part.at(candidate)!=0){
      double candidate_pt = particlesReduced[candidate].pt();
      int candidate_ptbin = int(candidate_pt*nbins_bkg/maxpt_bkgd)+1;
    //  if(ijet==0) cout << h_pt_constituents_estimate->GetBinContent(8) << endl;
      if (h_pt_constituents_estimate->GetBinContent(candidate_ptbin) < h_pt_constituents_bkg->GetBinContent(candidate_ptbin) && h_pt_constituents_bkg->GetBinContent(candidate_ptbin) > h_pt_constituents_bkg->GetBinCenter(candidate_ptbin)-binWidth/2){
    //    cout << "eeo" << endl;
      h_pt_constituents_estimate->Fill(candidate_pt,candidate_pt/binWidth);
      avail_part.at(candidate) = 0;
      maxPtCurrent+=candidate_pt;
      fastjet::PseudoJet partSel = particlesReduced[candidate];
      BkgParticles.push_back(partSel);
      selected.push_back(candidate);
    //  if(ijet==0) cout << candidate << endl;

      }
      else{
        fastjet::PseudoJet partSel = particlesReduced[candidate];
        notBkgParticles.push_back(partSel);
        avail_part.at(candidate) = 0;
        selected.push_back(candidate);
      }
    //  if(ijet==0) cout << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;
    }
  }
//  if(ijet==0){
//    TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 1200, 800);

//          h_pt_constituents_bkg->Draw("HIST");
//          h_pt_constituents_estimate->Draw("SAME HIST");
//        c1->SaveAs("Prueba_sampler.C");
//  }
//  if(ijet==1) cout << maxPtCurrent << endl;
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
// cout << "TotalPt: " << patch_pT << "Estimate: " << pTbkg_estimate << "Subtracted:" << maxPtCurrent << "Integral: " << h_pt_constituents_estimate->Integral(0,100,"width") << endl;
    //  if(ijet==1) cout << maxPtCurrent << endl;

      //       cout << std::accumulate(avail_part.begin(),avail_part.end(),0) << endl;
    for (int i=0; i<particlesReduced.size(); i++){
         if(avail_part.at(i)!=0){
           fastjet::PseudoJet partSel = particlesReduced[i];
           notBkgParticles.push_back(partSel);
         }
        }

       fjJetParticles_ = notBkgParticles;
  //  if(ijet==0) cout << "Number of particles final: "  << fjJetParticles_.size() << endl;

     }


    if(fjJetParticles_.size()>0) {
    csJetSub = new fastjet::ClusterSequenceArea(fjJetParticles_, jet_defSub,    area_def_sub);
      std::vector<fastjet::PseudoJet> jetSub = fastjet::sorted_by_pt(csJetSub->inclusive_jets());
    //   cout << jetSub[0].pt() << endl;
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

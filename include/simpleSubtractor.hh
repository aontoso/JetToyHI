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
#include "randomCones.hh"

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

fastjet::JetDefinition jet_estimate_bkgd(fastjet::antikt_algorithm, 0.4);
fastjet::AreaDefinition area_estimate_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4)*(!fastjet::SelectorNHardest(2));

fastjet::ClusterSequenceArea csKt(fjInputs_, jet_estimate_bkgd, area_estimate_bkgd);
std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

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

  // double binWidth_bkgd = h_pt_constituents->GetBinWidth(0);

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
  int sampling = 0;
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

      double cone_rapidity = -jet.eta()+jetRParam_;
      randomCones rigidcone(2, jetRParam_, cone_rapidity);
      rigidcone.setInputParticles(fjInputs_);
      std::vector<std::vector<double>> cones;
      cones = rigidcone.run();

  // Random cones constituents pT-spectrum
  //----------------------------------------------------------

  double maxdeltar_bkgd = 0;
  int nparticles_bkg = 0;
  //double total_deltaR_bkgd = 0;

  for(int i = 0; i<cones.size(); i++) {
   // Find the maximum-pt of the background
     for(int j = 0; j<cones[i].size(); j++) {
          double deltar = cones[i][j];
          if (deltar > maxdeltar_bkgd) maxdeltar_bkgd = deltar;
          nparticles_bkg++;
       }
     } // random_cone loop
//
      int nbins_bkg_deltar = 5;
      TH1D *h_cone = (TH1D *)gROOT->FindObject("deltaR const");
      delete h_cone;
      TH1D *h_deltar_constituents = new TH1D("deltaR const", "deltaR const", nbins_bkg_deltar, 0.,maxdeltar_bkgd);


  for(int i = 0; i<cones.size(); i++) {
    for(int j = 0; j<cones[i].size(); j++) {
         double deltar = cones[i][j];
         h_deltar_constituents->Fill(deltar);
      }
    }

    //      cout << maxpt_bkgd << endl;
    h_deltar_constituents->Scale(1./((double)h_deltar_constituents->GetEntries())); //normalize

    double deltar_binWidth = h_deltar_constituents->GetBinWidth(0); // will use it later
    int nbins_deltar = h_deltar_constituents->GetNbinsX();
    double deltarmax = h_deltar_constituents->GetBinCenter(nbins_deltar)+deltar_binWidth/2;

    //Jet constituents pT-spectrum
    //----------------------------------------------------------

    // Create a new list of particles whose pT is below the bkg max Pt
     std::vector<fastjet::PseudoJet> particlesReduced;
    // Create a new list of particles whose pT is above the bkg max Pt and, therefore, considered to be signal
     std::vector<fastjet::PseudoJet> signalParticles;

    // Histogram to save the pT-information of each jet

     TH1D *h_jet = (TH1D *)gROOT->FindObject("p_{T} const jet");
     delete h_jet;
     TH1D *h_pt_constituents_jet = new TH1D("p_{T} const jet", "p_{T} const jet", nbins, 0.,maxpt_bkgd);

     double pT_reduced = 0; // pT stored by the particles below the bkg max Pt
     double patch_pT = 0; // total pT of the patch
     int nconst = 0;

     TH1D *h_deltar = (TH1D *)gROOT->FindObject("deltaR const jet");
     delete h_deltar;
     TH1D *h_deltar_constituents_jet = new TH1D("deltaR const jet", "deltaR const jet", nbins_deltar, 0.,maxdeltar_bkgd);

     for(fastjet::PseudoJet p : particles) {
       double momentum = p.pt();
       double deltar = jet.delta_R(p);
       patch_pT+=momentum;
       if (momentum<=maxpt_bkgd && deltar<maxdeltar_bkgd){
       h_pt_constituents_jet->Fill(momentum,momentum);
       h_deltar_constituents_jet->Fill(deltar);
       particlesReduced.push_back(p);
       pT_reduced+=momentum;
       }
       else signalParticles.push_back(p);
     }

    h_deltar_constituents_jet->Scale(1./((double)h_deltar_constituents_jet->GetEntries()));



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

       }

   else if (pTbkg_estimate >= pT_reduced)
      {
     fjJetParticles_ = signalParticles;
     reduced++;
      }

    else {
     sampling++;
       h_pt_constituents->Scale(pTbkg_estimate);

    //   if(ijet==0){
    //  TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 1200, 800);

    //        h_deltar_constituents->Draw("HIST");
    //        h_deltar_constituents_jet->DrawNormalized(" HIST SAME");
    //        c1->SaveAs("Prueba_4.C");
    //      }

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

     int candidate = particleSelect(rndSeed);
     int ipSelected = -1;

       if(avail_part.at(candidate)!=0){

         double candidate_pt = particlesReduced[candidate].pt();
         double candidate_deltar=jet.delta_R(particlesReduced[candidate]);
      //   int candidate_ptbin = int(candidate_pt*nbins/ptmax)+1;
         int candidate_deltarbin = int(candidate_deltar*nbins_deltar/deltarmax)+1;
         double candidate_deltar_mean = h_deltar_constituents->GetBinContent(candidate_deltarbin);
        // double candidate_pt_mean = h_pt_constituents->GetBinContent(candidate_ptbin);
      //   if(candidate_pt_mean==0) candidate_pt_mean = (h_pt_constituents->GetBinContent(candidate_ptbin+1)+h_pt_constituents->GetBinContent(candidate_ptbin-1))/2;
      if(candidate_deltar_mean==0) candidate_deltar_mean = (h_deltar_constituents->GetBinContent(candidate_deltarbin+1)+h_deltar_constituents->GetBinContent(candidate_deltarbin-1))/2;


    //     double candidate_pt_error = h_pt_constituents->GetBinError(candidate_ptbin);
         double candidate_deltar_error = h_deltar_constituents->GetBinError(candidate_deltarbin);
    //     double candidate_pt_prob = candidate_pt_mean + candidate_pt_error;
         double candidate_deltar_prob = candidate_deltar_mean + candidate_deltar_error;

    //     double upper_limit_mean = h_pt_constituents_jet->GetBinContent(candidate_ptbin);
         double upper_limit_mean = h_deltar_constituents_jet->GetBinContent(candidate_deltarbin);
      //   double upper_limit_error = h_pt_constituents_jet->GetBinError(candidate_ptbin);
         double upper_limit_error = h_deltar_constituents_jet->GetBinError(candidate_deltarbin);

         double upper_limit = upper_limit_mean + upper_limit_error;



         if (candidate_deltar_prob >= upper_limit){ //When the signal+background is below the background

            ipSelected = candidate;
            avail_part.at(ipSelected) = 0;

       //     if(ijet==4) cout << "Particle accepted: " << candidate << endl;

         }

        else{
          std::uniform_real_distribution<double> distPt(0., upper_limit);

          double random_prob = distPt(rndSeed);
      //   if(ijet == 0) cout << "Bkg: " << candidate_deltar_prob << "Bkg+Sig" << random_prob << endl;
         if (candidate_deltar_prob >= random_prob){
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

#ifndef sharedLayerSubtractor_h
#define sharedLayerSubtractor_h

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

#include "Angularity.hh"
#include "rhoEstimator.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the jet-by-jet shared layer subtraction
// Author: A. Soto-Ontoso, M. Verweij, Y. Mehtar-Tani
//---------------------------------------------------------------

class sharedLayerSubtractor {

private :
  double jetRParam_;
  double ghostArea_;
  double ghostRapMax_;
  double jetRapMax_;
  double rho_;
  double rhoSigma_;
  double pTDbkg_;
  double pTDbkgSigma_;
  double massBkg_;
  double massBkgSigma_;


  Angularity pTD_;
  Angularity mass_;

  std::vector<fastjet::PseudoJet> fjInputs_;
  std::vector<fastjet::PseudoJet> fjJetInputs_;

  std::vector<fastjet::PseudoJet> fjJetParticles_;
  std::vector<std::vector<double>> fChi2s_;
  std::vector<std::vector<int>> fShare_;

  std::random_device rd_;
  int nInitCond_;
  int nTopInit_;

  fastjet::ClusterSequenceArea *csJetSub;

public :
  sharedLayerSubtractor(double rJet = 0.4,
                        double ghostArea = 0.001,
                        double ghostRapMax = 3.0,
                        double jetRapMax = 3.0,
                        int nInitCond = 1.,
                        int nTopInit = 1.) :
    jetRParam_(rJet),
    ghostArea_(ghostArea),
    ghostRapMax_(ghostRapMax),
    jetRapMax_(jetRapMax),
    nInitCond_(nInitCond),
    nTopInit_(nTopInit)
  {
    pTD_ = Angularity(0.,2.,0.4);
    mass_ = Angularity(2., 1., 0.4);

  }

  void setGhostArea(double a) { ghostArea_ = a; }

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }
  void setInputJets(std::vector<fastjet::PseudoJet> v)      { fjJetInputs_ = v; }

  double getRho()  const { return rho_; }
  double getRhoSigma() const { return rhoSigma_; }

  double getPTDBkg() const { return pTDbkg_; }
  double getPTDBkgSigma() const { return pTDbkgSigma_; }

  std::vector<std::vector<double>> getChi2s() const { return fChi2s_; }
  std::vector<std::vector<int>> getNShared() const { return fShare_; }
  void clearMemory() {
  fChi2s_.clear();
  fShare_.clear();
  if(csJetSub) delete csJetSub;
}

  std::vector<fastjet::PseudoJet> doSubtraction() {


      fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, 0.001);
    fastjet::GhostedAreaSpec ghost_spec_sub(ghostRapMax_, 1, 1.);
    std::vector<fastjet::PseudoJet> jets = fjJetInputs_;

    // Rho-correction estimate for the Background pT
   // -----------------------------------------------------------

   rhoEstimator rhoComputation(jetRParam_,ghostArea_,ghostRapMax_,jetRapMax_);
   rhoComputation.setInputParticles(fjInputs_);
   vector<double> rho_Estimate(rhoComputation.doEstimation());

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
    //UE metric
    //----------------------------------------------------------
    std::vector<double> pTD_bkgd;
    std::vector<double> mass_bkgd;
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

        int nbins_bkg = int(sqrt(nparticles_bkg));
        TH1D *h = (TH1D *)gROOT->FindObject("p_{T} const");
        delete h;
        TH1D *h_pt_constituents = new TH1D("p_{T} const", "p_{T} const", nbins_bkg, 0.,maxpt_bkgd);


    for(fastjet::PseudoJet& jet : bkgd_jets) {

      pTD_bkgd.push_back(pTD_.result(jet));
      mass_bkgd.push_back(mass_.result(jet));
      std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
        fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

        for(fastjet::PseudoJet p : particles) {
            h_pt_constituents->Fill(p.pt());
         }
      }

      int nentries = h_pt_constituents->GetEntries();
      h_pt_constituents->Sumw2();
      h_pt_constituents->Scale(1./(double)nentries); // normalize
      h_pt_constituents->Draw();
    //  exit(0);
      double pt_binWidth = h_pt_constituents->GetBinWidth(0); // will use it later
      int nbins = h_pt_constituents->GetNbinsX();
      double ptmax = h_pt_constituents->GetBinCenter(nbins)+pt_binWidth/2;

     std::nth_element(pTD_bkgd.begin(), pTD_bkgd.begin() + pTD_bkgd.size()/2, pTD_bkgd.end());
     double med_pTD = pTD_bkgd[pTD_bkgd.size()/2];

     std::nth_element(mass_bkgd.begin(), mass_bkgd.begin() + mass_bkgd.size()/2, mass_bkgd.end());
     double med_mass = mass_bkgd[mass_bkgd.size()/2];


     int nRMS = 0;
     double rms_pTD = 0.;
     for(int ip = 0; ip<(int)pTD_bkgd.size(); ++ip) {
      rms_pTD += (pTD_bkgd[ip]-med_pTD)*(pTD_bkgd[ip]-med_pTD);
      nRMS++;
     }
     if(nRMS>0.)
      rms_pTD = sqrt(rms_pTD/(double)nRMS);

      int nRMSm = 0;
      double rms_mass = 0.;
      for(int ip = 0; ip<(int)mass_bkgd.size(); ++ip) {
        rms_mass += (mass_bkgd[ip]-med_mass)*(mass_bkgd[ip]-med_mass);
        nRMSm++;
      }
      if(nRMSm>0.)
        rms_mass = sqrt(rms_mass/(double)nRMSm);

     pTDbkg_ = med_pTD;
     pTDbkgSigma_ = rms_pTD;

     massBkg_ = med_mass;
     massBkgSigma_ = rms_mass;

     std::mt19937 rndSeed(rd_()); //rnd number generator seed

     std::vector<fastjet::PseudoJet> subtracted_jets;
     subtracted_jets.reserve(jets.size());
     int ijet = -1;


    for(fastjet::PseudoJet& jet : jets) {
      ++ijet;
      if(jet.is_pure_ghost()) continue;

      //get ghosts and true particles (ghosts are distributed uniformly which we will use to create initial conditions)
      //----------------------------------------------------------
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
      if(particles.size()<1 || jet.pt()<1.) continue;
      std::random_shuffle(particles.begin(), particles.end()); // randomize the vector
      std::random_shuffle(ghosts.begin(),ghosts.end());

      // Use the background+signal information in each patch to do the pt-selection

      TH1D *h_jet = (TH1D *)gROOT->FindObject("p_{T} const jet");
      delete h_jet;
      TH1D *h_pt_constituents_jet = new TH1D("p_{T} const jet", "p_{T} const jet", nbins_bkg, 0.,maxpt_bkgd);
      double trueRho = 0;
    for(fastjet::PseudoJet p : particles) {
      double momentum = p.pt();
      if (momentum<=maxpt_bkgd){
      h_pt_constituents_jet->Fill(momentum);
      }
      if (p.user_info<PU14>().vertex_number() == 1) trueRho+=p.pt();
     }

     h_pt_constituents_jet->Sumw2();
     int nentries_jet = h_pt_constituents_jet->GetEntries();
     h_pt_constituents_jet->Scale(1./(double)nentries_jet);

      //set user_index of all particles to position particles vector
      for(int i = 0; i<(int)particles.size(); ++i) {
        particles[i].set_user_index(i);
      }

      // create map of ghost and closest particle
      std::vector<std::vector<int>> closestPartToGhost;
      for(int ighost = 0; ighost<ghosts.size(); ++ighost) {
        std::vector<int> ipSel = findClosestParticles(particles, ighost, ghosts, 5);
        closestPartToGhost.push_back(ipSel);
       }

      // create requested number of initial conditions
      //----------------------------------------------------------
      std::vector<std::vector<int>> collInitCond;
      for(int ii = 0; ii<nInitCond_; ++ii) {
        std::uniform_int_distribution<> distUni(0,ghosts.size()); //uniform distribution of ghosts in vector
        std::vector<int> initCondition;                           //list of particles in initial condition

      //  double maxPt = trueRho;
          double maxPt = rho_Estimate.at(ijet)-12.92;
      //    cout << rho_*jet.area() << endl;
      //    cout << maxPt << endl;
        //  cout <<trueRho << endl;
      //  cout <<maxPt <<endl;
      //  int rejection = 0;
        //make copy of particles so that a particle is not repeated inside the same initial condition
        std::vector<fastjet::PseudoJet> particlesNotUsed = particles;

        std::vector<std::vector<int>> closestPartToGhostNotUsed = closestPartToGhost;

        double maxPtCurrent = 0.;
        std::vector<int> avail(closestPartToGhost.size());
        std::vector<int> avail_part(particlesNotUsed.size());
        std::vector<int> part_accepted(particlesNotUsed.size());
        std::fill(avail.begin(),avail.end(),1);
        std::fill(avail_part.begin(),avail_part.end(),1);
        std::fill(part_accepted.begin(),part_accepted.end(),1);
      //  maxpt_bkgd = 3;
      //  int cutoff = int(particles.size());
      //  int cutoff_rejection = 100;
         while(maxPtCurrent<maxPt &&
         std::accumulate(avail_part.begin(),avail_part.end(),0)>0
         && std::accumulate(avail.begin(),avail.end(),0)>0){
          //cout << rejection << endl;
          //pick random ghost

        //      if (ijet==8 && std::accumulate(avail_part.begin(),avail_part.end(),0)<8) std::cout<< "iparticle: " << iparticle << avail_part.at(iparticle) << endl;
    //      int iparticle = -1;
    //pick random ghost
        int ighost = int(std::floor(distUni(rndSeed)));
        if(ighost>=ghosts.size()) continue;

        int iparticle = -1;
        if(closestPartToGhostNotUsed[ighost].size()>0) {
           iparticle = closestPartToGhostNotUsed[ighost][0];
        if(closestPartToGhostNotUsed[ighost].size()<2) avail[ighost] = 0;
        closestPartToGhostNotUsed[ighost].erase(closestPartToGhostNotUsed[ighost].begin()+0);
        } else
        continue;

     int ipSel = 0; // make sure you do not repeat a particle
     if (avail_part.at(iparticle)!=0){
      ipSel = iparticle;
      avail_part.at(iparticle) = 0;
      }
      else continue;
            // cout << rejection << endl;ontinue;
           //};

          fastjet::PseudoJet partSel = particlesNotUsed[ipSel];

          //---------------------------------------------------------------
          // Compare the particle pt with the average background spectrum
          //---------------------------------------------------------------

            double candidate_pt = partSel.pt();
            int candidate_ptbin = int(candidate_pt*nbins/ptmax)+1;

            if (candidate_pt > maxpt_bkgd) continue;
            double candidate_pt_mean = h_pt_constituents->GetBinContent(candidate_ptbin);
            double candidate_pt_error = h_pt_constituents->GetBinError(candidate_ptbin);
            double candidate_pt_prob = candidate_pt_mean + candidate_pt_error;

            // Do it on a jet-by-jet basis

               double upper_limit_mean = h_pt_constituents_jet->GetBinContent(candidate_ptbin);
               double upper_limit_error = h_pt_constituents_jet->GetBinError(candidate_ptbin);
               double upper_limit = upper_limit_mean + upper_limit_error;

               if (upper_limit <= candidate_pt_prob){ //When the signal+background is below the background
                //  avail_part.at(ipSel) = 0;
                  part_accepted.at(ipSel) = 0;
                  initCondition.push_back(partSel.user_index());
                  maxPtCurrent+=partSel.pt();
        //         if (ijet==0 && ii==0)    std::cout << "Downward Added new particle with pt = " << partSel.phi() <<  " iparticle: " << ipSel << " to init condition. total pt now " << maxPtCurrent << "/" << maxPt  << " Particles left: " << std::accumulate(part_accepted.begin(),part_accepted.end(),0) << endl;
               }
               else{

               std::uniform_real_distribution<double> distPt(0., upper_limit);

               double random_prob = distPt(rndSeed);
          //  cout << random_prob << endl;
              if (candidate_pt_prob > random_prob){
                 part_accepted.at(ipSel) = 0;
                 initCondition.push_back(partSel.user_index());
                 maxPtCurrent+=partSel.pt();
        //  if (ijet==0 && ii==0)    std::cout << "Added new particle with pt = " << partSel.pt() <<  " iparticle: " << ipSel << " to init condition. total pt now " << maxPtCurrent << "/" << maxPt  << " Particles left: " << std::accumulate(part_accepted.begin(),part_accepted.end(),0) << endl;
            } else continue;
           } continue;
        } // while loop

        //Complete the list

        if (maxPtCurrent<maxPt &&
           std::accumulate(part_accepted.begin(),part_accepted.end(),0)>0){
          while(maxPtCurrent<maxPt){
            int ipSelected = 0;
            int candidate = 0;
      //    if (ijet==0 && ii==0)  cout << particles.size() << endl;
            for (int i=0; i<particles.size(); i++){
        //    if (ijet==0 && ii==0)  cout << "Particula:"<< part_accepted.at(i)<< "pt " << particles[i].pt()<< endl;
              if(part_accepted.at(i)!=0 && particles[i].pt()<maxpt_bkgd)
               { double candidate_pt = particles[i].pt();
                 int candidate_ptbin = int(candidate_pt*nbins/ptmax)+1;
                 double candidate_pt_mean = h_pt_constituents->GetBinContent(candidate_ptbin);
                 double candidate_pt_error = h_pt_constituents->GetBinError(candidate_ptbin);
                 double candidate_pt_prob = candidate_pt_mean + candidate_pt_error;

                 double upper_limit_mean = h_pt_constituents_jet->GetBinContent(candidate_ptbin);
                 double upper_limit_error = h_pt_constituents_jet->GetBinError(candidate_ptbin);
                 double upper_limit = upper_limit_mean + upper_limit_error;

                 if (upper_limit <= candidate_pt_prob){ //When the signal+background is below the background
                  //  avail_part.at(ipSel) = 0;
                    candidate=i;
                    part_accepted.at(i) = 0;
                    break;
          //         if (ijet==0 && ii==0)    std::cout << "Downward Added new particle with pt = " << partSel.phi() <<  " iparticle: " << ipSel << " to init condition. total pt now " << maxPtCurrent << "/" << maxPt  << " Particles left: " << std::accumulate(part_accepted.begin(),part_accepted.end(),0) << endl;
                 }
               else{
                 std::uniform_real_distribution<double> distPt(0., upper_limit);

                 double random_prob = distPt(rndSeed);
            //  cout << random_prob << endl;
                if (candidate_pt_prob >= random_prob){
                  candidate = i;
                   part_accepted.at(i) = 0;
                break;
          //  if (ijet==0 && ii==0)    std::cout << "Added new particle with pt = " << partSel.pt() <<  " iparticle: " << ipSel << " to init condition. total pt now " << maxPtCurrent << "/" << maxPt  << " Particles left: " << std::accumulate(part_accepted.begin(),part_accepted.end(),0) << endl;
               }
              else continue;
            }
          }
          //    continue;}
                else continue;
              }

              ipSelected = candidate;
              fastjet::PseudoJet partSel = particlesNotUsed[ipSelected];
              maxPtCurrent+=partSel.pt();
      //        if (ijet==0 && ii==0)    std::cout << " Second way Added new particle with pt = " << partSel.phi() <<  " iparticle: " << ipSel << " to init condition. total pt now " << maxPtCurrent << "/" << maxPt  << " Particles left: " << std::accumulate(part_accepted.begin(),part_accepted.end(),0) << endl;
              initCondition.push_back(partSel.user_index());
          }
        }
        if(maxPtCurrent>maxPt) {
         int initConditionSize_ = initCondition.size();
         double maxPtPrev = 0;
         std::vector<fastjet::PseudoJet> initCondParticles;
          for(int ic = 0; ic<initConditionSize_; ++ic) {
            initCondParticles.push_back(particles[initCondition[ic]]);
            maxPtPrev+=initCondParticles.at(ic).pt();
          }
          double distance_one = sqrt((maxPt-maxPtCurrent)*(maxPt-maxPtCurrent));
          double distance_two = sqrt((maxPt-maxPtPrev)*(maxPt-maxPtPrev));

          if (distance_one > distance_two) {
        //    if(ijet==0) cout<<"yes"<< maxPtPrev<<endl;
            initCondition.pop_back();

          }
          collInitCond.push_back(initCondition); //avoid putting in a initial condition for which not enough particles were available anymore to get to the required pT. Might be an issue for sparse events.
        //}
        }
      //  cout << "ijet: " << ijet << "MaxPt: " << maxPtCurrent << endl;
      }//initial conditions loop
      //  cout << ijet << endl;
      //----------------------------------------------------------
      //Now we have the requested number of random initial condition

      //Next step: calc chi2 for each initial condition
      //----------------------------------------------------------
       std::vector<double> chi2s = calculateChi2s(collInitCond, particles, med_pTD, rms_pTD); //

       fChi2s_.push_back(chi2s);

      //sort the chi2s keeping track of indices
     //----------------------------------------------------------
     // initialize original index locations
     std::vector<size_t> idx(chi2s.size());
     iota(idx.begin(), idx.end(), 0);

    //  sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(),
               [&chi2s](size_t i1, size_t i2) {return chi2s[i1] < chi2s[i2];});


      //Next step: figure out how often each particle is shared in nTopInit_ initial conditions
      //----------------------------------------------------------
     std::vector<int> share_idx(particles.size(),0);
    for(int it = 0; it<std::min(nTopInit_,(int)collInitCond.size()); ++it) {
       int chi2Index = idx[it];
       std::vector<int> indices = collInitCond[chi2Index];
        for(int ic = 0; ic<(int)indices.size(); ++ic) {
      //    if (ijet==0)cout << particles[indices[ic]].user_index() << endl;
         share_idx[particles[indices[ic]].user_index()]++;

       }
     }

     fShare_.push_back(share_idx);

      //sort according to how often a particle is shared
      //----------------------------------------------------------
      // initialize original index locations
      std::vector<size_t> ish(share_idx.size());
      iota(ish.begin(), ish.end(), 0);
      std::sort(ish.begin(), ish.end(),
                [&share_idx](size_t i1, size_t i2) {return share_idx[i1] > share_idx[i2];});
      //create final UE object
      //----------------------------------------------------------
      double maxPtFinalUE = rho_Estimate.at(ijet)-12.92;
      double curPtFinalUE = 0.;
      double prevPtFinalUE = 0.;

      std::vector<fastjet::PseudoJet> bkgd_particles;
      fjJetParticles_.clear();
      for(auto userIndex : ish) {
    //for (int i=0; i<)
    //  if(ijet==0)  cout << userIndex<<endl;
      fastjet::PseudoJet part = particles[userIndex];
       if(curPtFinalUE<maxPtFinalUE) { //assign as bkgd particle
          curPtFinalUE+=part.pt();
  //       if (ijet==0) std::cout << "Final UE: Added new particle with pt = " << part.phi() <<  " iparticle: " << userIndex << " to init condition. total pt now " << curPtFinalUE << "/" << maxPtFinalUE << endl;
          bkgd_particles.push_back(part);
        }
       else {
         fjJetParticles_.push_back(part);}
      }
    //  cout << "ijet: " << ijet << "PtFinalUE: " << curPtFinalUE << endl;
      if(bkgd_particles.size()>0){
      int bkgd_particles_size = bkgd_particles.size();
      fastjet::PseudoJet last_bkg_part = bkgd_particles.at(bkgd_particles_size-1);

      for(int ipart = 0; ipart<bkgd_particles_size; ++ipart){
       prevPtFinalUE+=bkgd_particles.at(ipart).pt();
      }
      // Decide what is a better aproximation

      double distanceUE_one = sqrt((maxPtFinalUE-curPtFinalUE)*(maxPtFinalUE-curPtFinalUE));
      double distanceUE_two = sqrt((maxPtFinalUE-prevPtFinalUE)*(maxPtFinalUE-prevPtFinalUE));

      if (distanceUE_one > distanceUE_two)
      { fjJetParticles_.push_back(last_bkg_part); // add the last particle from the background to the signal
        bkgd_particles.pop_back();} // remove the last particle from the background
      }

      if(fjJetParticles_.size()>0) {
      csJetSub = new fastjet::ClusterSequenceArea(fjJetParticles_, jet_defSub, area_def_sub);
        std::vector<fastjet::PseudoJet> jetSub = fastjet::sorted_by_pt(csJetSub->inclusive_jets());
        if(jetSub[0].pt()>0.) subtracted_jets.push_back(jetSub[0]);
        if(subtracted_jets.size()>0 && subtracted_jets.size()<2)
        csJetSub->delete_self_when_unused();
      }

    }//jet loop
    std::cout << "\n n subtracted jets: " << subtracted_jets.size() << std::endl;
    return subtracted_jets;
  }

  int findClosestParticle(std::vector<fastjet::PseudoJet> particlesNotUsed, int ighost, std::vector<fastjet::PseudoJet> ghosts) {
    //find closest particle to ghost

    double dR2 = 1000.;
    fastjet::PseudoJet partSel;
    int ipSel = -1;
    for(int ip = 0; ip<(int)particlesNotUsed.size(); ++ip) {
      fastjet::PseudoJet part = particlesNotUsed[ip];
      double dR2Cur = ghosts[ighost].squared_distance(part);
      if(dR2Cur < dR2) {
        dR2 = dR2Cur;
        partSel = part;
        ipSel = ip;
      }
    }
    return ipSel;
  }

  std::vector<int> findClosestParticles(std::vector<fastjet::PseudoJet> particles, int ighost, std::vector<fastjet::PseudoJet> ghosts, int nClosest = 5) {
    //find closest particle to ghost

    std::vector<double> dR2(nClosest);
    std::fill(dR2.begin(), dR2.end(), 1000.);

    fastjet::PseudoJet partSel;
    std::vector<int> ipSel(nClosest);
    std::fill(ipSel.begin(), ipSel.end(), -1);

    for(int ip = 0; ip<(int)particles.size(); ++ip) {
      fastjet::PseudoJet part = particles[ip];
      double dR2Cur = ghosts[ighost].squared_distance(part);
      for(int i = 0; i<nClosest; ++i) {
        if(dR2Cur < dR2[i]) {
          dR2[i] = dR2Cur;
          ipSel[i] = ip;
        }
      }
    }

    //sort top nClosest from smallest to largest distance
    std::vector<size_t> idx(dR2.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
              [&dR2](size_t i1, size_t i2) {return dR2[i1] < dR2[i2];});

    std::vector<int> ipSelSorted(nClosest);
    for(int i = 0; i<nClosest; ++i) {
      ipSelSorted[i] = ipSel[idx[i]];
    }

    return ipSelSorted;
  }


  std::vector<double> calculateChi2s(std::vector<std::vector<int> > collInitCond, std::vector<fastjet::PseudoJet> particles, double med_pTD, double rms_pTD) {
    // calc chi2 for each initial condition

    std::vector<double> chi2s;
    for(int ii = 0; ii<collInitCond.size(); ++ii) {
      std::vector<int> indices = collInitCond[ii];
      double chi2 = 1e6;
      if(indices.size()>0) {
        std::vector<fastjet::PseudoJet> combinedparticles;
        for(int ic = 0; ic<(int)indices.size(); ++ic) {
          combinedparticles.push_back(particles[indices[ic]]);
        }
        fastjet::PseudoJet currInitJet = fastjet::PseudoJet(join(combinedparticles));
        double pTDCur = pTD_.result(currInitJet);
        chi2 = fabs(pTDCur-med_pTD)*(fabs(pTDCur-med_pTD))/rms_pTD/rms_pTD;
      }
      chi2s.push_back(chi2);
    }
    return chi2s;
  }

};

#endif

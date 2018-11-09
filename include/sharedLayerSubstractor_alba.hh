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

// Add root libraries for plotting

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "../PU14/PU14.hh"

#include "Angularity.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the jet-by-jet shared layer subtraction
// Author: M. Verweij, Y. Mehtar-Tani
//---------------------------------------------------------------

class sharedLayerSubtractor_alba {

private :
  double jetRParam_;
  double ghostArea_;
  double ghostRapMax_;
  double jetRapMax_;
  double rho_;
  double rhoSigma_;
  double pTDbkg_;
  double pTDbkgSigma_;

  Angularity pTD_;

  std::vector<fastjet::PseudoJet> fjInputs_;
  std::vector<fastjet::PseudoJet> fjJetInputs_;

  std::vector<fastjet::PseudoJet> fjJetParticles_;

  std::random_device rd_; // the seed for the random number generation
  int nInitCond_;
  int nTopInit_;

public :
  sharedLayerSubtractor(double rJet = 0.4,
                        double ghostArea = 0.001,
                        double ghostRapMax = 3.0,
                        double jetRapMax = 3.0,
                        int nInitCond = 200,
                        int nTopInit = 10) :
    jetRParam_(rJet),
    ghostArea_(ghostArea),
    ghostRapMax_(ghostRapMax),
    jetRapMax_(jetRapMax),
    nInitCond_(nInitCond),
    nTopInit_(nTopInit)
  {
    pTD_ = Angularity(0.,2.,0.4);
  }

  void setGhostArea(double a) { ghostArea_ = a; }

  void setInputParticles(std::vector<fastjet::PseudoJet> v) { fjInputs_ = v; }
  void setInputJets(std::vector<fastjet::PseudoJet> v)      { fjJetInputs_ = v; }

  double getRho()  const { return rho_; }
  double getRhoSigma() const { return rhoSigma_; }

  double getPTDBkg() const { return pTDbkg_; }
  double getPTDBkgSigma() const { return pTDbkgSigma_; }

  std::vector<fastjet::PseudoJet> doSubtraction() {

    //if(fjJetInputs_.size()==0 && fjInputs_.size()) {
    //  throw "You didn't give me input jets or particles. You should give me one of the two";
    //  return std::vector<fastjet::PseudoJet>();
    //}

    fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, ghostArea_);

    std::vector<fastjet::PseudoJet> jets = fjJetInputs_; // without the ghosts stuff
    //if(jets.size()==0) {

    // do the clustering with ghosts and get the jets
    //----------------------------------------------------------
    fastjet::JetDefinition jet_def(antikt_algorithm, jetRParam_);
    fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec); //Active areas add a uniform background of extremely soft massless ‘ghost’ particles to the event and allow them to participate in the clustering. The area of a given jet is proportional to the number of ghosts it contains.

    fastjet::ClusterSequenceArea cs(fjInputs_, jet_def, area_def);
    fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax_); // Keep jets with abs(y)<3
    jets = fastjet::sorted_by_pt(jet_selector(cs.inclusive_jets())); // output of the declustering procedure. Return a vector of jets sorted into decreasing transverse momentum that have abs(y)<3.

     fastjet::JetDefinition jet_defSub(antikt_algorithm, 999.); // the radius is huge becase we want all the particles to be accepted at this point as we have already label them as signal
    fastjet::ClusterSequenceArea csSub(fjInputs_, jet_defSub, area_def);

    // create what we need for the background estimation
    //----------------------------------------------------------
    fastjet::JetDefinition jet_def_bkgd(fastjet::kt_algorithm, 0.4);
    fastjet::AreaDefinition area_def_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
    fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4) * (!fastjet::SelectorNHardest(2)); // Select the jets with a different value of maxium rapidity and discard the 2 hardest ones
  //  fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4); // to run with only background

    fastjet::ClusterSequenceArea csKt(fjInputs_, jet_def_bkgd, area_def_bkgd);
    std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

    //Write the rho values for each jet in a tree
  //  TFile *f = new TFile("rho_dist_10events.root", "RECREATE");
  //  double rho_value = 0;
  //  double ptd_value = 0;
  //  TTree *rho_dist = new TTree("rho_dist", "N1");
  //  rho_dist->Branch("rho_value", &rho_value, "rho_value/D");
  //  rho_dist->Branch("ptd_value", &ptd_value, "ptd_value/D");
  //  for(fastjet::PseudoJet& jet : bkgd_jets) {
  //     rho_value= jet.pt()/jet.area();
    //  trw.addCollection("rho", prueba);
  //     rho_dist->Fill();
  //  std::cout << jet.pt()/jet.area() << std::endl;
    //std::cout << jet.area() << std::endl;
//  }

    fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    bkgd_estimator.set_particles(fjInputs_);
    bkgd_estimator.set_jets(bkgd_jets);

    rho_ = bkgd_estimator.rho();
    rhoSigma_ = bkgd_estimator.sigma();

    if(rho_ < 0)    rho_ = 0;

    std::cout << "rho: " << rho_ << "  rhoSigma: " << rhoSigma_ << std::endl;

    //initial gaus with mean=rho_ and width=rhoSigma_
    std::normal_distribution<> gausDist(rho_, rhoSigma_);

    //UE metric
    //----------------------------------------------------------
    std::vector<double> pTD_bkgd;
    for(fastjet::PseudoJet& jet : bkgd_jets) {
    //  ptd_value = pTD_.result(jet);
    //  rho_dist->Fill();
      pTD_bkgd.push_back(pTD_.result(jet));
    }
  //  rho_dist->Write();
    //rho_dist->Scan();
  //  f->Write();
  //  f->Close();
    std::nth_element(pTD_bkgd.begin(), pTD_bkgd.begin() + pTD_bkgd.size()/2, pTD_bkgd.end()); // place the numbers in value order
    double med_pTD = pTD_bkgd[pTD_bkgd.size()/2]; // find the middle

    int nRMS = 0;
    double rms_pTD = 0.;
    for(int ip = 0; ip<(int)pTD_bkgd.size(); ++ip) {
      rms_pTD += (pTD_bkgd[ip]-med_pTD)*(pTD_bkgd[ip]-med_pTD);
      nRMS++;
    }
    if(nRMS>0.)
      rms_pTD = sqrt(rms_pTD/(double)nRMS);

    pTDbkg_ = med_pTD; // median of the background
    pTDbkgSigma_ = rms_pTD; // rms of the background

   //std::cout << "pTDbkg: " << pTDbkg_ << "  pTDbkgSigma: " << pTDbkgSigma_ << std::endl;

    std::mt19937 rndSeed(rd_()); //rnd number generator seed

    std::vector<fastjet::PseudoJet> subtracted_jets;
    subtracted_jets.reserve(jets.size()); //Requests that the vector capacity be at least enough to contain jets.size() elements.
    int ijet = -1;
    for(fastjet::PseudoJet& jet : jets) { // Here we use the full set of jets coming from the initial clustering procedure
      ++ijet; // needed?
      if(jet.is_pure_ghost()) continue; // skip the rest
      //      if(ijet>0) continue;
      //std::cout << "start jet loop. entry: " << ijet << "/" << jets.size() << " pt: " << jet.pt() << " eta: " << jet.eta() << std::endl;

      //get ghosts and true particles (ghosts are distributed uniformly which we will use to create initial conditions)
      //----------------------------------------------------------
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles); // sift the constituents of the input jets into two vectors -- those that pass the selector and those that do not. In this case the selector tells you if it is a ghost or not.
      if(particles.size()<1 || jet.pt()<1.) continue; // ignore jets with only one constituent or with a pt smaller than 1GeV

      //set user_index of all particles to position particles vector
      for(int i = 0; i<(int)particles.size(); ++i) {
        particles[i].set_user_index(i);
      }

      // create map of ghost and closest particle
      std::vector<std::vector<int>> closestPartToGhost;
      for(int ighost = 0; ighost<ghosts.size(); ++ighost) {
        std::vector<int> ipSel = findClosestParticles(particles, ighost, ghosts, 5);
        closestPartToGhost.push_back(ipSel); // indexes of the 5 closest particles to the ghost. For each ghost one has 5 particles.
      }

      // create requested number of initial conditions
      //----------------------------------------------------------
      std::vector<std::vector<int>> collInitCond;
      for(int ii = 0; ii<nInitCond_; ++ii) {
        std::uniform_int_distribution<> distUni(0,ghosts.size()-1); //it's going to be equally probable to pick any of the ghosts
        std::vector<int> initCondition;                           //list of particles in initial condition

        //get random maxPt for this initial condition
        double maxPt = gausDist(rndSeed)*jet.area(); // the rho value is obtained from the gaussian distribution

        //make copy of particles so that a particle is not repeated inside the same initial condition
        std::vector<fastjet::PseudoJet> particlesNotUsed = particles;

        std::vector<std::vector<int>> closestPartToGhostNotUsed = closestPartToGhost;
        double maxPtCurrent = 0.;
        std::vector<int> avail(closestPartToGhost.size()); // New vector as big as the closestPartToGhost one
        std::fill(avail.begin(),avail.end(),1); // Fill it wih ones
        while(maxPtCurrent<maxPt &&  std::accumulate(avail.begin(),avail.end(),0)>0 ) {
          //std::cout << "avail: " << std::accumulate(avail.begin(),avail.end(),0) << std::endl;
          // While we didn't reach the maximum pt and still have ghosts available

          //pick random ghost
          int ighost = int(std::floor(distUni(rndSeed)));
          if(ighost>=ghosts.size()) continue; // skip if there are no more ghosts

          int ipSel = -1;
          if(closestPartToGhostNotUsed[ighost].size()>0) { // if it has particles next to it
            ipSel = closestPartToGhostNotUsed[ighost][0]; // pick the first one
            if(closestPartToGhostNotUsed[ighost].size()<2) avail[ighost] = 0; // If by chance I pick the same ghost it might have that I have already used all the closest particles to it.
            closestPartToGhostNotUsed[ighost].erase(closestPartToGhostNotUsed[ighost].begin()+0); // erase the ighost element from the vector
          } else
            continue;
          
          fastjet::PseudoJet partSel = particles[ipSel]; // new vector that contains the selected particles
          initCondition.push_back(partSel.user_index());
          maxPtCurrent+=partSel.pt();
    //      std::cout << "Added new particle with pt = " << partSel.pt() << " to init condition. total pt now " << maxPtCurrent << "/" << maxPt << std::endl;
        }
        if(maxPtCurrent>maxPt) {
          collInitCond.push_back(initCondition); //avoid putting in a initial condition for which not enough particles were available anymore to get to the required pT. Might be an issue for sparse events.
        }
      }//initial conditions loop.

      //----------------------------------------------------------
      //Now we have the requested number of random initial condition

      //Next step: calc chi2 for each initial condition
      //----------------------------------------------------------
      std::vector<double> chi2s = calculateChi2s(collInitCond, particles, med_pTD, rms_pTD);

      //sort the chi2s keeping track of indices
      //----------------------------------------------------------
      // initialize original index locations
      std::vector<size_t> idx(chi2s.size());
      std::iota(idx.begin(), idx.end(), 0); // Assigns to every element in the range a zero value

      // sort indexes based on comparing values in v
      std::sort(idx.begin(), idx.end(),
                [&chi2s](size_t i1, size_t i2) {return chi2s[i1] < chi2s[i2];}); // sort in decresing order of chi2 values


      //Next step: figure out how often each particle is shared in nTopInit_ initial conditions
      //----------------------------------------------------------
      std::vector<int> share_idx(particles.size(),0);
      for(int it = 0; it<std::min(nTopInit_,(int)collInitCond.size()); ++it) {
        int chi2Index = idx[it];
        std::vector<int> indices = collInitCond[chi2Index];
        for(int ic = 0; ic<(int)indices.size(); ++ic) {
          share_idx[particles[indices[ic]].user_index()]++; // if the same particle appears in more than on initial condition add one
        }
      }

      //sort according to how often a particle is shared
      //----------------------------------------------------------
      // initialize original index locations
      std::vector<size_t> ish(share_idx.size()); // particles size
      iota(ish.begin(), ish.end(), 0);
      std::sort(ish.begin(), ish.end(),
                [&share_idx](size_t i1, size_t i2) {return share_idx[i1] > share_idx[i2];}); // sort in increasing order of times shared values


      //create final UE object
      //----------------------------------------------------------
      double maxPtFinalUE = gausDist(rndSeed)*jet.area(); // not the same as maxPt
      double curPtFinalUE = 0.;
      std::vector<fastjet::PseudoJet> bkgd_particles;
      fjJetParticles_.clear();
      for(auto userIndex : ish) {
        fastjet::PseudoJet part = particles[userIndex];
        if(curPtFinalUE<maxPtFinalUE) { //assign as bkgd particle
          curPtFinalUE+=part.pt();
          bkgd_particles.push_back(part);
        } else { //assign as jet particle
          fjJetParticles_.push_back(part);
        }
      } //Now I have all the particles of the jet classified as bkgd or signal

      if(fjJetParticles_.size()>0) {
        fastjet::ClusterSequenceArea *csSub = new fastjet::ClusterSequenceArea(fjJetParticles_, jet_defSub, area_def);
        std::vector<fastjet::PseudoJet> jetSub = fastjet::sorted_by_pt(csSub->inclusive_jets());
        if(jetSub[0].pt()>0.) subtracted_jets.push_back(jetSub[0]);
        if(subtracted_jets.size()>0 && subtracted_jets.size()<2) csSub->delete_self_when_unused();
      }
    }//jet loop

    std::cout << "\n n subtracted jets: " << subtracted_jets.size() << std::endl;
    //return subtracted_jets;
    return bkgd_jets;
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
    //find 5 closest particles to ghost

    std::vector<double> dR2(nClosest);
    std::fill(dR2.begin(), dR2.end(), 1000.); // Assigns 1000. to the elements in the range [first, last).

    fastjet::PseudoJet partSel;
    std::vector<int> ipSel(nClosest);
    std::fill(ipSel.begin(), ipSel.end(), -1); // Assigns -1 to the elements in the range [first, last).

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
        double ptDCur = pTD_.result(currInitJet);
        chi2 = fabs(ptDCur-med_pTD)*(fabs(ptDCur-med_pTD))/rms_pTD/rms_pTD;
      }
      chi2s.push_back(chi2);
    }
    return chi2s;
  }

};

#endif

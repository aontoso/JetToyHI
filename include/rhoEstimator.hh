#ifndef rhoEstimator_h
#define rhoEstimator_h

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
#include "TF1.h"
#include "TCanvas.h"


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "../PU14/PU14.hh"

#include "Angularity.hh"

using namespace std;
using namespace fastjet;

//---------------------------------------------------------------
// Description
// This class runs the jet-by-jet rho Estimator
// Author: M. Verweij, Y. Mehtar-Tani, Alba Soto-Ontoso
//---------------------------------------------------------------

class rhoEstimator{

private :
  double jetRParam_;
  double ghostArea_;
  double ghostRapMax_;
  double jetRapMax_;
  double rho_;
  double rhoSigma_;

  std::vector<fastjet::PseudoJet> fjInputs_;
  std::vector<fastjet::PseudoJet> fjJetInputs_;

//  std::random_device rd_;

  //fastjet::ClusterSequenceArea *csJetSub;

public :
  rhoEstimator(double rJet = 0.4,
                        double ghostArea = 0.001,
                        double ghostRapMax = 3.0,
                        double jetRapMax = 3.0
                      ):
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

  // By calling this function we can get the value of the "real" pt given the linear assumption
   double correlation_line(double slope, int n_constituents){
   return slope*(n_constituents+1);
   };

  std::vector<double> doEstimation() {

    fastjet::GhostedAreaSpec ghost_spec(ghostRapMax_, 1, ghostArea_);
    fastjet::GhostedAreaSpec ghost_spec_sub(ghostRapMax_, 1, 1.);
    std::vector<fastjet::PseudoJet> jets = fjJetInputs_;

    // Create what we need for the background estimation
    //----------------------------------------------------------
    fastjet::JetDefinition jet_def_bkgd(fastjet::kt_algorithm, jetRParam_);
    fastjet::AreaDefinition area_def_bkgd(fastjet::active_area_explicit_ghosts,ghost_spec);
    fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4)* (!fastjet::SelectorNHardest(2));
  // fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetRapMax_-0.4);

    fastjet::ClusterSequenceArea csKt(fjInputs_, jet_def_bkgd, area_def_bkgd);
    std::vector<fastjet::PseudoJet> bkgd_jets = fastjet::sorted_by_pt(selector(csKt.inclusive_jets()));

    fastjet::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    bkgd_estimator.set_particles(fjInputs_);
    bkgd_estimator.set_jets(bkgd_jets);

    double rhoMedian_ = 0;
    double rhoSigma_ = 0;

    // Compute the rho median for each event and its standard deviation

    rhoMedian_ = bkgd_estimator.rho();
    rhoSigma_ = bkgd_estimator.sigma();

    // Determine the ptCut a la Soft Killer
    //-------------------------------------------------------------

    std::vector<double> ptConst_bkgd;

    for(fastjet::PseudoJet& jet : bkgd_jets) {
        double maxpt_bkgd = 0;
        std::vector<fastjet::PseudoJet> particles, ghosts; // make sure that the ghosts do not screw up the pt spectrum
        fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
     // Find the maximum-pt of patch in the background
       for(fastjet::PseudoJet p : particles) {
            double momentum = p.pt();
            if (momentum > maxpt_bkgd) maxpt_bkgd = momentum;
          }
         ptConst_bkgd.push_back(maxpt_bkgd);
       } // bkgd_jets loop

       std::nth_element(ptConst_bkgd.begin(), ptConst_bkgd.begin() + ptConst_bkgd.size()/2, ptConst_bkgd.end());
       double pTsoftKiller = ptConst_bkgd[ptConst_bkgd.size()/2];

    //   cout << pTsoftKiller << endl;

    // Determination of the slope ("temperature") of the correlation line
    //----------------------------------------------------------

    double temperature = 1.2;

    // Do the anti-kt clustering
    //------------------------------------------------------------

    fastjet::JetDefinition jet_def(antikt_algorithm, jetRParam_);
    fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,ghost_spec);

    fastjet::ClusterSequenceArea cs(fjInputs_, jet_def, area_def);
    fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax_);
    jets = fastjet::sorted_by_pt(jet_selector(cs.inclusive_jets()));

    std::vector<double> rho_estimate; // one for each patch;
    rho_estimate.reserve(jets.size());

    int ijet = -1;

    for(fastjet::PseudoJet& jet : jets) {
      ++ijet;
      if(jet.is_pure_ghost()) continue;

      //Get ghosts and true particles
      //----------------------------------------------------------
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);
      if(particles.size()<1 || jet.pt()<1.) continue;

      // Raw values
      double truePatchPt = 0;
      int truePatchN = 0;
      for(fastjet::PseudoJet p : particles) {
      // Obtain the true rho for each patch by using only the background particles
       truePatchPt+=p.pt();
       truePatchN++;
       }
    //   if (ijet==0)

  //  if (truePatchPt-jet.pt()>5) cout << "Jet.pt(): " << jet.pt() << " Patch Pt: " << truePatchPt-jet.pt() << endl;

      // Information to plot the trajectories

        int tamano = particles.size();
	      Double_t pt_progress[tamano];
	      Double_t n_steps[tamano];

        for (int j=0; j<tamano; j++){
	          pt_progress[j]=0;
	          n_steps[j] = 0;
	        }
      // Order particles from soft to high pT

      std::vector<size_t> idx(particles.size());
      iota(idx.begin(), idx.end(), 0);

      std::sort(idx.begin(), idx.end(),
                [&particles](size_t i1, size_t i2) {return particles[i1].pt() < particles[i2].pt();});
      std::vector<fastjet::PseudoJet> particles_pTOrdered;
      for (int i = 0; i<particles.size(); i++){
         particles_pTOrdered.push_back(particles[idx[i]]);

      }

      // Set user_index of all particles to position particles vector

       for(int i = 0; i<(int)particles_pTOrdered.size(); ++i) {
         particles_pTOrdered[i].set_user_index(i);
       }

      std::vector<fastjet::PseudoJet> bkgEstimate;  //list of particles until we reach the correlation line
      std::vector<int> avail_part(particles_pTOrdered.size()); // list of particles still available
      double pTcurrent = 0; // to make sure it goes into the while loop
      int nparticles = -1;

      double pTcorrelation = 1;

      std::fill(avail_part.begin(),avail_part.end(),1);

    //  if (correlation_line(temperature, truePatchN) > truePatchPt){
    //   pTcurrent = truePatchPt;

//      }
      double pTcut = pTsoftKiller; // new trial
    //  else {
      while(pTcurrent<rhoMedian_*jet.area() && std::accumulate(avail_part.begin(),avail_part.end(),0)>0 && particles_pTOrdered[nparticles+1].pt()<pTcut){

  //   while(std::accumulate(avail_part.begin(),avail_part.end(),0)>0 && particles_pTOrdered[nparticles+1].pt()<pTcut){

      nparticles++;
      pTcorrelation = correlation_line(temperature, nparticles);
      fastjet::PseudoJet partSel = particles_pTOrdered[nparticles];
      bkgEstimate.push_back(partSel);
      pTcurrent+= partSel.pt();
      avail_part[nparticles] = 0;
      pt_progress[nparticles] = pTcurrent;
      n_steps[nparticles] = nparticles+1;

    //  if (ijet == 9) cout << "Signal/Bkg particle: " << partSel.user_info<PU14>().vertex_number() << "nparticle: " << nparticles+1 << "pTcurrent: " << pTcurrent << endl;
      }
    //  pTcurrent = correlation_line(temperature,nparticles);

    //}

    //Compute the shortest distance to the correlation line
    //----------------------------------------------------------

    //  double minus_slope_perp = 1/temperature;
    //  double intercept_perp = pTcurrent + (nparticles+1)*minus_slope_perp;
    //  int n_shortest = int(intercept_perp/(temperature+minus_slope_perp));
    //  double pTshortest = correlation_line(temperature, n_shortest);
    //  pTcurrent = pTshortest;
  //  }


      //Compute delta_rho as given by Eq. 4.36 in Yacine's notes
      //----------------------------------------------------------
  //   if(pTcurrent>pTcorrelation && std::accumulate(avail_part.begin(),avail_part.end(),0)>0)
    //  {
    //  double rho_star = pTcurrent;
    //  double n_star = nparticles+1;
    //  double rho_raw = truePatchPt;
    //  double n_raw = truePatchN;
    //  double deltaRho = temperature*pow(n_raw-n_star,2)*rho_star/((rho_raw-rho_star)*n_star);

      //if(ijet==0) { cout << truePatchPt << endl;
      //  cout << truePatchN << endl;
      //  cout <<deltaRho <<endl;
      //}
    //  cout << deltaRho << endl;
  //    pTcurrent-=deltaRho;
      // }

      //Compute delta_rho as given by Eq. 2.6 in Yacine's notes
      //----------------------------------------------------------
  //    if(pTcurrent>pTcorrelation && std::accumulate(avail_part.begin(),avail_part.end(),0)>0)
    //     {
  //    double rho_star = particles_pTOrdered[nparticles].pt();
  //    double n_star = nparticles+1;
    //  double n_raw = truePatchN;
  //    double delta = n_star*(rho_star+temperature)*exp(-rho_star/temperature);
      //if(ijet==1) cout << n_star << endl;
  //      pTcurrent-=delta;
  //  }

  //  if()

  //  if(pTcurrent<pTcorrelation && std::accumulate(avail_part.begin(),avail_part.end(),0)>0)
  //     {
  //  double rho_star = particles_pTOrdered[nparticles].pt();
  //  double n_star = nparticles+1;
  //  double n_raw = truePatchN;
  //  double delta = n_star*exp(-rho_star/temperature);
    //if(ijet==1) cout << n_star << endl;
  //    pTcurrent+=delta;
  //}

  //Compute delta_rho as given by Eq. 2.14 in Yacine's notes
  //----------------------------------------------------------
     if(pTcurrent>pTcorrelation && std::accumulate(avail_part.begin(),avail_part.end(),0)>0)
     {

      double rho_star = particles_pTOrdered[nparticles].pt();
      double n_star = nparticles+1;
      double deltaRho = n_star*exp(-2*rho_star/temperature)*(pow(temperature,2)+2*rho_star*temperature+2*pow(rho_star,2)*temperature);
  //    cout << deltaRho << endl;
      pTcurrent-=deltaRho;

     }
  //
    if (std::accumulate(avail_part.begin(),avail_part.end(),0)>0){
    double minus_slope_perp = 1/temperature;
    double intercept_perp = pTcurrent + (nparticles+1)*minus_slope_perp;
    int n_shortest = int(intercept_perp/(temperature+minus_slope_perp));
    double pTshortest = correlation_line(temperature, n_shortest);
    pTcurrent = pTshortest;
    }

  //   if (pTcurrent > truePatchPt) pTcurrent = truePatchPt;
 //  }

      //Check the effect of the last particle in a fancy way
      //----------------------------------------------------------

  //    if (pTcurrent>pTcorrelation && std::accumulate(avail_part.begin(),avail_part.end(),0)>3){

        // Find the intersection between the two lines
    //     double x1 = nparticles-1;
  //       double x2= nparticles;
  //       double y1 = pTcurrent;
  //       double y2 = pTcurrent - particles_pTOrdered[nparticles].pt();
  //       double slope_two = (y2-y1)/(x2-x1);
  //       double intercept_two = y2 - slope_two*x2;
  //       double x_intersection = intercept_two/(temperature-slope_two);
  //       pTcurrent=correlation_line(temperature, int(x_intersection));
      //   double distance_two = ;

    //  }
      rho_estimate.push_back(pTcurrent);
    //  cout << "Estimate: " << pTcurrent << "Rho median: " << rhoMedian_*jet.area()<< endl;

  //  if (ijet == 9) {
    //     int nsteps = nparticles+1;
    //     TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 1200, 800);
    //     TGraph *pt_evolution = new TGraph(nsteps,n_steps, pt_progress);
    //     pt_evolution->Draw("AP");
    //     TF1 *facorrelation_1 = new TF1("correlation_","1.2*x",0,nsteps);
    //     facorrelation_1->Draw("SAME");
    //     c1->SaveAs("Prueba_rho9.C");
    //   }
    } // jets_loop

    return rho_estimate;


};

};

#endif

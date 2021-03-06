#ifndef softDropCounter_h
#define softDropCounter_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Recluster.hh"

//---------------------------------------------------------------
// Description
// This class runs iterative SoftDrop on a set of jets
// Author: Y. Chen
//---------------------------------------------------------------

class softDropCounter
{
private :
   double zcut_;   // ZCut used in SD
   double beta_;   // Beta parameter
   double r0_;     // Jet radius
   double rcut_;   // Termination Criteria

   std::vector<fastjet::PseudoJet> fjInputs_;     // ungroomed jets
   std::vector<std::vector<double>> zgs_;         // all the zg's in the algorithm
   std::vector<std::vector<double>> drs_;         // and the angles in the algorithm
   std::vector<std::vector<double>> pts_;         // and the angles in the algorithm
   std::vector<std::vector<double>> masses_;         // and the mass in the algorithm
   std::vector<std::vector<fastjet::PseudoJet>> nodes_4momenta_; // returns the sum of the daughters 4 momenta at every node
public :
   softDropCounter(double z = 0.1, double beta = 0.0, double r0 = 0.4, double rcut = 0.1);
   void setZCut(double c);
   void setBeta(double b);
   void setR0(double r) ;
   void setRCut(double r);
   void setInputJets(const std::vector<fastjet::PseudoJet> &v);
   void run(const jetCollection &c);
   void run(const std::vector<fastjet::PseudoJet> &v);
   void run();
   std::vector<double> calculateNSD(double Kappa, double AngleKappa = 0);
   std::vector<std::vector<double>> calculateTf();
   std::vector<std::vector<double>> calculateDeltaR();
   std::vector<std::vector<double>> calculateZg();
   std::vector<std::vector<double>> calculatekT();
   std::vector<std::vector<double>> calculateMass();
   std::vector<std::vector<double>> calculatepT();
   std::vector<std::vector<fastjet::PseudoJet>> getDipole4momenta();
   std::vector<std::vector<double>> getDipolePx();
   std::vector<std::vector<double>> getDipolePy();
   std::vector<std::vector<double>> getDipolePz();
   std::vector<std::vector<double>> getDipoleE();
};

softDropCounter::softDropCounter(double z, double beta, double r0, double rcut)
   : zcut_(z), beta_(beta), r0_(r0), rcut_(rcut)
{
}

void softDropCounter::setZCut(double c)
{
   zcut_ = c;
}

void softDropCounter::setBeta(double b)
{
   beta_ = b;
}

void softDropCounter::setR0(double r)
{
   r0_   = r;
}

void softDropCounter::setRCut(double r)
{
   rcut_ = r;
}

void softDropCounter::setInputJets(const std::vector<fastjet::PseudoJet> &v)
{
   fjInputs_ = v;
}

void softDropCounter::run(const jetCollection &c)
{
   run(c.getJet());
}

void softDropCounter::run(const std::vector<fastjet::PseudoJet> &v)
{
   setInputJets(v);
   run();
}

void softDropCounter::run()
{
   //int N = fjInputs_.size();

   for(fastjet::PseudoJet &jet: fjInputs_)
   {
      if(jet.has_constituents() == false)
      {
         zgs_.push_back(vector<double>());
         drs_.push_back(vector<double>());
         pts_.push_back(vector<double>());
         masses_.push_back(vector<double>());
         nodes_4momenta_.push_back(std::vector<fastjet::PseudoJet>());
         continue;
      }

      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm,1.);
      fastjet::ClusterSequence cs(particles, jet_def);
      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());

      if(tempJets.size() == 0)
      {
         zgs_.push_back(vector<double>());
         drs_.push_back(vector<double>());
         pts_.push_back(vector<double>());
         masses_.push_back(vector<double>());
         nodes_4momenta_.push_back(std::vector<fastjet::PseudoJet>());
         continue;
      }

      PseudoJet CurrentJet = tempJets[0];
      PseudoJet Part1, Part2;

      std::vector<double> z;
      std::vector<double> dr;
      std::vector<double> pt;
      std::vector<double> mass;
      std::vector<fastjet::PseudoJet> node_4momenta;

      while(CurrentJet.has_parents(Part1, Part2))
      {
         if(CurrentJet.pt2() <= 0)
            break;

         double DeltaR = std::sqrt(Part1.squared_distance(Part2));
         if(DeltaR < rcut_)
            break;

         double PT1 = Part1.pt();
         double PT2 = Part2.pt();
         double zg = -1;

         if(PT1 + PT2 > 0)
            zg = min(PT1, PT2) / (PT1 + PT2);
         else
            break;

         double Threshold = zcut_ * std::pow(DeltaR / r0_, beta_);

         if(zg >= Threshold)   // yay
         {  fastjet::PseudoJet transformedJet;
            z.push_back(zg);
            dr.push_back(DeltaR);
            pt.push_back(PT1+PT2);
            transformedJet = fastjet::PseudoJet(Part1.px()+Part2.px(),Part1.py()+Part2.py(),Part1.pz()+Part2.pz(),Part1.E()+Part2.E());
            node_4momenta.push_back(transformedJet);
            mass.push_back(transformedJet.m());
         }

         if(PT1 > PT2)
            CurrentJet = Part1;
         else
            CurrentJet = Part2;
      }

      zgs_.push_back(z);
      drs_.push_back(dr);
      pts_.push_back(pt);
      masses_.push_back(mass);
      nodes_4momenta_.push_back(node_4momenta);

   }
}

std::vector<std::vector<fastjet::PseudoJet>> softDropCounter::getDipole4momenta()
{

   return nodes_4momenta_;
}

std::vector<std::vector<double>> softDropCounter::getDipolePx()
{

  std::vector<std::vector<double>> Result_three;

  for(int i = 0; i < (int)zgs_.size(); i++)
  {
     std::vector<double> onejet_px;
     onejet_px.reserve(zgs_[i].size());
     for(int j = 0; j < (int)zgs_[i].size(); j++){
        double px = nodes_4momenta_[i][j].px();
        onejet_px.push_back(px);
      }
     Result_three.push_back(onejet_px);
  }

  return Result_three;
}

std::vector<std::vector<double>> softDropCounter::getDipolePy()
{

  std::vector<std::vector<double>> Result_three;

  for(int i = 0; i < (int)zgs_.size(); i++)
  {
     std::vector<double> onejet_py;
     onejet_py.reserve(zgs_[i].size());
     for(int j = 0; j < (int)zgs_[i].size(); j++){
        double py = nodes_4momenta_[i][j].py();
        onejet_py.push_back(py);
      }
     Result_three.push_back(onejet_py);
  }

  return Result_three;
}

std::vector<std::vector<double>> softDropCounter::getDipolePz()
{

  std::vector<std::vector<double>> Result_three;

  for(int i = 0; i < (int)zgs_.size(); i++)
  {
     std::vector<double> onejet_pz;
     onejet_pz.reserve(zgs_[i].size());
     for(int j = 0; j < (int)zgs_[i].size(); j++){
        double pz = nodes_4momenta_[i][j].pz();
        onejet_pz.push_back(pz);
      }
     Result_three.push_back(onejet_pz);
  }

  return Result_three;
}

std::vector<std::vector<double>> softDropCounter::getDipoleE()
{

  std::vector<std::vector<double>> Result_three;

  for(int i = 0; i < (int)zgs_.size(); i++)
  {
     std::vector<double> onejet_energy;
     onejet_energy.reserve(zgs_[i].size());
     for(int j = 0; j < (int)zgs_[i].size(); j++){
        double energy = nodes_4momenta_[i][j].E();
        onejet_energy.push_back(energy);
      }
     Result_three.push_back(onejet_energy);
  }

  return Result_three;
}

std::vector<double> softDropCounter::calculateNSD(double Kappa, double AngleKappa)
{
   std::vector<double> Result;

   for(int i = 0; i < (int)zgs_.size(); i++)
   {
      double Total = 0;
      for(int j = 0; j < (int)zgs_[i].size(); j++)
         Total = Total + pow(zgs_[i][j], Kappa) * pow(drs_[i][j], AngleKappa);
      Result.push_back(Total);
   }

   return Result;
}

std::vector<std::vector<double>> softDropCounter::calculateTf()
{
   std::vector<std::vector<double>> Result_two;

   for(int i = 0; i < (int)zgs_.size(); i++)
   {
      std::vector<double> onejet_tf;
      onejet_tf.reserve(zgs_[i].size());
      for(int j = 0; j < (int)zgs_[i].size(); j++){
         double mass2 = (zgs_[i][j]*(1-zgs_[i][j])*pow(pts_[i][j],2)*pow(drs_[i][j],2))/(pow(r0_,2));
         double GeVtofm = 5.068;
         double formation_time = (2*pts_[i][j])/(mass2*GeVtofm);
         onejet_tf.push_back(formation_time);
       }
      Result_two.push_back(onejet_tf);
   }

   return Result_two;
}

std::vector<std::vector<double>> softDropCounter::calculateDeltaR()
{
   std::vector<std::vector<double>> Result_three;

   for(int i = 0; i < (int)zgs_.size(); i++)
   {
      std::vector<double> onejet_deltaR;
      onejet_deltaR.reserve(zgs_[i].size());
      for(int j = 0; j < (int)zgs_[i].size(); j++){
         double delta_r = drs_[i][j];
         onejet_deltaR.push_back(delta_r);
       }
      Result_three.push_back(onejet_deltaR);
   }

   return Result_three;
}

std::vector<std::vector<double>> softDropCounter::calculateZg()
{
   std::vector<std::vector<double>> Result_four;

   for(int i = 0; i < (int)zgs_.size(); i++)
   {
      std::vector<double> onejet_zg;
      onejet_zg.reserve(zgs_[i].size());
      for(int j = 0; j < (int)zgs_[i].size(); j++){
         double zg = zgs_[i][j];
         onejet_zg.push_back(zg);
       }
      Result_four.push_back(onejet_zg);
   }

   return Result_four;
}

std::vector<std::vector<double>> softDropCounter::calculatekT()
{
   std::vector<std::vector<double>> Result_five;

   for(int i = 0; i < (int)zgs_.size(); i++)
   {
      std::vector<double> onejet_kts;
      onejet_kts.reserve(zgs_[i].size());
      for(int j = 0; j < (int)zgs_[i].size(); j++){
         double kt = zgs_[i][j]*(1-zgs_[i][j])*pts_[i][j]*drs_[i][j]/r0_;
         onejet_kts.push_back(kt);
       }
      Result_five.push_back(onejet_kts);
   }

   return Result_five;
}

std::vector<std::vector<double>> softDropCounter::calculateMass()
{
   std::vector<std::vector<double>> Result_six;

   for(int i = 0; i < (int)masses_.size(); i++)
   {
      std::vector<double> onejet_mass;
      onejet_mass.reserve(masses_[i].size());
      for(int j = 0; j < (int)masses_[i].size(); j++){
         double mass = masses_[i][j];
         onejet_mass.push_back(mass);
       }
      Result_six.push_back(onejet_mass);
   }

   return Result_six;
}

std::vector<std::vector<double>> softDropCounter::calculatepT()
{
   std::vector<std::vector<double>> Result_seven;

   for(int i = 0; i < (int)pts_.size(); i++)
   {
      std::vector<double> onejet_pt;
      onejet_pt.reserve(pts_[i].size());
      for(int j = 0; j < (int)pts_[i].size(); j++){
         double pt = pts_[i][j];
         onejet_pt.push_back(pt);
       }
      Result_seven.push_back(onejet_pt);
   }

   return Result_seven;
}


#endif

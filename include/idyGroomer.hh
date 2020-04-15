#ifndef idyGroomer_h
#define idyGroomer_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/WrappedStructure.hh"
#include "fastjet/contrib/RecursiveSymmetryCutBase.hh"

#include "jetCollection.hh"

//---------------------------------------------------------------
// Description
// This class runs iterative Dynamical grooming on a set of jets
// Author: Alba Soto-Ontoso
//---------------------------------------------------------------

class idyGroomer {

private :
  double a_;
  int niter_; // number of times that we want to iterate (-1 == infinite)
  std::vector<fastjet::PseudoJet> fjInputs_;   //ungroomed jets
  std::vector<fastjet::PseudoJet> fjOutputs_;  //groomed jets
  std::vector<std::vector<double>>             zg_;         //zg of groomed jets
  std::vector<std::vector<int>>              drBranches_; //dropped branches
  std::vector<std::vector<double>>            dr12_;       //distance between the two subjet
  std::vector<double>             tau32_; //n-subjettiness ratio 32
public :
  idyGroomer(double a = 1, int niter = -1);
  void setInputJets(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> getGroomedJets() const;
  std::vector<std::vector<double>> getZgs() const;
  std::vector<std::vector<double>> getDR12() const;
  std::vector<std::vector<int>> getNDroppedSubjets() const;
  std::vector<double> getTau32() const;
  fastjet::PseudoJet getSoftBranch(fastjet::PseudoJet piece1, fastjet::PseudoJet piece2);
  double getKappa(double pt, double theta, double z);
  std::vector<fastjet::PseudoJet> doGrooming(jetCollection &c);
  std::vector<fastjet::PseudoJet> doGrooming(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> doGrooming();
};

idyGroomer::idyGroomer(double a, int niter)
   : a_(a), niter_(niter)
{
}

void idyGroomer::setInputJets(const std::vector<fastjet::PseudoJet> v)
{
   fjInputs_ = v;
}

std::vector<double> idyGroomer::getTau32() const
{
   return tau32_;
}

std::vector<fastjet::PseudoJet> idyGroomer::getGroomedJets() const
{
   return fjOutputs_;
}

fastjet::PseudoJet idyGroomer::getSoftBranch(fastjet::PseudoJet piece1, fastjet::PseudoJet piece2)
{
  fastjet::PseudoJet softBranch;
  if(piece1.pt() < piece2.pt())
  softBranch = piece1;
  else
  softBranch = piece2;
  return softBranch;
}

std::vector<std::vector<double>> idyGroomer::getZgs() const
{
   return zg_;
}

std::vector<std::vector<double>> idyGroomer::getDR12() const
{
   return dr12_;
}

std::vector<std::vector<int>> idyGroomer::getNDroppedSubjets() const
{
   return drBranches_;
}

std::vector<fastjet::PseudoJet> idyGroomer::doGrooming(jetCollection &c)
{
   return doGrooming(c.getJet());
}

std::vector<fastjet::PseudoJet> idyGroomer::doGrooming(std::vector<fastjet::PseudoJet> v)
{
   setInputJets(v);
   return doGrooming();
}

double idyGroomer::getKappa(double pt, double theta, double z)
{
   double kappa = 1/(z * (1-z) * pt * pow(theta,a_));
   return kappa;
}

std::vector<fastjet::PseudoJet> idyGroomer::doGrooming()
{
   fjOutputs_.reserve(fjInputs_.size());
   tau32_.reserve(fjInputs_.size());

   int ijet = -1;

   for(fastjet::PseudoJet& jet : fjInputs_) {
      ++ijet;
      std::vector<fastjet::PseudoJet> particles, ghosts;
      fastjet::SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

      fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 1.);
      fastjet::ClusterSequence cs(particles, jet_def);
      vector<fastjet::ClusterSequence::history_element> cs_history = cs.history();
      std::vector<fastjet::PseudoJet> tempJets = fastjet::sorted_by_pt(cs.inclusive_jets());

      // To make the jetp_index work we need to define our jets like this
      std::vector<fastjet::PseudoJet> tempJets_two = cs.jets();

      if(tempJets.size()<1) {
         fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
         zg_.push_back(vector<double>());
         dr12_.push_back(vector<double>());
         drBranches_.push_back(vector<int>());
         break;
      }

      fastjet::PseudoJet CurrentJet = tempJets[0];
      fastjet::PseudoJet piece1, piece2;

      //int current_in_ca_tree = -1; // (history) index of the current particle in the C/A tree

    //  tuple <double, double, int> tuple_kappa_deltar_node;
      std::vector<tuple<double, double, int>> vector_tuple;

      // now recurse into the primary branch and find the hardness at each node. Saving it together with the angular separation among the subjets and the clustering history index into a tuple.

  while (CurrentJet.has_parents(piece1, piece2)) {

     if (CurrentJet.pt2() <= 0) break;

     if(piece1.pt() + piece2.pt() > 0 && piece1.E()>0. && piece2.E()>0. && piece1.m()>0. && piece2.m()>0.){

    double pt = piece1.pt() + piece2.pt();
    double zg = min(piece1.pt(), piece2.pt()) / pt;
    double deltaR = piece1.delta_R(piece2);
    double kappa = getKappa(pt,deltaR,zg);
    int current_in_ca_tree = CurrentJet.cluster_hist_index();
    tuple <double, double, int> tuple_kappa_deltar_node = make_tuple(kappa, deltaR, current_in_ca_tree);
    vector_tuple.push_back(tuple_kappa_deltar_node);
      }
    if(piece1.pt() > piece2.pt())
      CurrentJet = piece1;
    else
      CurrentJet = piece2;
  }
  // Sort the  list of tuples in the dropped variable keeping the original position in the CA declustering fixed
    sort(vector_tuple.begin(), vector_tuple.end());

  // This list contains the nodes of the C/A history that survive after grooming

  //  std::vector<int> accepted_branches;
  //  accepted_branches.reserve(vector_tuple.size());

  // Keep the branches that are angular ordered
    if(vector_tuple.size()>0){

   double deltaR_hardest = get<1>(vector_tuple[0]);
   int ca_tree_node_hardest = get<2>(vector_tuple[0]);
   fastjet::PseudoJet primary,secondary;

   primary = tempJets_two[cs_history[ca_tree_node_hardest].jetp_index];
   secondary = fastjet::PseudoJet(0.,0.,0.,0.);

   if(vector_tuple.size()>1){

   double deltaR_nxthardest = get<1>(vector_tuple[1]);
   int ca_tree_node_nxthardest = get<2>(vector_tuple[1]);
   primary = fastjet::PseudoJet(0.,0.,0.,0.);
   secondary = fastjet::PseudoJet(0.,0.,0.,0.);

   if(deltaR_nxthardest > deltaR_hardest){

   primary = tempJets_two[cs_history[ca_tree_node_hardest].jetp_index];
   fastjet::PseudoJet child1, child2,dipole;
   dipole = tempJets_two[cs_history[ca_tree_node_nxthardest].jetp_index];
   dipole.has_parents(child1, child2);
   secondary = getSoftBranch(child1, child2);

 }
   else{

   primary = tempJets_two[cs_history[ca_tree_node_nxthardest].jetp_index];
   fastjet::PseudoJet child1, child2, dipole;
   dipole = tempJets_two[cs_history[ca_tree_node_hardest].jetp_index];
   dipole.has_parents(child1, child2);
   secondary = getSoftBranch(child1, child2);

 }
 }
   fastjet::PseudoJet groomed_jet = join(primary, secondary);
   fjOutputs_.push_back(groomed_jet);

   // Compute the n-subjettiness ratio
   double beta = 2;

   fastjet::contrib::NsubjettinessRatio nSub32_beta2(3,2, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));

   double tau32_beta2 = nSub32_beta2(groomed_jet);

   tau32_.push_back(tau32_beta2);

} // if vector_tuple.size() ==0
   else fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
 }
   return fjOutputs_;
}


#endif

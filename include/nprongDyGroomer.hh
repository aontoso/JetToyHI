#ifndef nprongDyGroomer_h
#define nprongDyGroomer_h

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
// This class runs 3-prong Dynamical grooming on a set of jets
// Author: Alba Soto-Ontoso
//---------------------------------------------------------------

class nprongDyGroomer {

private :
  double a_;
  int niter_; // number of times that we want to iterate (-1 == infinite)
  std::vector<fastjet::PseudoJet> fjInputs_;   //ungroomed jets
  std::vector<fastjet::PseudoJet> fjOutputs_;  //groomed jets
  std::vector<std::vector<double>>             zg_;         //zg of groomed jets
  std::vector<std::vector<int>>              drBranches_; //dropped branches
  std::vector<std::vector<double>>            dr12_;       //distance between the two subjet
  std::vector<double>             tau32_; //n-subjettiness ratio 32
  std::vector<double> tau21_; //n-subjettiness ratio 21
public :
  nprongDyGroomer(double a = 1, int niter = -1);
  void setInputJets(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> getGroomedJets() const;
  std::vector<std::vector<double>> getZgs() const;
  std::vector<std::vector<double>> getDR12() const;
  std::vector<std::vector<int>> getNDroppedSubjets() const;
  std::vector<double> getTau32() const;
  std::vector<double> getTau21() const;
  fastjet::PseudoJet getSoftBranch(fastjet::PseudoJet piece1, fastjet::PseudoJet piece2);
  fastjet::PseudoJet getHardBranch(fastjet::PseudoJet piece1, fastjet::PseudoJet piece2);
  double getKappa(double pt, double theta, double z);
  std::vector<fastjet::PseudoJet> doGrooming(jetCollection &c);
  std::vector<fastjet::PseudoJet> doGrooming(std::vector<fastjet::PseudoJet> v);
  std::vector<fastjet::PseudoJet> doGrooming();
};

nprongDyGroomer::nprongDyGroomer(double a, int niter)
   : a_(a), niter_(niter)
{
}

void nprongDyGroomer::setInputJets(const std::vector<fastjet::PseudoJet> v)
{
   fjInputs_ = v;
}

std::vector<double> nprongDyGroomer::getTau21() const
{
   return tau21_;
}

std::vector<double> nprongDyGroomer::getTau32() const
{
   return tau32_;
}

std::vector<fastjet::PseudoJet> nprongDyGroomer::getGroomedJets() const
{
   return fjOutputs_;
}

fastjet::PseudoJet nprongDyGroomer::getSoftBranch(fastjet::PseudoJet piece1, fastjet::PseudoJet piece2)
{
  fastjet::PseudoJet softBranch;
  if(piece1.pt() < piece2.pt())
  softBranch = piece1;
  else
  softBranch = piece2;
  return softBranch;
}

fastjet::PseudoJet nprongDyGroomer::getHardBranch(fastjet::PseudoJet piece1, fastjet::PseudoJet piece2)
{
  fastjet::PseudoJet hardBranch;
  if(piece1.pt() > piece2.pt())
  hardBranch = piece1;
  else
  hardBranch = piece2;
  return hardBranch;
}

std::vector<std::vector<double>> nprongDyGroomer::getZgs() const
{
   return zg_;
}

std::vector<std::vector<double>> nprongDyGroomer::getDR12() const
{
   return dr12_;
}

std::vector<std::vector<int>> nprongDyGroomer::getNDroppedSubjets() const
{
   return drBranches_;
}

std::vector<fastjet::PseudoJet> nprongDyGroomer::doGrooming(jetCollection &c)
{
   return doGrooming(c.getJet());
}

std::vector<fastjet::PseudoJet> nprongDyGroomer::doGrooming(std::vector<fastjet::PseudoJet> v)
{
   setInputJets(v);
   return doGrooming();
}

double nprongDyGroomer::getKappa(double pt, double theta, double z)
{
   double kappa = 1/(z * (1-z) * pt * pow(theta,a_));
   return kappa;
}

std::vector<fastjet::PseudoJet> nprongDyGroomer::doGrooming()
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


    std::vector<tuple<double, double, int>> vector_tuple_lp;

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
    vector_tuple_lp.push_back(tuple_kappa_deltar_node);
      }
    if(piece1.pt() > piece2.pt())
      CurrentJet = piece1;
    else
      CurrentJet = piece2;
  }
  // Sort the  list of tuples in the dropped variable keeping the original position in the CA declustering fixed
    sort(vector_tuple_lp.begin(), vector_tuple_lp.end());


  // Keep the branches that are angular ordered
    if(vector_tuple_lp.size()>0){

   int ca_tree_node_hardest = get<2>(vector_tuple_lp[0]);
   fastjet::PseudoJet primary,secondary;

   fastjet::PseudoJet leading_dipole, leading_child1, leading_child2, softest_child;
   leading_dipole = tempJets_two[cs_history[ca_tree_node_hardest].jetp_index];;
   leading_dipole.has_parents(leading_child1, leading_child2);
   softest_child  = getSoftBranch(leading_child1, leading_child2);

   //////////////////////////////////////////////////////////////////////
   // Navigate through the secondary Lund Plane of the softest daughter
   /////////////////////////////////////////////////////////////////////

   fastjet::PseudoJet CurrentJet_slp = softest_child;
   fastjet::PseudoJet piece1_slp, piece2_slp;
   std::vector<tuple<double, double, int>> vector_tuple_slp;

   while (CurrentJet_slp.has_parents(piece1_slp, piece2_slp)) {

    if (CurrentJet_slp.pt2() <= 0) break;

    if(piece1_slp.pt() + piece2_slp.pt() > 0 && piece1_slp.E()>0. && piece2_slp.E()>0. && piece1_slp.m()>0. && piece2_slp.m()>0.){

     double pt = piece1_slp.pt() + piece2_slp.pt();
     double zg = min(piece1_slp.pt(), piece2_slp.pt()) / pt;
     double deltaR = piece1_slp.delta_R(piece2_slp);

     double kappa = getKappa(pt,deltaR,zg);
     int current_in_ca_tree = CurrentJet_slp.cluster_hist_index();
     tuple <double, double, int> tuple_kappa_deltar_node = make_tuple(kappa, deltaR, current_in_ca_tree);
     vector_tuple_slp.push_back(tuple_kappa_deltar_node);
       }
     if(piece1_slp.pt() > piece2_slp.pt())
       CurrentJet_slp = piece1_slp;
     else
       CurrentJet_slp = piece2_slp;
   }
    //New list with the secondary lp nodes

    sort(vector_tuple_slp.begin(), vector_tuple_slp.end());

   // SLP empty, nxt hardest splitting on the LP
   if(vector_tuple_slp.size()<1 && vector_tuple_lp.size()>1){

   int ca_tree_node_nxthardest = get<2>(vector_tuple_lp[1]);


   if(ca_tree_node_nxthardest > ca_tree_node_hardest){

   primary = tempJets_two[cs_history[ca_tree_node_hardest].jetp_index];
   fastjet::PseudoJet child1, child2, dipole;
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

  // nxt hardest splitting  on the SLP, no nxt hardest splitting on the LP

  else if(vector_tuple_slp.size()>0 && vector_tuple_lp.size()<1){

    int ca_tree_node_nxthardest = get<2>(vector_tuple_slp[0]);
    primary = tempJets_two[cs_history[ca_tree_node_nxthardest].jetp_index];
    fastjet::PseudoJet child1, child2, dipole;
    dipole = tempJets_two[cs_history[ca_tree_node_hardest].jetp_index];
    dipole.has_parents(child1, child2);
    secondary = getHardBranch(child1, child2);
  }

  // Check which splitting is the next to hardest in both planes

  else if (vector_tuple_slp.size()>0 && vector_tuple_lp.size()>1){
    double kappa_nxthardest_lp =  get<0>(vector_tuple_lp[1]);
    double kappa_nxthardest_slp = get<0>(vector_tuple_slp[0]);

     // next to hardest is on the primary lund plane
    if(kappa_nxthardest_lp < kappa_nxthardest_slp){

      int ca_tree_node_nxthardest = get<2>(vector_tuple_lp[1]);

      if(ca_tree_node_nxthardest > ca_tree_node_hardest){

      primary = tempJets_two[cs_history[ca_tree_node_hardest].jetp_index];
      fastjet::PseudoJet child1, child2, dipole;
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
     // Next-to-hardest is on the secondary lund plane
    else{
      //  cout << "eo" << endl;
      int ca_tree_node_nxthardest = get<2>(vector_tuple_slp[0]);
      primary = tempJets_two[cs_history[ca_tree_node_nxthardest].jetp_index];
      fastjet::PseudoJet child1, child2, dipole;
      dipole = tempJets_two[cs_history[ca_tree_node_hardest].jetp_index];
      dipole.has_parents(child1, child2);
      secondary = getHardBranch(child1, child2);
    }

  }

    // there is no next-to-hardest splitting
  else{

    int ca_tree_node_hardest = get<2>(vector_tuple_lp[0]);
    primary =  tempJets_two[cs_history[ca_tree_node_hardest].jetp_index];;
    secondary = fastjet::PseudoJet(0.,0.,0.,0.);

  }
   fastjet::PseudoJet groomed_jet = join(primary, secondary);
   fjOutputs_.push_back(groomed_jet);

   // Compute the n-subjettiness ratio
   if (groomed_jet.pt()>0){
   double beta = 2;

   fastjet::contrib::NsubjettinessRatio nSub32_beta2(3,2, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));

   fastjet::contrib::NsubjettinessRatio nSub21_beta2(2,1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));

   double tau32_beta2 = nSub32_beta2(groomed_jet);
   double tau21_beta2 = nSub21_beta2(groomed_jet);

   tau32_.push_back(tau32_beta2);
   tau21_.push_back(tau21_beta2);
}
else{
  tau32_.push_back(0);
  tau21_.push_back(0);
}

} // if vector_tuple.size() ==0
   else fjOutputs_.push_back(fastjet::PseudoJet(0.,0.,0.,0.));
 }
   return fjOutputs_;
}


#endif

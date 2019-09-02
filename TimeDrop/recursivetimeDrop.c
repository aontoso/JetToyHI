

#include <iostream>
#include <vector>
#include <algorithm> // std::min_element
#include <iterator>
#include <tuple>

void recursivetimeDrop(){

TFile *file = TFile::Open("/Users/albasotoontoso/Work/Jet_substraction/JetToyHI/TimeDrop/JetToyHIResultSD_PYTHIA_pp_13TeV_hadronLevel_recursive_ps_timeDrop.root");
TTree *treef = (TTree *)file->Get("jetTree");
Long64_t nevents_pp = treef->GetEntriesFast();


Long64_t njets_real_pp = 0;

TH1D* h_mass_infinity = new TH1D("mass spectrum inf", "mass spectrum inf", 50., -8. ,0.);



std::vector<double> * sigJetPt = 0;
std::vector<double> * sigJetEta = 0;
std::vector<double> * sigJetPhi = 0;
std::vector<double> * sigJetM = 0;
std::vector<double> * sigJetArea = 0;

std::vector<double> * sdJetPt = 0;
std::vector<double> * sdJetEta = 0;
std::vector<double> * sdJetPhi = 0;
std::vector<double> * sdJetM = 0;
std::vector<double> * sdJetArea = 0;
//std::vector<int> * ndropped = 0;

std::vector<std::vector<double>> * tf = 0;
std::vector<std::vector<double>> * deltaR = 0;
std::vector<std::vector<double>> * zg = 0;
std::vector<std::vector<double>> * kt = 0;
std::vector<std::vector<double>> * mass = 0;
std::vector<std::vector<double>> * pt = 0;
std::vector<std::vector<double>> * dipolepx_p = 0;
std::vector<std::vector<double>> * dipolepy_p = 0;
std::vector<std::vector<double>> * dipolepz_p = 0;
std::vector<std::vector<double>> * dipoleE_p = 0;
std::vector<std::vector<double>> * dipolepx_s = 0;
std::vector<std::vector<double>> * dipolepy_s = 0;
std::vector<std::vector<double>> * dipolepz_s = 0;
std::vector<std::vector<double>> * dipoleE_s = 0;

std::vector<double> * eventWeight = 0;
std::vector<double> * nsd = 0;
// Real Signal JET

treef->SetBranchAddress("sigJetPt",&sigJetPt);
treef->SetBranchAddress("sigJetEta",&sigJetEta);
treef->SetBranchAddress("sigJetPhi",&sigJetPhi);
treef->SetBranchAddress("sigJetM",&sigJetM);

// SoftDropped Jet

treef->SetBranchAddress("sdJetPt",&sdJetPt);
treef->SetBranchAddress("sdJetPhi",&sdJetPhi);
treef->SetBranchAddress("sdJetM",&sdJetM);
treef->SetBranchAddress("sdJetArea",&sdJetArea);

treef->SetBranchAddress("tf",&tf);
treef->SetBranchAddress("deltaR",&deltaR);
treef->SetBranchAddress("zg",&zg);
treef->SetBranchAddress("kt",&kt);
treef->SetBranchAddress("mass", &mass);
treef->SetBranchAddress("pt", &pt);
treef->SetBranchAddress("dipolepx_p",&dipolepx_p);
treef->SetBranchAddress("dipolepy_p",&dipolepy_p);
treef->SetBranchAddress("dipolepz_p", &dipolepz_p);
treef->SetBranchAddress("dipoleE_p", &dipoleE_p);
treef->SetBranchAddress("dipolepx_s",&dipolepx_s);
treef->SetBranchAddress("dipolepy_s",&dipolepy_s);
treef->SetBranchAddress("dipolepz_s", &dipolepz_s);
treef->SetBranchAddress("dipoleE_s", &dipoleE_s);

treef->SetBranchAddress("eventWeight",&eventWeight);
treef->SetBranchAddress("nsd",&nsd);

int njets_pp_above = 0;
double ptcut = 450;
double factorR2 = 0.8*0.8; // I have to divide because of how I created the tree
double factorR = 0.8;

// We define "a" as a power of theta that selects the variable that we want

double a = 2; // timeDrop. if a=0 zDrop. if a=1 kTDrop.

for (Long64_t k=0; k<nevents_pp; k++){ // Loop sobre el total de entradas del trees
      treef->GetEntry(k);
      njets_real_pp = sigJetPt->size();
      for (int i=0; i<njets_real_pp; i++){
        double pt_value_sig = sigJetPt->at(i);
  //      cout << pt_value_sig<<endl;
        if(pt_value_sig > ptcut){
          njets_pp_above++;
          int n_splittings = nsd->at(i);

if(n_splittings>0){
        tuple <double, int> info_dropvariable;
        std::vector<tuple<double, int>> dropvariable_and_j;
       for(int j=0; j<n_splittings; j++){
      //   cout << deltaR->at(i).at(j) << endl;
      double delta_r = deltaR->at(i).at(j);
      double momentum_share = zg->at(i).at(j);
      double piti = pt->at(i).at(j);
      double drop_variable = 1/(momentum_share*(1-momentum_share)*piti*pow(delta_r,a));

       info_dropvariable = make_tuple(drop_variable, j);
       dropvariable_and_j.push_back(info_dropvariable);
    } // Loop over splittings

    //Sort the  list of tuples in the dropped variable keeping the original position in the CA declustering fixed
    std::vector<int> accepted_branches;
    sort(dropvariable_and_j.begin(), dropvariable_and_j.end());

    // Find which branches are angular and dropped variable ordered
    int deltar_maximal = get<1>(dropvariable_and_j[0]);
    accepted_branches.push_back(deltar_maximal);

     for (int i = 0; i < dropvariable_and_j.size(); i++) {
         int deltar_subsequent = get<1>(dropvariable_and_j[i]);
         if(deltar_subsequent > deltar_maximal){ accepted_branches.push_back(deltar_subsequent);
         deltar_maximal = deltar_subsequent;
     }
   }

       double secondaries_energy = 0;
       double secondaries_px = 0;
       double secondaries_py = 0;
       double secondaries_pz = 0;
     for (int k = 0; k < accepted_branches.size(); k++){
       int position = accepted_branches.at(k);
       secondaries_energy+=dipoleE_s->at(i).at(position);
       secondaries_px+=dipolepx_s->at(i).at(position);
       secondaries_py+=dipolepy_s->at(i).at(position);
       secondaries_pz+=dipolepz_s->at(i).at(position);
   }
   int position = accepted_branches.at(accepted_branches.size()-1);
   double primary_energy = dipoleE_p->at(i).at(position);
   double primary_px = dipolepx_p->at(i).at(position);
   double primary_py = dipolepy_p->at(i).at(position);
   double primary_pz = dipolepz_p->at(i).at(position);

   double groomed_energy = primary_energy+secondaries_energy;
   double groomed_px = primary_px+secondaries_px;
   double groomed_py = primary_py+secondaries_py;
   double groomed_pz = primary_pz+secondaries_pz;
   double groomed_mass_dipole = sqrt(pow(groomed_energy,2)-pow(groomed_px,2)-pow(groomed_py,2)-pow(groomed_pz,2));
   double rho_inf = pow(groomed_mass_dipole,2)/(pow(sigJetPt->at(i),2)*pow(factorR,2));
  h_mass_infinity->Fill(log10(rho_inf));

} // close the if (nsplittings>-)condition


     } // close ptcut condition
    } // Loop over jets

  } // Loop over events

h_mass_infinity->Sumw2();
h_mass_infinity->Scale(1./(h_mass_infinity->GetBinWidth(0)*njets_pp_above));

  TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 800, 600);
  Int_t azul;
  azul = TColor::GetColor("#034F84");
  Int_t rojo;
  rojo = TColor::GetColor("#B93A32");
  Int_t gris;
  gris = TColor::GetColor("#838487");

  TPaveText *pt1 = new TPaveText(0.2,0.7,0.3,1,"blNDC");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->SetFillStyle(0);
  pt1->SetTextSize(0.04);
  TText *AText1 = pt1->AddText("PYTHIA8 p+p @13 TeV");
  pt1->Draw();
  TPaveText *pt2 = new TPaveText(0.1,0.7,0.3,1,"blNDC");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  pt2->SetFillStyle(0);
  pt2->SetTextSize(0.04);
  TText *AText2 = pt2->AddText("p_{T,jet}>450 GeV/c, AK08, CA1");

  TPaveText *pt4 = new TPaveText(0.4,0.7,0.3,1,"blNDC");
  pt4->SetBorderSize(0);
  pt4->SetFillColor(0);
  pt4->SetFillStyle(0);
  pt4->SetTextSize(0.04);
  TText *AText4 = pt4->AddText("MPI off, ISR off, HAD off");

  h_mass_infinity->GetXaxis()->SetTitle("log_{10}#rho");
  h_mass_infinity->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dlog_{10}#rho}");

  h_mass_infinity->SetTitle("");
  h_mass_infinity->GetYaxis()->SetLabelSize(0.04);
  h_mass_infinity->GetXaxis()->SetLabelSize(0.04);
  h_mass_infinity->GetYaxis()->SetTitleSize(0.06);
  h_mass_infinity->GetXaxis()->SetTitleSize(0.06);
  //h_zg_pp->SetStats(0);
  h_mass_infinity->SetMarkerColor(azul);
  h_mass_infinity->SetLineColor(azul);
  h_mass_infinity->SetMarkerStyle(20);
  h_mass_infinity->SetMarkerSize(1.4);

  h_mass_infinity->Draw("HIST");

  pt2->Draw();
  pt1->Draw();
  pt4->Draw();

}

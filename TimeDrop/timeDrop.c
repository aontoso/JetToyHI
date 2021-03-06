

#include <iostream>
#include <vector>
#include <algorithm> // std::min_element
#include <iterator>
#include <tuple>

void timeDrop(){

TFile *file = TFile::Open("/Users/albasotoontoso/Work/Jet_substraction/JetToyHI/TimeDrop/JetToyHIResultSD_PYTHIA_pp_13TeV_hadronLevel_recursive_ps_timeDrop.root");
TTree *treef = (TTree *)file->Get("jetTree");
Long64_t nevents_pp = treef->GetEntriesFast();


Long64_t njets_real_pp = 0;
TH1D* h_zg_pp = new TH1D("zg spectrum pp", "zg spectrum pp", 50., 0. ,0.5);

TH1D* h_mass_pp = new TH1D("mass spectrum pp", "mass spectrum pp", 50., -10. ,0.);

TH1D* h_mass_sd = new TH1D("mass spectrum sd", "mass spectrum sd", 50., -10. ,0.);

TH1D* h_mass_tag = new TH1D("mass spectrum tag", "mass spectrum tag", 50., -6. ,0.);

TH1I* h_ntot = new TH1I("ntot pp", "ntot pp", 20., 0.,20.);

TH1I* h_nsd = new TH1I("nsd pp", "nsd pp", 20., 0.,20.);

TH1I* h_ndropped = new TH1I("ndropped pp", "ndropped pp", 20., 0.,20.);

TH2D* h_lund = new TH2D("lund plane pp", "lund plane pp", 25., 0.,7.,25.,-4.,6.);

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

double a = 0; // timeDrop. if a=0 zDrop. if a=1 kTDrop.

for (Long64_t k=0; k<nevents_pp; k++){ // Loop sobre el total de entradas del trees
      treef->GetEntry(k);
      njets_real_pp = sigJetPt->size();
      for (int i=0; i<njets_real_pp; i++){
        double pt_value_sig = sigJetPt->at(i);
  //      cout << pt_value_sig<<endl;
        if(pt_value_sig > ptcut){
          njets_pp_above++;
          int n_splittings = nsd->at(i);
      //    cout << "Jet: " << i << endl;

if(n_splittings>0){
       double min_drop_variable = 1e4;
        int jmin = -1;

        h_ntot->Fill(n_splittings);
        int nsd_pp = 0;
       for(int j=0; j<n_splittings; j++){
      //   cout << deltaR->at(i).at(j) << endl;
       double delta_r = deltaR->at(i).at(j);
       double momentum_share = zg->at(i).at(j);
       double piti = pt->at(i).at(j);
       double drop_variable = 1/(momentum_share*(1-momentum_share)*piti*pow(delta_r,a));

       if(drop_variable < min_drop_variable) {
         //next_to_shortest = min_formation_time;
         min_drop_variable = drop_variable;
        jmin = j;
       }

    } // Loop over splittings

    h_nsd->Fill(nsd_pp);
    //////////////////////////////////////////////////////////////////
    // Compute the Lund plane
    ///////////////////////////////////////////////////////////////////
    double deltar_shortest_drop_variable = deltaR->at(i).at(jmin);
    double kt_shortest_drop_variable = kt->at(i).at(jmin)*factorR;
    int ndropped = 0;
    for(int j=0; j<n_splittings; j++){

     if(deltaR->at(i).at(j) > deltar_shortest_drop_variable) ndropped++;

    }
    h_ndropped->Fill(ndropped);
    h_lund->Fill(log(1/deltar_shortest_drop_variable),log(kt_shortest_drop_variable));

    //////////////////////////////////////////////////////////////////
    // Compute the plain, tagged and groomed mass
    ///////////////////////////////////////////////////////////////////

    double momentum_share = zg->at(i).at(jmin);
    double mass_subjets = mass->at(i).at(jmin);

    double theta= deltaR->at(i).at(jmin);
    double piti = pt->at(i).at(jmin);
    double mass_tagging = momentum_share*(1-momentum_share)*pow(piti,2)*pow(theta,2);
    double rho_sd = pow(mass_subjets,2)/(pow(sigJetPt->at(i),2)*pow(factorR,2));
    double rho = pow(sigJetM->at(i),2)/(pow(sigJetPt->at(i),2)*pow(factorR,2));
    double rho_tag = mass_tagging/(pow(sigJetPt->at(i),2)*pow(factorR,2));

    h_mass_pp->Fill(log10(rho));
    h_mass_sd->Fill(log10(rho_sd));
    h_mass_tag->Fill(log10(rho_tag));
    h_zg_pp->Fill(momentum_share);

  }

} // close ptcut condition
    } // Loop over jets

  } // Loop over events
//  cout << njets_pp_above << endl;
  h_zg_pp->Sumw2();
  h_zg_pp->Scale(1./(h_zg_pp->GetBinWidth(0)*njets_pp_above));
   h_mass_pp->Sumw2();
  h_mass_pp->Scale(1./(h_mass_pp->GetBinWidth(0)*njets_pp_above));
  h_mass_sd->Sumw2();
 h_mass_sd->Scale(1./(h_mass_sd->GetBinWidth(0)*njets_pp_above));
 h_mass_tag->Sumw2();
 h_mass_tag->Scale(1./(h_mass_tag->GetBinWidth(0)*njets_pp_above));

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

  TPaveText *pt3 = new TPaveText(0.5,0.7,0.3,1,"blNDC");
  pt3->SetBorderSize(0);
  pt3->SetFillColor(0);
  pt3->SetFillStyle(0);
  pt3->SetTextSize(0.04);
  //TText *AText3 = pt3->AddText("log(k_{T})>0");

  h_zg_pp->GetXaxis()->SetTitle("z");
  h_zg_pp->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dz}");

  h_zg_pp->SetTitle("");
  c1->SetLogy();
  //c1->SetLogx();
  h_zg_pp->GetYaxis()->SetLabelSize(0.04);
  h_zg_pp->GetXaxis()->SetLabelSize(0.04);
  h_zg_pp->GetYaxis()->SetTitleSize(0.06);
  h_zg_pp->GetXaxis()->SetTitleSize(0.06);
  //h_zg_pp->SetStats(0);
  h_zg_pp->SetMarkerColor(azul);
  h_zg_pp->SetLineColor(azul);
  h_zg_pp->SetMarkerStyle(20);
  h_zg_pp->SetMarkerSize(1.4);

  h_zg_pp->Draw("HIST");

  pt2->Draw();
  pt1->Draw();
//  pt3->Draw();
  pt4->Draw();

  TCanvas *c2 = new TCanvas ("c2", "c2", 65, 52, 800, 600);
//  c2->SetLogy();
//  c2->SetLogx();
  h_mass_sd->GetXaxis()->SetTitle("log_{10}#rho");
  h_mass_sd->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dlog_{10}#rho}");

  h_mass_sd->SetTitle("");
  h_mass_sd->GetYaxis()->SetLabelSize(0.04);
  h_mass_sd->GetXaxis()->SetLabelSize(0.04);
  h_mass_sd->GetYaxis()->SetTitleSize(0.06);
  h_mass_sd->GetXaxis()->SetTitleSize(0.06);
  h_mass_sd->SetMarkerColor(azul);

  h_mass_sd->SetLineColor(azul);
  h_mass_sd->SetMarkerStyle(20);
  h_mass_sd->SetMarkerSize(1.4);

  h_mass_sd->Draw("HIST");

  pt2->Draw();
  pt1->Draw();
//  pt3->Draw();
  pt4->Draw();

  TCanvas *c3 = new TCanvas ("c3", "c3", 65, 52, 800, 600);

  h_nsd->GetXaxis()->SetTitle("n");
  h_nsd->GetYaxis()->SetTitle("#frac{dN}{dn}");

  h_nsd->SetTitle("");
  h_nsd->GetYaxis()->SetLabelSize(0.04);
  h_nsd->GetXaxis()->SetLabelSize(0.04);
  h_nsd->GetYaxis()->SetTitleSize(0.06);
  h_nsd->GetXaxis()->SetTitleSize(0.06);
  h_nsd->SetStats(0);
  h_ntot->SetMarkerColor(azul);
  h_ntot->SetLineColor(rojo);
  h_nsd->SetLineColor(azul);
  h_nsd->SetLineColor(kGreen+2);
  h_nsd->SetMarkerStyle(20);
  h_nsd->SetMarkerSize(1.4);
  //h_ntot->Draw("HIST");
  //h_nsd->Draw("SAME HIST");
  h_nsd->Draw("HIST");

   auto leg = new TLegend(0.6,0.6,0.75,0.75);

   leg->AddEntry(h_ntot, "n_{tot}", "l");
  leg->AddEntry(h_nsd, "n_{SD}", "l");
  leg->AddEntry(h_ndropped, "n_{TDropped}", "l");
  //  leg->AddEntry(h_mass_tag, "SD tagging mode", "l");
    leg->SetTextSize(0.05);
    leg->Draw();

    pt2->Draw();
    pt1->Draw();
  //  pt3->Draw();
    pt4->Draw();

    TCanvas *c4 = new TCanvas ("c4", "c4", 65, 52, 800, 600);
    h_lund->GetYaxis()->SetTitle("ln(k_{T})[GeV/c]");
    h_lund->GetXaxis()->SetTitle("ln(1/#Delta R)");

    h_lund->SetTitle("");
    h_lund->GetYaxis()->SetLabelSize(0.04);
    h_lund->GetXaxis()->SetLabelSize(0.04);
    h_lund->GetYaxis()->SetTitleSize(0.06);
    h_lund->GetXaxis()->SetTitleSize(0.06);
    h_lund->Draw("COLZ");

}

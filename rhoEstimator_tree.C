
void rhoEstimator_tree(){

TFile *file =TFile::Open("/Users/albasotoontoso/Work/Jet_substraction/JetToyHI/Bkg_resolution/Figures6and7/JetToyHIResultRhoEstimator_ModifiedSK_JEWEL_Medium_deltaRho2_plusSK.root");
TTree *treef = (TTree *)file->Get("jetTree");
Long64_t nevents = treef->GetEntriesFast();
Long64_t npatches = 0;
Long64_t npatches_sk = 0;

TH1D* h_ptbkg_response_median = new TH1D("p_{T} bkg response", "p_{T} bkg response", 50., -80. , 70);
TH1D* h_ptbkg_response_rhoN = new TH1D("Intersection-#rho(N)", "Intersection-#rho(N)", 50., -80. , 70);
TH1D* h_ptbkg_response_rholine = new TH1D("p_{T} bkg response THREE", "p_{T} bkg response THREE", 50., -80. , 70);
//TH1D* h_ptbkg_response_rhosk = new TH1D("p_{T} bkg response FOUR", "p_{T} bkg response FOUR", 50., -80. , 70);

TH1D* h_ptsig_response_median = new TH1D("p_{T} sig response", "p_{T} bkg response", 50., -80. , 80);
TH1D* h_ptsig_response_rhoN = new TH1D("p_{T} sig response TWO", "p_{T} sig response TWO", 50., -80. , 80.);
TH1D* h_ptsig_response_rholine = new TH1D("p_{T} sig response THREE", "p_{T} sig response THREE", 50., -80. , 80);
TH1D* h_ptsig_response_rhosk = new TH1D("p_{T} sig response FOUR", "p_{T} sig response FOUR", 50., -100. , 80);

TH1D* h_gamma_bkg = new TH1D("Gamma bkg", "gamma bkg", 100., 0. , 60);
TH1D* h_gamma_sig = new TH1D("Gamma sig", "gamma sig", 100., -60. , 120.);
TH1D* h_deltaR_sig = new TH1D("Delta R sig", "Delta R sig", 20., 0. , 0.6);
TH1D* h_deltaR_bkg = new TH1D("Delta R bkg", "Delta R sig", 20., 0. , 0.6);

std::vector<double> * rho_median = 0;
std::vector<double> * rho_True = 0;
std::vector<double> * pT_patch = 0;
std::vector<double> * rho_Estimate = 0;
std::vector<double> * nConstituents_True = 0;
std::vector<double> * rho_line = 0;
std::vector<double> * rho_SK = 0;

std::vector<double> * bkgJetPt = 0;
std::vector<double> * bkgJetEta = 0;
std::vector<double> * bkgJetPhi = 0;
std::vector<double> * bkgJetM = 0;
std::vector<double> * bkgJetArea = 0;

std::vector<std::vector<double>> * bkgJetConstPt = 0;
std::vector<std::vector<double>> * bkgJetConstEta = 0;
std::vector<std::vector<double>> * bkgJetConstPhi = 0;
std::vector<std::vector<double>> * bkgJetConstM = 0;
std::vector<std::vector<double>> * bkgJetConstId = 0;
std::vector<std::vector<double>> * bkgJetConstDeltaR = 0;

treef->SetBranchAddress("rho_True",&rho_True);
treef->SetBranchAddress("pT_patch",&pT_patch);
treef->SetBranchAddress("rho_median",&rho_median);
treef->SetBranchAddress("nConstituents_True", &nConstituents_True);
treef->SetBranchAddress("rho_Estimate",&rho_Estimate);
treef->SetBranchAddress("rho_line",&rho_line);

treef->SetBranchAddress("bkgJetPt",&bkgJetPt);
treef->SetBranchAddress("bkgJetEta",&bkgJetEta);
treef->SetBranchAddress("bkgJetPhi",&bkgJetPhi);
treef->SetBranchAddress("bkgJetM",&bkgJetM);
treef->SetBranchAddress("bkgJetArea",&bkgJetArea);
treef->SetBranchAddress("bkgJetConstPt",&bkgJetConstPt);
treef->SetBranchAddress("bkgJetConstEta",&bkgJetConstEta);
treef->SetBranchAddress("bkgJetConstPhi",&bkgJetConstPhi);
treef->SetBranchAddress("bkgJetConstM",&bkgJetConstM);
treef->SetBranchAddress("bkgJetConstId",&bkgJetConstId);
treef->SetBranchAddress("bkgJetConstDeltaR",&bkgJetConstDeltaR);
treef->SetBranchAddress("rho_SK",&rho_SK);

Long64_t nconstituents = 0;

for (Long64_t k=0; k<nevents; k++){ // Loop sobre el total de entradas del trees
      treef->GetEntry(k);
      npatches = rho_True->size();

    //  cout << bkgJetPt->size() << endl;
     npatches_sk = rho_SK->size();
      for (int i=0; i<npatches; i++){
        double pt_sig_true = pT_patch->at(i)-rho_True->at(i);
      //  if (i==0) cout <<rho_True->at(i) <<endl;
        if(pt_sig_true > 120){
         double pt_signal_reco_estimate = pT_patch->at(i)-rho_Estimate->at(i);
      //    double pt_signal_reco_estimate = rho_Estimate->at(i);
          double pt_signal_reco_median = pT_patch->at(i)-rho_median->at(i);
          double pt_signal_reco_line = pT_patch->at(i)-rho_line->at(i);
    //      double pt_signal_reco_sk = rho_SK->at(i);
          //if (pt_signal_reco_estimate == 0) cout << "eo" << endl;
         h_ptsig_response_median->Fill(pt_signal_reco_median-pt_sig_true);
        h_ptsig_response_rhoN->Fill(pt_signal_reco_estimate-pt_sig_true);
    //     h_ptsig_response_rhoN->Fill(pt_signal_reco_estimate);
         h_ptsig_response_rholine->Fill(pt_signal_reco_line-pt_sig_true);
  //       h_ptsig_response_rhosk->Fill(pt_signal_reco_sk-pt_sig_true);

      //  nconstituents = bkgJetConstPt->at(i).size();
        double ptSignal = 0;
    //    for (int j=0; j<nconstituents; j++){
    //         int iD = bkgJetConstId->at(i).at(j);
    //         double momentum = bkgJetConstPt->at(i).at(j);
    //         double deltaR = bkgJetConstDeltaR->at(i).at(j);
          //   cout << momentum << endl;
      //       double lambda = 0.5;
      //       double gamma = momentum/pow(deltaR,lambda);
          //   cout << iD << endl;
          //   if (iD == -1 && momentum < 2.4) {ptSignal+=momentum;
          //   h_deltaR_sig->Fill(deltaR);
          //   }
    //         if (iD == 0 && momentum < 2.4) {ptSignal+=momentum;
    //         h_deltaR_sig->Fill(deltaR);
    //         }
    //         if (abs(iD) == 1 && momentum < 2.4) {
    //         h_deltaR_bkg->Fill(deltaR);
    //         }
             //else h_gamma_bkg->Fill(gamma);

            // h_ptconst_shared->Fill(pt_const_value_shared);
      //    } // loop over anti-kt constituents
        //  cout << ptSignal << endl;
    //      h_gamma_sig->Fill(ptSignal);
       }

       double pt_bkg_true = rho_True->at(i);
    // if(pt_bkg_true>100){
       double pt_bkg_reco_median = rho_median->at(i);
       double pt_bkg_reco_estimate = rho_Estimate->at(i);
       double pt_bkg_reco_line = rho_line->at(i);
    //   if (i==0){
    //     cout << rho_True->at(i) << endl;
    //     cout << rho_line->at(i) << endl;
    //   }
    //   if(pt_bkg_reco_line-pt_bkg_reco_estimate>1e-8)
  //  if (abs(pt_bkg_reco_estimate-pt_bkg_true)>0){
  //  cout << pt_bkg_reco_estimate-pt_bkg_true << endl;
    h_ptbkg_response_rhoN->Fill(pt_bkg_reco_estimate-pt_bkg_true);
    //}

  //     cout << rho_Estimate->at(0) << endl;
    h_ptbkg_response_rholine->Fill(pt_bkg_reco_line-pt_bkg_true);
    h_ptbkg_response_median->Fill(pt_bkg_reco_median-pt_bkg_true);


        }

       //}

    }


   TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 800, 600);
    Int_t azul;
    azul = TColor::GetColor("#034F84");
    Int_t rojo;
    rojo = TColor::GetColor("#B93A32");
    Int_t gris;
    gris = TColor::GetColor("#838487");
    // c1->SetLogy();
     h_ptbkg_response_rhoN->GetXaxis()->SetTitle("p_{T,reco}^{bkg}-p_{T,truth}^{bkg}[GeV/c]");
    //  h_ptbkg_response_rhoN->GetXaxis()->SetTitle("#delta p_{T}[GeV/c]");
      h_ptbkg_response_rhoN->GetYaxis()->SetTitle("Probability density");
      h_ptbkg_response_rhoN->SetTitle("");
      h_ptbkg_response_rhoN->GetYaxis()->SetLabelSize(0.04);
      h_ptbkg_response_rhoN->GetXaxis()->SetLabelSize(0.04);
      h_ptbkg_response_rhoN->GetYaxis()->SetTitleSize(0.06);
      h_ptbkg_response_rhoN->GetXaxis()->SetTitleSize(0.06);
    //  h_ptbkg_response_median->GetYaxis()->SetRangeUser(0.,0.2);
    //  h_ptbkg_response_rhoN->SetStats(0);
      h_ptbkg_response_median->SetMarkerColor(azul);
      h_ptbkg_response_median->SetLineColor(azul);
      h_ptbkg_response_median->SetMarkerStyle(20);
      h_ptbkg_response_rhoN->SetMarkerColor(rojo);
      h_ptbkg_response_rhoN->SetLineColor(rojo);
      h_ptbkg_response_rhoN->SetMarkerStyle(21);
      h_ptbkg_response_rholine->SetMarkerColor(gris);
      h_ptbkg_response_rholine->SetLineColor(gris);
       h_ptbkg_response_rholine->SetMarkerStyle(21);
      h_ptbkg_response_rhoN->DrawNormalized("HIST");

  //     h_ptbkg_response_median->DrawNormalized("SAME HIST");
  //    h_ptbkg_response_rholine->DrawNormalized("SAME HIST");


      auto leg3 = new TLegend(0.6,0.6,0.75,0.75);
    //  leg3->AddEntry(h_ptbkg_response_rhoN, "YASM+#delta#rho^{sig}_{<}", "l");
      leg3->AddEntry(h_ptbkg_response_rhoN, "med(#rho)A", "l");
  //    leg3->AddEntry(h_ptbkg_response_rholine, "#rho_{corr}(N^{truth})", "l");
      leg3->SetTextSize(0.05);
      leg3->Draw("SAME");
      TPaveText *pt1 = new TPaveText(0.1,0.7,0.3,1,"blNDC");
      pt1->SetBorderSize(0);
      pt1->SetFillColor(0);
      pt1->SetFillStyle(0);
      pt1->SetTextSize(0.04);
      TText *AText1 = pt1->AddText("p^{SK}_{T,cut}");
    //  pt1->Draw();


 TCanvas *c2 = new TCanvas ("c2", "c2", 65, 52, 800, 600);
  c2->SetLogy();
  h_ptsig_response_rhoN->GetXaxis()->SetTitle("p_{T,reco}^{jet}-p_{T,truth}^{jet}[GeV/c]");
//  h_ptsig_response_rhoN->GetXaxis()->SetTitle("p_{T,<}^{sig}[GeV/c]");
  h_ptsig_response_rhoN->GetYaxis()->SetTitle("Probability density");
  h_ptsig_response_rhoN->SetTitle("");
  h_ptsig_response_rhoN->GetYaxis()->SetLabelSize(0.04);
  h_ptsig_response_rhoN->GetXaxis()->SetLabelSize(0.04);
  h_ptsig_response_rhoN->GetYaxis()->SetTitleSize(0.06);
  h_ptsig_response_rhoN->GetXaxis()->SetTitleSize(0.06);
//   h_ptsig_response_rhoN->SetStats(0);
  h_ptsig_response_median->SetMarkerColor(azul);
  h_ptsig_response_median->SetLineColor(azul);
  h_ptsig_response_median->SetMarkerStyle(20);
  h_ptsig_response_rhoN->SetMarkerColor(rojo);
  h_ptsig_response_rhoN->SetLineColor(rojo);
  h_ptsig_response_rhoN->SetMarkerStyle(21);

  h_ptsig_response_rhosk->SetMarkerColor(gris);
  h_ptsig_response_rhosk->SetLineColor(gris);
  h_ptsig_response_rhosk->SetMarkerStyle(21);

   h_ptsig_response_rhoN->DrawNormalized("HIST");
  // h_ptsig_response_median->DrawNormalized("SAME HIST");
  // h_ptsig_response_rhosk->DrawNormalized("SAME HIST");

  auto leg4 = new TLegend(0.6,0.6,0.75,0.75);
  leg4->AddEntry(h_ptsig_response_rhoN, "YAM,p_{T}^{cut}=2#mu", "l");
  leg4->AddEntry(h_ptsig_response_median, "med(#rho)A", "l");
//  leg4->AddEntry(h_ptsig_response_rholine, "#rho_{corr}(N^{truth})", "l");

  leg4->SetTextSize(0.05);
  leg4->Draw("SAME");

  TPaveText *pt6 = new TPaveText(0.1,0.4,0.3,1,"blNDC");
  pt6->SetBorderSize(0);
  pt6->SetFillColor(0);
  pt6->SetFillStyle(0);
  pt6->SetTextSize(0.04);
  TText *AText7 = pt6->AddText("<med(#rho)>=-0.56, #sigma(med(#rho))=15.44");
//  pt6->Draw();

  TPaveText *pt8 = new TPaveText(0.3,0.7,0.8,1,"blNDC");
  pt8->SetBorderSize(0);
  pt8->SetFillColor(0);
  pt8->SetFillStyle(0);
  pt8->SetTextSize(0.04);
  TText *AText8 = pt8->AddText("<#rho_{corr}(N^{cross})>=-13.02, #sigma(#rho_{corr}(N^{cross}))=16.4");
//  pt8->Draw();

  TPaveText *pt = new TPaveText(0.1,0.7,0.3,1,"blNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextSize(0.04);
  TText *AText = pt->AddText("p_{T,truth}^{jet}>120 GeV/c");
  pt->Draw();

  TPaveText *pt80 = new TPaveText(0.1,0.7,0.3,1,"blNDC");
  pt80->SetBorderSize(0);
  pt80->SetFillColor(0);
  pt80->SetFillStyle(0);
  pt80->SetTextSize(0.04);
  TText *AText80 = pt80->AddText("PYTHIA8");
  pt80->Draw();


//h_ptbkg_response_TWO->DrawNormalized();
//h_ptbkg_response->DrawNormalized("SAME");

TCanvas *c3 = new TCanvas ("c3", "c3", 65, 52, 800, 600);
// c3->SetLogy();
 h_gamma_sig->GetXaxis()->SetTitle("p_{T,<}^{sig}[GeV/c]");
 h_gamma_sig->GetYaxis()->SetTitle("Probability density");
 h_gamma_sig->SetTitle("");
 //h_gamma_sig->SetStats(0);
 h_gamma_sig->GetYaxis()->SetLabelSize(0.04);
 h_gamma_sig->GetXaxis()->SetLabelSize(0.04);
 h_gamma_sig->GetYaxis()->SetTitleSize(0.06);
 h_gamma_sig->GetXaxis()->SetTitleSize(0.06);
// h_gamma_bkg->SetMarkerColor(azul);
// h_gamma_bkg->SetLineColor(azul);
// h_gamma_bkg->SetMarkerStyle(20);
 h_gamma_sig->SetMarkerColor(rojo);
 h_gamma_sig->SetLineColor(rojo);
 h_gamma_sig->SetMarkerStyle(21);

  h_gamma_sig->DrawNormalized("HIST");
  pt80->Draw();
  pt->Draw();

  TCanvas *c4 = new TCanvas ("c4", "c4", 65, 52, 800, 600);
   c4->SetLogy();
   h_deltaR_sig->GetXaxis()->SetTitle("#DeltaR");
   h_deltaR_sig->GetYaxis()->SetTitle("Probability density");
   h_deltaR_sig->SetTitle("");
   h_deltaR_sig->SetStats(0);
   h_deltaR_sig->GetYaxis()->SetLabelSize(0.04);
   h_deltaR_sig->GetXaxis()->SetLabelSize(0.04);
   h_deltaR_sig->GetYaxis()->SetTitleSize(0.06);
   h_deltaR_sig->GetXaxis()->SetTitleSize(0.06);
   h_deltaR_bkg->SetMarkerColor(azul);
   h_deltaR_bkg->SetLineColor(azul);
   h_deltaR_bkg->SetMarkerStyle(20);
   h_deltaR_sig->SetMarkerColor(rojo);
   h_deltaR_sig->SetLineColor(rojo);
   h_deltaR_sig->SetMarkerStyle(21);

  //  h_deltaR_sig->DrawNormalized("HIST");
  //  h_deltaR_bkg->DrawNormalized("SAME HIST");
//   h_ptsig_response_median->DrawNormalized("SAME HIST");
//   h_ptsig_response_rholine->DrawNormalized("SAME HIST");

 auto leg5 = new TLegend(0.6,0.6,0.75,0.75);

 leg5->AddEntry(h_deltaR_bkg, "Background", "l");
 leg5->AddEntry(h_deltaR_sig, "Signal", "l");

 leg5->SetTextSize(0.05);
 leg5->Draw("SAME");

 TPaveText *pt21 = new TPaveText(0.1,0.7,0.3,1,"blNDC");
 pt21->SetBorderSize(0);
 pt21->SetFillColor(0);
 pt21->SetFillStyle(0);
 pt21->SetTextSize(0.04);
 TText *AText21 = pt21->AddText("p_{T,truth}^{jet}>120 GeV/c");
 pt21->Draw();

 TPaveText *pt800 = new TPaveText(0.1,0.7,0.3,1,"blNDC");
 pt800->SetBorderSize(0);
 pt800->SetFillColor(0);
 pt800->SetFillStyle(0);
 pt800->SetTextSize(0.04);
 TText *AText800 = pt800->AddText("p_{T,cut}=2#mu GeV/c");
 pt800->Draw();


}

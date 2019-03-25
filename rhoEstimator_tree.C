
void rhoEstimator_tree(){

TFile *file =TFile::Open("/Users/albasotoontoso/Work/Jet_substraction/JetToyHI/Macros_paperPLB/JetToyHIResultRhoEstimator_ModifiedSK_PYTHIA_3Mu_YAM.root");
TTree *treef = (TTree *)file->Get("jetTree");
Long64_t nevents = treef->GetEntriesFast();
Long64_t npatches = 0;
Long64_t npatches_sk = 0;

TH1D* h_ptbkg_response_median = new TH1D("p_{T} bkg response", "p_{T} bkg response", 50., -80. , 70);
TH1D* h_ptbkg_response_rhoN = new TH1D("Intersection-#rho(N)", "Intersection-#rho(N)", 50., -8. , 8.);
TH1D* h_ptbkg_response_rholine = new TH1D("p_{T} bkg response THREE", "p_{T} bkg response THREE", 50., -80. , 70);
//TH1D* h_ptbkg_response_rhosk = new TH1D("p_{T} bkg response FOUR", "p_{T} bkg response FOUR", 50., -80. , 70);

TH1D* h_ptsig_response_median = new TH1D("p_{T} sig response", "p_{T} bkg response", 50., -80. , 80);
TH1D* h_ptsig_response_rhoN = new TH1D("p_{T} sig response TWO", "p_{T} sig response TWO", 40., -120. , 120.);
TH1D* h_ptsig_response_rholine = new TH1D("p_{T} sig response THREE", "p_{T} sig response THREE", 50., -80. , 80);
TH1D* h_ptsig_response_rhosk = new TH1D("p_{T} sig response FOUR", "p_{T} sig response FOUR", 100., 0. , 2.);

TH1D* h_ptsig_response_SK = new TH1D("p_{T} sig response SK", "p_{T} sig response SK", 40., -120., 120.);

TH1D* h_gamma_bkg = new TH1D("Gamma bkg", "gamma bkg", 100., 0. , 500.);
TH1D* h_gamma_sig = new TH1D("Gamma sig", "gamma sig", 100., 0. , 300.);
TH1D* h_deltaR_sig = new TH1D("Delta R sig", "Delta R sig", 20., 0. , 0.6);
TH1D* h_deltaR_bkg = new TH1D("Delta R bkg", "Delta R sig", 20., 0. , 0.6);
TH1D* h_nsignal = new TH1D("N sig", "n sig", 100., 0. , 60);
TH1D* h_skPt = new TH1D("sk Pt", "sk Pt", 100., 0. , 6.);
TH1D* h_pt_sig = new TH1D("pT sig", "pT sig", 50., 0. , 100.);

TProfile *ptsignal_ptjet = new TProfile("hprof","Profile of psignal versus pt", 50., 0., 250.);
TProfile *ptsignal_ptabove = new TProfile("hprof3","Profile of psignal versus pTabove", 599., 1., 600., "s" );
TProfile *ptsignal_nabove = new TProfile("hproftwo","Profile of psignal versus nabove", 49., 1., 50., "s" );

TProfile2D *ptsignal_nabove_ptabove  = new TProfile2D("hprof2d","Profile of psignal versus nabove and ptabove",49.,1.,50.,599.,1.,600.,0.,200.);

//TH1D *h_rhoBelow_nBelow = new TH2D("h_rhoBelow_nBelow","Rho-n",80,0,180,80,0,180);
TH1D* h_rhoBelow = new TH1D("Rho Below", "Rho Below", 50., -80. , 80.);

std::vector<double> * rho_median = 0;
std::vector<double> * rho_True = 0;
std::vector<double> * pT_patch = 0;
std::vector<double> * rho_Estimate = 0;
std::vector<double> * nConstituents_True = 0;
std::vector<double> * rho_line = 0;
std::vector<double> * rho_SK = 0;
std::vector<double> * skPtThreshold = 0;
//std::vector<double> * rhoBelow_bkgd = 0;
//std::vector<double> * nparticles_below_bkgd = 0;

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

std::vector<double> * sigJetPt = 0;
std::vector<double> * sigJetEta = 0;
std::vector<double> * sigJetPhi = 0;
std::vector<double> * sigJetM = 0;
std::vector<double> * sigJetArea = 0;

std::vector<std::vector<double>> * sigJetConstPt = 0;
std::vector<std::vector<double>> * sigJetConstEta = 0;
std::vector<std::vector<double>> * sigJetConstPhi = 0;
std::vector<std::vector<double>> * sigJetConstM = 0;
std::vector<std::vector<double>> * sigJetConstId = 0;
std::vector<std::vector<double>> * sigJetConstDeltaR = 0;

std::vector<double> * skJetPt = 0;
std::vector<double> * skJetEta = 0;
std::vector<double> * skJetPhi = 0;
std::vector<double> * skJetM = 0;
std::vector<double> * skJetArea = 0;

std::vector<std::vector<double>> * skJetConstPt = 0;
std::vector<std::vector<double>> * skJetConstEta = 0;
std::vector<std::vector<double>> * skJetConstPhi = 0;
std::vector<std::vector<double>> * skJetConstM = 0;
std::vector<std::vector<double>> * skJetConstId = 0;
std::vector<std::vector<double>> * skJetConstDeltaR = 0;

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
treef->SetBranchAddress("skPtThreshold", &skPtThreshold);

treef->SetBranchAddress("sigJetPt",&sigJetPt);
treef->SetBranchAddress("sigJetEta",&sigJetEta);
treef->SetBranchAddress("sigJetPhi",&sigJetPhi);
treef->SetBranchAddress("sigJetM",&sigJetM);
treef->SetBranchAddress("sigJetArea",&sigJetArea);
treef->SetBranchAddress("sigJetConstPt",&sigJetConstPt);
treef->SetBranchAddress("sigJetConstEta",&sigJetConstEta);
treef->SetBranchAddress("sigJetConstPhi",&sigJetConstPhi);
treef->SetBranchAddress("sigJetConstM",&sigJetConstM);
treef->SetBranchAddress("sigJetConstId",&sigJetConstId);
treef->SetBranchAddress("sigJetConstDeltaR",&sigJetConstDeltaR);


treef->SetBranchAddress("skJetPt",&skJetPt);
treef->SetBranchAddress("skJetEta",&skJetEta);
treef->SetBranchAddress("skJetPhi",&skJetPhi);
treef->SetBranchAddress("skJetM",&skJetM);
treef->SetBranchAddress("skJetArea",&skJetArea);
treef->SetBranchAddress("skJetConstPt",&skJetConstPt);
treef->SetBranchAddress("skJetConstEta",&skJetConstEta);
treef->SetBranchAddress("skJetConstPhi",&skJetConstPhi);
treef->SetBranchAddress("skJetConstM",&skJetConstM);
treef->SetBranchAddress("skJetConstId",&skJetConstId);
treef->SetBranchAddress("skJetConstDeltaR",&skJetConstDeltaR);


Long64_t nconstituents = 0;
int nconstituents2 = 0;
int njets=0;

// Antes de hacer el analisis vamos a obtener el TProfile con signal por debajo de un cut en funcion de particulas por encima del cut

double mu = 1.2;
double pTcut_pTbelow = 3*mu;
double pTcut_pTabove = 6*mu;
double pTcut_nAbove = 6*mu;

for (Long64_t k=0; k<nevents; k++){ // Loop sobre el total de entradas del trees
      treef->GetEntry(k);
      npatches = rho_True->size();

      for (int i=0; i<npatches; i++){
        double pt_sig_true = pT_patch->at(i)-rho_True->at(i);
        if(pt_sig_true > 120){
       nconstituents2 = bkgJetConstPt->at(i).size();
       double ptSignal = 0;
       int nabove = 0;
       double ptAbove = 0;
       for (int j=0; j<nconstituents2; j++){
             int iD = bkgJetConstId->at(i).at(j);
             double momentum = bkgJetConstPt->at(i).at(j);
             if (momentum < pTcut_pTbelow && iD==0){
             ptSignal+=momentum;}
             if (momentum > pTcut_nAbove && iD==0) {
               nabove++;
             }
             if (momentum > pTcut_pTabove && iD==0) ptAbove+=momentum;
           }
      ptsignal_ptabove->Fill(ptAbove,ptSignal);
      ptsignal_nabove->Fill(nabove,ptSignal);
      ptsignal_nabove_ptabove->Fill(nabove,ptAbove,ptSignal);
  //    if (nabove == 6 && ptAbove > 180) cout << "yes" << ptAbove << endl;
    }
  }
}
//cout << ptsignal_nabove->GetBinContent(24) << endl;

for (Long64_t k=0; k<nevents; k++){ // Loop sobre el total de entradas del trees
      treef->GetEntry(k);
      npatches = rho_True->size();
  //    npatches_sk = skJetPt->size();
  //    for (int i = 0; i<npatches_sk; i++) {
  //     int nconst_True = sigJetConstPt->at(i).size();
  //     double scalarPt_True = 0;
  //     for(int j = 0; j<nconst_True; j++){
  //      scalarPt_True+=sigJetConstPt->at(i).at(j);
  //     }
  //     if(scalarPt_True > 120){
  //     int nconst_SK = skJetConstPt->at(i).size();
  //     double scalarPt_SK = 0;
  //     for(int j = 0; j<nconst_SK; j++){
  //      scalarPt_SK+=skJetConstPt->at(i).at(j);
  //     }
  //     h_ptsig_response_SK->Fill(scalarPt_SK-scalarPt_True);
    //   if (-7.3<scalarPt_SK-scalarPt_True && scalarPt_SK-scalarPt_True<7.3) njets++;
  //      }
  //  }
  //    for (int i = 0; i<npatches_sk; i++){
//         double pTSK = skJetPt->at(i);
//         double pTSig = sigJetPt->at(i);
  //     double pTSoftKiller = skPtThreshold->at(0);
  //       if(pTSig>120){
  //      h_ptsig_response_SK->Fill(pTSK-pTSig);
  //     h_skPt->Fill(pTSoftKiller);
   //    cout << pt_signal_reco_estimate << endl;
  //    }

//  }

  //   for (int i = 0; i<npatches; i++){
      //   double pTSoftKiller = skPtThreshold->at(0);
  //       double pTOurSoftKiller = rho_Estimate->at(0);
  //       h_skPt->Fill(pTOurSoftKiller);
  //    }
      for (int i=0; i<npatches; i++){
        double pt_sig_true = pT_patch->at(i)-rho_True->at(i);
        //cout << pt_sig_true << endl;

      //  if (pt_sig_true>10) h_gamma_bkg-> Fill(pt_sig_true);

      //  if (i==0) cout <<rho_True->at(i) <<endl;

  //         }
        if(pt_sig_true > 120){
        // Find the number of particles above the cut
        nconstituents = bkgJetConstPt->at(i).size();
        int nabove_jet = 0;
        double ptabove_jet = 0;
        for (int j=0; j<nconstituents; j++){
        double momentum = bkgJetConstPt->at(i).at(j);
        if (momentum > pTcut_pTabove) ptabove_jet+=momentum;
        if (momentum > pTcut_nAbove) nabove_jet++;
        }
        int ptabove_bin = int(ptabove_jet);
        int nabove_bin = nabove_jet;
        double pT_correction = ptsignal_ptabove->GetBinContent(ptabove_bin);
        double pT_correction_one = ptsignal_nabove->GetBinContent(nabove_bin);
        double pT_correction_two = ptsignal_nabove_ptabove->GetBinContent(nabove_bin,ptabove_bin);
      //  if (pT_correction_one != 0) cout << "N" << nabove_bin << endl;

         double pt_signal_reco_estimate = pT_patch->at(i)-rho_Estimate->at(i)+pT_correction_two;
    //      double pt_signal_reco_estimate = rho_Estimate->at(i);
          double pt_signal_reco_median = pT_patch->at(i)-rho_median->at(i);
          double pt_signal_reco_line = pT_patch->at(i)-rho_line->at(i);
          h_pt_sig->Fill(pT_correction_two);
          //cout << pT_correction_two << endl;
        //  if (pt_signal_reco_estimate-pt_sig_true > -10.4 && pt_signal_reco_estimate-pt_sig_true < 10.4) njets++;
    //    cout << pt_signal_reco_estimate << endl;

        // double pt_signal_reco_sk = rho_SK->at(i);
         h_ptsig_response_median->Fill(pt_signal_reco_median-pt_sig_true);
         h_ptsig_response_rhoN->Fill(pt_signal_reco_estimate-pt_sig_true);
      //   h_ptsig_response_rhoN->Fill(pt_signal_reco_estimate);
         h_ptsig_response_rholine->Fill(pt_signal_reco_line-pt_sig_true);
        // h_ptsig_response_rhosk->Fill(pt_signal_reco_sk-pt_sig_true);
         nconstituents = bkgJetConstPt->at(i).size();

        int nbelow2 = 0;
        double ptSignal2 = 0;
//        int nsignal = 0;
      for (int j=0; j<nconstituents; j++){
             int iD = bkgJetConstId->at(i).at(j);
             double momentum = bkgJetConstPt->at(i).at(j);
             if (momentum < 2.4 && iD==0){
             ptSignal2+=momentum;
            nbelow2++;
           }
           }
        //   cout << ptSignal2 << endl;
          h_gamma_sig->Fill(ptSignal2);

       }

    //   double pt_bkg_true = rho_True->at(i);
    // if(pt_bkg_true>100){
       //double pt_bkg_reco_median = rho_median->at(i);
       //double pt_bkg_reco_estimate = rho_Estimate->at(i);
      // double pt_bkg_reco_line = rho_line->at(i);

    //h_ptbkg_response_rhoN->Fill(pt_bkg_reco_estimate-pt_bkg_true);
    //}

  //     cout << rho_Estimate->at(0) << endl;
  //  h_ptbkg_response_rholine->Fill(pt_bkg_reco_line-pt_bkg_true);
//    h_ptbkg_response_median->Fill(pt_bkg_reco_median-pt_bkg_true);
  //h_ptbkg_response_sk->Fill(pt_bkg_reco_median-pt_bkg_true);

        }

    }
        // cout << ptsignal_nabove->GetBinError(2) << endl;
  //  h_ptsig_response_SK->Sumw2();
  //   double binWidth = h_ptsig_response_SK->GetBinWidth(0);
    // binWidth=1;
  //    h_ptsig_response_SK->Scale(1./((double)nevents*binWidth));

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
  h_ptsig_response_rhoN->SetMarkerColor(azul);
  h_ptsig_response_rhoN->SetLineColor(azul);
  h_ptsig_response_rhoN->SetMarkerStyle(20);
  h_ptsig_response_rhoN->SetMarkerColor(rojo);
  h_ptsig_response_rhoN->SetLineColor(rojo);
  h_ptsig_response_rhoN->SetMarkerStyle(21);

  h_ptsig_response_SK->SetMarkerColor(gris);
  h_ptsig_response_SK->SetLineColor(gris);
  h_ptsig_response_SK->SetMarkerStyle(21);

  // h_ptsig_response_rhoN->DrawNormalized("HIST");
  h_ptsig_response_rhoN->DrawNormalized("HIST");
   //h_ptsig_response_SK->DrawNormalized("SAME HIST");
  // h_rhoBelow->DrawNormalized("SAME HIST");
  // h_ptsig_response_median->DrawNormalized("SAME HIST");
  // h_ptsig_response_rhosk->DrawNormalized("SAME HIST");

  auto leg4 = new TLegend(0.6,0.6,0.75,0.75);
  leg4->AddEntry(h_ptsig_response_rhoN, "#rho-correction,p_{T}^{cut}=3#mu", "l");
//  leg4->AddEntry(h_ptsig_response_median, "med(#rho)A", "l");
//  leg4->AddEntry(h_ptsig_response_SK, "SoftKiller,a=0.27", "l");

  leg4->SetTextSize(0.05);
  leg4->Draw("SAME");

  TPaveText *pt6 = new TPaveText(0.1,0.4,0.3,1,"blNDC");
  pt6->SetBorderSize(0);
  pt6->SetFillColor(0);
  pt6->SetFillStyle(0);
  pt6->SetTextSize(0.04);
  TText *AText7 = pt6->AddText("<N>=7000, #mu=1.2 GeV/c");
  pt6->Draw();

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
 //c3->SetLogy();
 //c3->SetLeftMargin(50);
 ptsignal_nabove_ptabove->GetXaxis()->SetTitle("N_{const,>}");
 ptsignal_nabove_ptabove->GetZaxis()->SetTitle("<p_{T,<}^{sig}>[GeV/c]");
 ptsignal_nabove_ptabove->GetYaxis()->SetTitle("p_{T,>}[GeV/c]");
 ptsignal_nabove_ptabove->SetTitle("");
 //ptsignal_nabove_ptabove->SetStats(0);
 ptsignal_nabove_ptabove->GetYaxis()->SetLabelSize(0.04);
 ptsignal_nabove_ptabove->GetXaxis()->SetLabelSize(0.04);
 ptsignal_nabove_ptabove->GetYaxis()->SetTitleSize(0.06);
 ptsignal_nabove_ptabove->GetXaxis()->SetTitleSize(0.06);
 ptsignal_nabove_ptabove->GetZaxis()->SetLabelSize(0.04);
 ptsignal_nabove_ptabove->GetZaxis()->SetTitleSize(0.06);

// h_gamma_bkg->SetMarkerColor(azul);
// h_gamma_bkg->SetLineColor(azul);
// h_gamma_bkg->SetMarkerStyle(20);
 ptsignal_nabove_ptabove->SetMarkerColor(rojo);
 ptsignal_nabove_ptabove->SetLineColor(rojo);
 ptsignal_nabove_ptabove->SetMarkerStyle(21);

 ptsignal_nabove_ptabove->Draw("LEGO");
//  h_gamma_bkg->DrawNormalized("HIST");


  TCanvas *c4 = new TCanvas ("c4", "c4", 65, 52, 800, 600);
  // c4->SetLogy();
   h_pt_sig->GetXaxis()->SetTitle("p_{T}^{corr}[GeV/c]");
   h_pt_sig->GetYaxis()->SetTitle("Probability density");
   h_pt_sig->SetTitle("");
   //ptsignal_ptjet->SetStats(0);
   h_pt_sig->GetYaxis()->SetLabelSize(0.04);
   h_pt_sig->GetXaxis()->SetLabelSize(0.04);
   h_pt_sig->GetYaxis()->SetTitleSize(0.06);
   h_pt_sig->GetXaxis()->SetTitleSize(0.06);
   h_pt_sig->SetMarkerColor(azul);
   h_pt_sig->SetLineColor(azul);
   h_pt_sig->SetMarkerStyle(20);
  // h_deltaR_sig->SetMarkerColor(rojo);
  // h_deltaR_sig->SetLineColor(rojo);
  // h_deltaR_sig->SetMarkerStyle(21);
//  h_rhoBelow_nBelow->GetXaxis()->SetTitle("N_{<}");
//h_rhoBelow_nBelow->GetYaxis()->SetTitle("p^{T}_{patch,<}[GeV/c]");
//h_rhoBelow_nBelow->GetYaxis()->SetLabelSize(0.04);
//h_rhoBelow_nBelow->GetXaxis()->SetLabelSize(0.04);
//h_rhoBelow_nBelow->GetYaxis()->SetTitleSize(0.06);
// h_rhoBelow_nBelow->GetXaxis()->SetTitleSize(0.06);
 //h_rhoBelow_nBelow->SetTitle("PYTHIA8");
//h_rhoBelow_nBelow->Draw("COLZ");
//      h_nsignal->DrawNormalized();
  //  h_deltaR_sig->DrawNormalized("HIST");
  //  h_deltaR_bkg->DrawNormalized("SAME HIST");
//   h_ptsig_response_median->DrawNormalized("SAME HIST");
//   h_ptsig_response_rholine->DrawNormalized("SAME HIST");

  h_pt_sig->DrawNormalized();


}

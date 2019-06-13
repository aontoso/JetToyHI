
void rhoEstimator_tree(){

TFile *file =TFile::Open("/Users/albasotoontoso/Work/Jet_substraction/JetToyHI/RHIC_results/JetToyHIResultRhoEstimator_RHIC_pp_16mu_2k.root");
TTree *treef = (TTree *)file->Get("jetTree");
Long64_t nevents = treef->GetEntriesFast();
Long64_t npatches = 0;
Long64_t npatches_sk = 0;

TH1D* h_ptbkg_response_median = new TH1D("p_{T} bkg response", "p_{T} bkg response", 50., -80. , 70);
TH1D* h_ptbkg_response_rhoN = new TH1D("Intersection-#rho(N)", "Intersection-#rho(N)", 50., -8. , 8.);
TH1D* h_ptbkg_response_rholine = new TH1D("p_{T} bkg response THREE", "p_{T} bkg response THREE", 50., -80. , 70);

TH1D* h_ptsig_response_median = new TH1D("p_{T} sig response", "p_{T} bkg response", 50., 0. , 5.);
TH1D* h_ptsig_response_rhoN = new TH1D("p_{T} sig response TWO", "p_{T} sig response TWO", 40., -20. , 20.);
TH1D* h_ptsig_response_rholine = new TH1D("p_{T} sig response THREE", "p_{T} sig response THREE", 50., -80. , 80);
TH1D* h_ptsig_response_rhosk = new TH1D("p_{T} sig response FOUR", "p_{T} sig response FOUR", 100., 0. , 2.);

TH1D* h_ptsig_response_SK = new TH1D("p_{T} sig response SK", "p_{T} sig response SK", 40., -120., 120.);

TH1D* h_gamma_bkg = new TH1D("Gamma bkg", "gamma bkg", 100., 0. , 500.);
TH1D* h_gamma_sig = new TH1D("Gamma sig", "gamma sig", 100., 0. , 300.);
TH1D* h_deltaR_sig = new TH1D("Delta R sig", "Delta R sig", 20., 0. , 0.6);
TH1D* h_deltaR_bkg = new TH1D("Delta R bkg", "Delta R sig", 20., 0. , 0.6);
TH1D* h_nsignal = new TH1D("N sig", "n sig", 100., 0. , 60);
TH1D* h_skPt = new TH1D("sk Pt", "sk Pt", 100., 0. , 6.);
TH1D* h_pt_sig = new TH1D("pT sig", "pT sig", 20., 10. , 600.);

double pTMax_hist = 150.;
double pTMin_hist = 0;
double nbins_pT = 10.;
double nbins_pT_giant = 2.;

double nMax_hist = 100.;
double nMin_hist = 0.;
double nbins_n = 10.;

double nbins_width = 10.;
double widthMax_hist = 1.;

TProfile *ptsignal_ptjet = new TProfile("hprof","Profile of psignal versus pt", 30., 0., 600.);
TProfile *ptsignal_ptabove = new TProfile("Pt-pT","Profile of psignal versus pTabove", nbins_pT, pTMin_hist, pTMax_hist, "s" );

TProfile *ptsignal_ptabove_giant = new TProfile("Pt-pT-giant","Profile of psignal versus pTabove", nbins_pT_giant, pTMax_hist, 350., "s" );

TProfile *ptsignal_nabove = new TProfile("Pt-N","Profile of psignal versus nabove", nbins_n, nMin_hist, nMax_hist, "s" );

TProfile *ptsignal_widthAbove = new TProfile("Pt-Width","Profile of psignal versus Width above", nbins_width, 0., widthMax_hist, "s" );

TProfile2D *ptsignal_nabove_ptabove  = new TProfile2D("Pt-pTN","Profile of psignal versus nabove and ptabove",nbins_n,nMin_hist,nMax_hist,nbins_pT,pTMin_hist,pTMax_hist,0.,200.,"s");

TProfile2D *ptsignal_nabove_ptabove_giant  = new TProfile2D("Pt-pTN-giant","Profile of psignal versus nabove and ptabove",nbins_n,nMin_hist,nMax_hist,nbins_pT_giant,pTMax_hist,600.,0.,200.,"s");

TProfile2D *ptsignal_widthabove_ptabove  = new TProfile2D("Pt-pTWidth","Profile of psignal versus Widthabove and ptabove",nbins_width,0.,widthMax_hist,nbins_pT,pTMin_hist,pTMax_hist,0.,200.,"s");

TProfile2D *ptsignal_widthabove_ptabove_giant = new TProfile2D("Pt-pTWidth giant","Profile of psignal versus Widthabove and ptabove giant",nbins_width,0.,widthMax_hist,nbins_pT_giant,pTMax_hist,600.,0.,200.,"s");


TH1D* h_rhoBelow = new TH1D("Rho Below", "Rho Below", 50., -80. , 80.);

std::vector<double> * rho_median = 0;
std::vector<double> * rho_True = 0;
std::vector<double> * pT_patch = 0;
std::vector<double> * rho_Estimate = 0;
//std::vector<double> * nConstituents_True = 0;
//std::vector<double> * rho_line = 0;
//std::vector<double> * rho_SK = 0;
//std::vector<double> * skPtThreshold = 0;

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

//std::vector<double> * sigJetPt = 0;
//std::vector<double> * sigJetEta = 0;
//std::vector<double> * sigJetPhi = 0;
//std::vector<double> * sigJetM = 0;
//std::vector<double> * sigJetArea = 0;

//std::vector<std::vector<double>> * sigJetConstPt = 0;
//std::vector<std::vector<double>> * sigJetConstEta = 0;
//std::vector<std::vector<double>> * sigJetConstPhi = 0;
//std::vector<std::vector<double>> * sigJetConstM = 0;
//std::vector<std::vector<double>> * sigJetConstId = 0;
//std::vector<std::vector<double>> * sigJetConstDeltaR = 0;

//std::vector<double> * skJetPt = 0;
//std::vector<double> * skJetEta = 0;
//std::vector<double> * skJetPhi = 0;
//std::vector<double> * skJetM = 0;
//std::vector<double> * skJetArea = 0;

//std::vector<std::vector<double>> * skJetConstPt = 0;
//std::vector<std::vector<double>> * skJetConstEta = 0;
//std::vector<std::vector<double>> * skJetConstPhi = 0;
//std::vector<std::vector<double>> * skJetConstM = 0;
//std::vector<std::vector<double>> * skJetConstId = 0;
//std::vector<std::vector<double>> * skJetConstDeltaR = 0;

treef->SetBranchAddress("rho_True",&rho_True);
treef->SetBranchAddress("pT_patch",&pT_patch);
treef->SetBranchAddress("rho_median",&rho_median);
//treef->SetBranchAddress("nConstituents_True", &nConstituents_True);
treef->SetBranchAddress("rho_Estimate",&rho_Estimate);
//treef->SetBranchAddress("rho_line",&rho_line);

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
//treef->SetBranchAddress("rho_SK",&rho_SK);
//treef->SetBranchAddress("skPtThreshold", &skPtThreshold);

//treef->SetBranchAddress("sigJetPt",&sigJetPt);
//treef->SetBranchAddress("sigJetEta",&sigJetEta);
//treef->SetBranchAddress("sigJetPhi",&sigJetPhi);
//treef->SetBranchAddress("sigJetM",&sigJetM);
//treef->SetBranchAddress("sigJetArea",&sigJetArea);
//treef->SetBranchAddress("sigJetConstPt",&sigJetConstPt);
//treef->SetBranchAddress("sigJetConstEta",&sigJetConstEta);
//treef->SetBranchAddress("sigJetConstPhi",&sigJetConstPhi);
//treef->SetBranchAddress("sigJetConstM",&sigJetConstM);
//treef->SetBranchAddress("sigJetConstId",&sigJetConstId);
//treef->SetBranchAddress("sigJetConstDeltaR",&sigJetConstDeltaR);


//treef->SetBranchAddress("skJetPt",&skJetPt);
//treef->SetBranchAddress("skJetEta",&skJetEta);
//treef->SetBranchAddress("skJetPhi",&skJetPhi);
//treef->SetBranchAddress("skJetM",&skJetM);
//treef->SetBranchAddress("skJetArea",&skJetArea);
//treef->SetBranchAddress("skJetConstPt",&skJetConstPt);
//treef->SetBranchAddress("skJetConstEta",&skJetConstEta);
//treef->SetBranchAddress("skJetConstPhi",&skJetConstPhi);
//treef->SetBranchAddress("skJetConstM",&skJetConstM);
//treef->SetBranchAddress("skJetConstId",&skJetConstId);
//treef->SetBranchAddress("skJetConstDeltaR",&skJetConstDeltaR);


Long64_t nconstituents = 0;
int nconstituents2 = 0;
int njets=0;

Long64_t nevents_training = int(nevents/2);
// Antes de hacer el analisis vamos a obtener el TProfile con signal por debajo de un cut en funcion de particulas por encima del cut

double mu = 0.6;
double pTcut_soft = 1.6*mu;
double pTcut_hard = 5*mu;

double pTAbove_cut = 20.;
double Radius = 0.4;

double kappa = 1.;
double beta = 1.;

//ofstream myfile_width;
  // myfile_width.open ("ptsignal_width_pt_nconst_test.txt");
//   myfile_width << "ptSignal  widthAbove\n";

for (Long64_t k=0; k<nevents_training; k++){ // Loop sobre el total de entradas del trees
      treef->GetEntry(k);
      npatches = rho_True->size();

      for (int i=0; i<npatches; i++){
        double pt_sig_true = pT_patch->at(i)-rho_True->at(i);


        if (pt_sig_true>pTAbove_cut){
       nconstituents2 = bkgJetConstPt->at(i).size();

       double ptSignal = 0;
       int nabove = 0;
       double ptAbove = 0;
       double widthAbove = 0;


       for (int j=0; j<nconstituents2; j++){
        //  cout << bkgJetConstId->at(i).at(j) << endl;
             int iD= bkgJetConstId->at(i).at(j);
             double momentum = bkgJetConstPt->at(i).at(j);

             if (momentum < pTcut_soft && iD==0){
             ptSignal+=momentum;}
             if (momentum > pTcut_hard && iD==0) {
               nabove++;
             }
             if (momentum > pTcut_hard && iD==0) ptAbove+=momentum;

           }

        //   if (k==0) cout << nabove << endl;
           //if (ptSignal ==0) cout << ptAbove << endl;
      for (int j=0; j<nconstituents2; j++){
        int iD = bkgJetConstId->at(i).at(j);
      //  cout << iD << endl;
        double momentum = bkgJetConstPt->at(i).at(j);
        double delta_R = bkgJetConstDeltaR->at(i).at(j);
      //  if(iD!=0) cout << "eo " << endl;

        if (momentum > pTcut_hard && iD==0){
          widthAbove += pow(momentum/ptAbove, kappa)*std::pow(delta_R/Radius,beta); //use scalar pT for jet?
        }
      }

     if (ptAbove>pTMax_hist) ptsignal_ptabove_giant->Fill(ptAbove,ptSignal);
     else ptsignal_ptabove->Fill(ptAbove,ptSignal);

      ptsignal_nabove->Fill(nabove,ptSignal);
      if (ptAbove>pTMax_hist) ptsignal_nabove_ptabove_giant->Fill(nabove,ptAbove,ptSignal);
      else ptsignal_nabove_ptabove->Fill(nabove,ptAbove,ptSignal);

     ptsignal_widthAbove->Fill(widthAbove,ptSignal);

      if (ptAbove>pTMax_hist) ptsignal_widthabove_ptabove_giant->Fill(widthAbove,ptAbove,ptSignal);
      else ptsignal_widthabove_ptabove->Fill(widthAbove,ptAbove,ptSignal);

//     myfile_width << ptSignal << " " << widthAbove << " " << ptAbove << " " << nabove << " \n";
    }

  }
}
//myfile_width.close();
//cout << ptsignal_nabove->GetBinEntries(7) << endl;

for (Long64_t k=nevents_training; k<nevents; k++){ // Loop sobre el total de entradas del trees
      treef->GetEntry(k);
      npatches = rho_True->size();
    //  npatches_sk = skJetPt->size();
  //    for (int i = 0; i<npatches_sk; i++) {
  //     int nconst_True = sigJetConstPt->at(i).size();
  //     double scalarPt_True = 0;
  //      for(int j = 0; j<nconst_True; j++){
  //      scalarPt_True+=sigJetConstPt->at(i).at(j);
  //     }
  //     if (scalarPt_True>80) h_pt_sig->Fill(scalarPt_True);
//     }
      // if(scalarPt_True > 130){
      // int nconst_SK = skJetConstPt->at(i).size();
      // double scalarPt_SK = 0;
      // for(int j = 0; j<nconst_SK; j++){
      //  scalarPt_SK+=skJetConstPt->at(i).at(j);
       //}
       //h_ptsig_response_SK->Fill(scalarPt_SK-scalarPt_True);
    //   if (-7.3<scalarPt_SK-scalarPt_True && scalarPt_SK-scalarPt_True<7.3) njets++;
       //}
  //  }

      for (int i=0; i<npatches; i++){
        double pt_sig_true = pT_patch->at(i)-rho_True->at(i);

      //  if(pt_sig_true>pTAbove_cut && pT_patch->at(i)-rho_median->at(i)>pTAbove_cut){
        if(pt_sig_true>pTAbove_cut){
        nconstituents = bkgJetConstPt->at(i).size();
        int nabove_jet = 0;
        double ptabove_jet = 0;
        double widthAbove_jet = 0;

        for (int j=0; j<nconstituents; j++){
        double momentum = bkgJetConstPt->at(i).at(j);
        if (momentum > pTcut_hard) ptabove_jet+=momentum;
        if (momentum > pTcut_hard) nabove_jet++;
        }

        for (int j=0; j<nconstituents; j++){
        double momentum = bkgJetConstPt->at(i).at(j);
        double delta_R = bkgJetConstDeltaR->at(i).at(j);

        if (momentum > pTcut_hard){
          widthAbove_jet += std::pow(momentum/ptabove_jet, kappa)*std::pow(delta_R/Radius,beta); //use scalar pT for jet?
        }
        }


        int ptabove_bin_giant = ptsignal_ptabove_giant->FindFixBin(ptabove_jet);
        int ptabove_bin = ptsignal_ptabove->FindFixBin(ptabove_jet);
        int nabove_bin = ptsignal_nabove->FindFixBin(nabove_jet);

        int widthabove_bin = ptsignal_widthAbove->FindFixBin(widthAbove_jet);

      //  cout << ptsignal_widthAbove->FindLastBinAbove(0) << endl;
        if (widthabove_bin > ptsignal_widthAbove->FindLastBinAbove(0)) widthabove_bin = ptsignal_widthAbove->FindLastBinAbove(0);

        if (ptabove_bin_giant > ptsignal_ptabove_giant->FindLastBinAbove(0)) ptabove_bin_giant = ptsignal_ptabove_giant->FindLastBinAbove(0);

      //  cout << ptsignal_widthAbove->GetBinContent(12) << endl;
  //    cout << ptabove_jet << endl;
         int ptn_bin = ptsignal_nabove_ptabove->FindFixBin(nabove_jet,ptabove_jet);
        //if(ptabove_jet>350) ptabove_jet =300;
        int ptn_bin_giant = ptsignal_nabove_ptabove_giant->FindFixBin(nabove_jet,ptabove_jet);


        int ptwidth_bin = ptsignal_widthabove_ptabove->FindFixBin(widthAbove_jet, ptabove_jet);

        int ptwidth_bin_giant =  ptsignal_widthabove_ptabove_giant->FindFixBin(widthAbove_jet, ptabove_jet);
      //  cout << ptabove_giant << endl;

////////////////////////////////////////////////////////////////////////////
// Corrections/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// pT ///
        double pt_first_term = ptsignal_ptabove->GetBinContent(ptabove_bin);
        double pt_second_term = ptsignal_ptabove_giant->GetBinContent(ptabove_bin_giant);

        if (pt_first_term == 0 && ptabove_bin < nbins_pT) pt_first_term=(ptsignal_ptabove->GetBinContent(ptabove_bin+1)+ptsignal_ptabove->GetBinContent(ptabove_bin-1))/2;


        if (pt_second_term == 0 && ptabove_bin_giant > nbins_pT) pt_second_term=(ptsignal_ptabove_giant->GetBinContent(ptabove_bin_giant+1)+ptsignal_ptabove_giant->GetBinContent(ptabove_bin_giant-1))/2;

    //    cout << ptsignal_ptabove_giant->GetBinContent(3) << endl;

        double pT_correction = pt_first_term+pt_second_term;

// Nparticles ////

        double pT_correction_one = ptsignal_nabove->GetBinContent(nabove_bin);

        if (pT_correction_one == 0 && nabove_bin < nbins_n)
        pT_correction_one=(ptsignal_nabove->GetBinContent(nabove_bin+1)+ptsignal_nabove->GetBinContent(nabove_bin-1))/2;

// Nparticles and pT ///

        double npt_first_term = ptsignal_nabove_ptabove->GetBinContent(ptn_bin);
        double npt_second_term = ptsignal_nabove_ptabove_giant->GetBinContent(ptn_bin_giant);


       if (npt_first_term == 0) npt_first_term=(ptsignal_nabove_ptabove->GetBinContent(ptn_bin+1)+ptsignal_nabove_ptabove->GetBinContent(ptn_bin-1))/2;

        if (npt_second_term == 0) npt_second_term=(ptsignal_nabove_ptabove_giant->GetBinContent(ptn_bin_giant+1)+ptsignal_nabove_ptabove_giant->GetBinContent(ptn_bin_giant-1))/2;

        double pT_correction_two = npt_first_term + npt_second_term;

// Width ///

        double pT_correction_three = ptsignal_widthAbove->GetBinContent(widthabove_bin);

        if (pT_correction_three == 0 && widthabove_bin < nbins_width)
        pT_correction_three=(ptsignal_widthAbove->GetBinContent(widthabove_bin+1)+ptsignal_widthAbove->GetBinContent(widthabove_bin-1))/2;

// Width and pT ////

        double wpt_first_term = ptsignal_widthabove_ptabove->GetBinContent(ptwidth_bin);

        double wpt_second_term = ptsignal_widthabove_ptabove_giant->GetBinContent(ptwidth_bin_giant);

        if (wpt_first_term == 0) wpt_first_term=(ptsignal_widthabove_ptabove->GetBinContent(ptwidth_bin+1)+ptsignal_widthabove_ptabove->GetBinContent(ptwidth_bin-1))/2;

         if (wpt_second_term == 0) wpt_second_term=(ptsignal_nabove_ptabove_giant->GetBinContent(ptwidth_bin_giant+1)+ptsignal_nabove_ptabove_giant->GetBinContent(ptwidth_bin_giant-1))/2;

        double pT_correction_four = wpt_first_term + wpt_second_term;

/// Check /////

       if (pT_correction_three == 0) cout << "w" << ptabove_jet << "pt: " << nabove_jet << "Bin: " << npt_second_term  <<  "Bin Giant: " << npt_second_term << endl;

/////////////////// Reconstructed pT ////////////////////////////////////

         double pt_signal_reco_estimate = pT_patch->at(i)-rho_Estimate->at(i);
    //      double pt_signal_reco_estimate = rho_Estimate->at(i);
          double pt_signal_reco_median = pT_patch->at(i)-rho_median->at(i);
      //    double pt_signal_reco_line = pT_patch->at(i)-rho_line->at(i);
      double pt_signal_reco_line = 0;
          //cout << pT_correction_two << endl;
        //  if (pt_signal_reco_estimate-pt_sig_true ==0) cout << "eo" << endl;
    //    cout << pt_signal_reco_estimate << endl;

         h_ptsig_response_median->Fill(rho_Estimate->at(i));
         h_ptsig_response_rhoN->Fill(pt_signal_reco_estimate-pt_sig_true);
         h_ptsig_response_rholine->Fill(pt_signal_reco_line-pt_sig_true);
        // h_ptsig_response_rhosk->Fill(pt_signal_reco_sk-pt_sig_true);
         nconstituents = bkgJetConstPt->at(i).size();

        int nbelow2 = 0;
        double ptSignal2 = 0;
//        int nsignal = 0;
      for (int j=0; j<nconstituents; j++){
             int iD = bkgJetConstId->at(i).at(j);
             double momentum = bkgJetConstPt->at(i).at(j);
             if (momentum < pTcut_soft && iD==0){
             ptSignal2+=momentum;
            nbelow2++;
           }
           }
        //   cout << ptSignal2 << endl;
          h_gamma_sig->Fill(ptSignal2);

       }


        }

    }
  // Truquis


//double nsig = 3.;
//double nsig2 = 0.5;
 //TF1 *fg1 = NULL;
  //int maxIter = 20;

//  double meanold = h_ptsig_response_rhoN->GetMean();
//  double mean = h_ptsig_response_rhoN->GetMean();
//  double sigma = h_ptsig_response_rhoN->GetRMS();
//  fg1 = new TF1(h_ptsig_response_rhoN->GetName(),"gaus",-nsig*h_ptsig_response_rhoN->GetRMS(),nsig2*h_ptsig_response_rhoN->GetRMS());
//  for(int i = 0; i<maxIter; ++i) {
//    h_ptsig_response_rhoN->Fit(fg1,"R0","",mean-nsig*sigma,mean+nsig2*sigma);
//    meanold = mean;
//    sigma = fg1->GetParameter(2);
//    mean  = fg1->GetParameter(1);
    //cout << nsig2*h_ptsig_response_rhoN->GetRMS() << endl;
//    if(fabs(mean-meanold)<0.001) cout << "mean" << endl;
//  }

//   TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 800, 600);
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
  //    h_ptbkg_response_rhoN->DrawNormalized("HIST");

  //     h_ptbkg_response_median->DrawNormalized("SAME HIST");
  //    h_ptbkg_response_rholine->DrawNormalized("SAME HIST");


      auto leg3 = new TLegend(0.6,0.6,0.75,0.75);
    //  leg3->AddEntry(h_ptbkg_response_rhoN, "YASM+#delta#rho^{sig}_{<}", "l");
      leg3->AddEntry(h_ptbkg_response_rhoN, "med(#rho)A", "l");
  //    leg3->AddEntry(h_ptbkg_response_rholine, "#rho_{corr}(N^{truth})", "l");
      leg3->SetTextSize(0.05);
    //  leg3->Draw("SAME");
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
//h_ptsig_response_rhoN->SetOptStat(000002211);
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
//  leg4->AddEntry(h_ptsig_response_rhoN, "#rho-correction,5K", "l");
//  leg4->AddEntry(h_ptsig_response_median, "med(#rho)A", "l");
  leg4->AddEntry(h_ptsig_response_rhoN, "#rho-correction+(#lambda_{1}^{1},p_{T})_{>}", "l");

  leg4->SetTextSize(0.04);
  leg4->Draw("SAME");

  TPaveText *pt6 = new TPaveText(0.1,0.4,0.3,1,"blNDC");
  pt6->SetBorderSize(0);
  pt6->SetFillColor(0);
  pt6->SetFillStyle(0);
  pt6->SetTextSize(0.04);
  TText *AText7 = pt6->AddText("<N>=7000, #mu=1.2 GeV/c");
  pt6->Draw();

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

  TPaveText *pt800 = new TPaveText(0.1,0.7,0.3,1,"blNDC");
  pt800->SetBorderSize(0);
  pt800->SetFillColor(0);
  pt800->SetFillStyle(0);
  pt800->SetTextSize(0.04);
  TText *AText800 = pt800->AddText("p_{T}^{cut,hard}=5#mu,p_{T}^{cut,soft}=2#mu");
//  pt800->Draw();

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

 ptsignal_widthAbove->Draw();
 pt->Draw();
 pt80->Draw();
 pt800->Draw();
 pt6->Draw();
//  h_gamma_bkg->DrawNormalized("HIST");


  TCanvas *c4 = new TCanvas ("c4", "c4", 65, 52, 800, 600);
//   c4->SetLogx();
   ptsignal_widthAbove->GetXaxis()->SetTitle("Width_{>,hard}");
   ptsignal_widthAbove->GetYaxis()->SetTitle("<p_{T,<}^{sig}>[GeV/c]");
   ptsignal_widthAbove->SetTitle("");
  // ptsignal_widthAbove->SetStats(0);
   ptsignal_widthAbove->GetYaxis()->SetLabelSize(0.04);
   ptsignal_widthAbove->GetXaxis()->SetLabelSize(0.04);
   ptsignal_widthAbove->GetYaxis()->SetTitleSize(0.06);
   ptsignal_widthAbove->GetXaxis()->SetTitleSize(0.06);
   ptsignal_widthAbove->SetMarkerColor(azul);
   ptsignal_widthAbove->SetLineColor(azul);
   ptsignal_widthAbove->SetMarkerStyle(20);
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

  h_gamma_sig->Draw();
  pt->Draw();
  pt80->Draw();
  pt6->Draw();
  pt800->Draw();
  //pt800->Draw();
  //pt6->Draw();

}

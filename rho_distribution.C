//---------------------------------------------------------------
// Description
// Simple macro to analyze the rho-distribution and test the
// Gaussian assumption.
// Author: M. Verweij, Y. Mehtar-Tani, Alba Soto-Ontoso
//---------------------------------------------------------------

Double_t fitf(Double_t *x,Double_t *par) {
    Double_t arg = 0;
    if (par[1]!=0) arg = (x[0] - par[0])/par[1];
    Double_t fitval = TMath::Exp(-0.5*arg*arg)/(pow(2*TMath::Pi()*par[1],0.5));
    return fitval;
 }
 double Median(const TH1D * h1) {

    int n = h1->GetXaxis()->GetNbins();
    std::vector<double>  x(n);
    h1->GetXaxis()->GetCenter( &x[0] );
    const double * y = h1->GetArray();
    // exclude underflow/overflows from bin content array y
    return TMath::Median(n, &x[0], &y[1]);
 }

void rho_distribution() {

//--------- Creat an histogram

TH1D* h_rho = new TH1D("#rho", "Rho dist", 50, 150, 400);

//-------- Read the root file
double rho_value;
double ptd_value;
TFile *f = TFile::Open("rho_dist.root");
TTree *data_tree = (TTree *)f->Get("rho_dist");
data_tree->SetBranchAddress("rho_value", &rho_value);
data_tree->SetBranchAddress("ptd_value", &ptd_value);
Long64_t nentries = data_tree->GetEntries();
//double rho_value_array[nentries];
cout << nentries << endl;
for (Long64_t k = 0; k < nentries; k++){
data_tree->GetEntry(k);
h_rho->Fill(rho_value);
//rho_value_array[k] = rho_value;
}
TF1 *func = new TF1("fit",fitf,160,400,2);
func->SetParameters(50,h_rho->GetMean(),h_rho->GetRMS());
// give the parameters meaningful names
      func->SetParNames ("Mean_value","Sigma");

  //    func->FixParameter(1,12.3201);
      // call TH1::Fit with the name of the TF1 object
      h_rho->Fit("fit");
gStyle->SetOptFit(1111);
TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 700, 600);
Int_t azul;
azul = TColor::GetColor("#034F84");
Int_t rojo;
rojo = TColor::GetColor("#B93A32");
h_rho->GetXaxis()->SetTitle("#rho[GeV]");
h_rho->GetYaxis()->SetTitle("P(#rho)");
h_rho->SetTitle("");
h_rho->GetYaxis()->SetLabelSize(0.04);
h_rho->GetXaxis()->SetLabelSize(0.04);
h_rho->GetYaxis()->SetTitleSize(0.06);
h_rho->GetXaxis()->SetTitleSize(0.06);
h_rho->SetLineWidth(3);
h_rho->SetLineColor(azul);
h_rho->DrawNormalized();
TF1 *fa = new TF1("fa","TMath::Gaus(x, 2.50864e+02, 4.41549e+01)/(pow(2*TMath::Pi()*4.41549e+01,0.5))",150,400);
fa->SetLineWidth(3);
fa->SetLineColor(rojo);
fa->Draw("SAME l");
TF1 *fb = new TF1("fb","TMath::Gaus(x,  251.499, 12.3201)/(pow(2*TMath::Pi()*12.3201,0.5))",150,400);
fb->SetLineWidth(3);
fb->SetLineColor(kBlack);
fb->Draw("SAME l");
auto leg = new TLegend(0.6,0.6,0.8,0.8);
leg->AddEntry(h_rho, "#rho from data", "l");
leg->AddEntry(fa, "Fit", "l");
leg->AddEntry(fb, "Our estimate", "l");
leg->SetTextSize(0.05);
leg->Draw("SAME");
}

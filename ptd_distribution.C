//---------------------------------------------------------------
// Description
// Simple macro to analyze the ptD-distribution
// Author: M. Verweij, Y. Mehtar-Tani, Alba Soto-Ontoso
//---------------------------------------------------------------
Double_t fitf(Double_t *x,Double_t *par) {
    Double_t fitval =TMath::Exp(par[0]+x[0]*par[1]);
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

void ptd_distribution() {

//--------- Creat an histogram

TH1D* h_ptd = new TH1D("p_{TD}", "Rho dist", 50, 0.01, 0.04);

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
for (Long64_t k = 106; k < nentries; k++){
data_tree->GetEntry(k);
h_ptd->Fill(ptd_value);
//rho_value_array[k] = rho_value;
}
TF1 *func = new TF1("fit",fitf,0.01,0.04,2);
func->SetParameters(50,h_ptd->GetMean(),h_ptd->GetRMS());
// give the parameters meaningful names
      func->SetParNames ("Constant","Slope");

  //    func->FixParameter(1,12.3201);
      // call TH1::Fit with the name of the TF1 object
      h_ptd->Fit("fit");
//gStyle->SetOptFit(1111);
TCanvas *c1 = new TCanvas ("c1", "c1", 65, 52, 700, 600);
Int_t azul;
azul = TColor::GetColor("#034F84");
Int_t rojo;
rojo = TColor::GetColor("#B93A32");
h_ptd->GetXaxis()->SetTitle("p_{TD}");
h_ptd->GetYaxis()->SetTitle("P(p_{TD})");
h_ptd->SetTitle("");
h_ptd->GetYaxis()->SetLabelSize(0.04);
h_ptd->GetXaxis()->SetLabelSize(0.04);
h_ptd->GetYaxis()->SetTitleSize(0.06);
h_ptd->GetXaxis()->SetTitleSize(0.06);
h_ptd->SetLineWidth(3);
h_ptd->SetLineColor(azul);
h_ptd->DrawNormalized();
TF1 *fa = new TF1("fa","TMath::Exp(-1.98390e+00-x*8.98919e+01)",0.01,0.04);
fa->SetLineWidth(3);
fa->SetLineColor(rojo);
fa->Draw("SAME l");
TF1 *fb = new TF1("fb","TMath::Gaus(x,  0.0151409, 0.00687818)/(pow(2*TMath::Pi()*0.00687818,0.5))",0.01,0.04);
fb->SetLineWidth(3);
fb->SetLineColor(kBlack);
//fb->Draw("SAME l");
auto leg = new TLegend(0.6,0.6,0.8,0.8);
leg->AddEntry(h_ptd, "p_{TD} from data", "l");
leg->AddEntry(fa, "Exponential fit", "l");
//leg->AddEntry(fb, "Our estimate", "l");
leg->SetTextSize(0.05);
leg->Draw("SAME");
}

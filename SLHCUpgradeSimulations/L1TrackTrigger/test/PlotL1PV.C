#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include <TError.h>


#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text);


// ----------------------------------------------------------------------------------------------------------------
// Main script
// ----------------------------------------------------------------------------------------------------------------


void PlotL1PV(TString name) {

  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;
  
  SetPlotStyle();


  // ----------------------------------------------------------------------------------------------------------------
  TChain* tree = new TChain("L1TrackNtuple/eventTree");
  tree->Add(name+".root");


  // ----------------------------------------------------------------------------------------------------------------
  // leafs & branches

  vector<float>*  vtx_MC;
  vector<float>*  vtx_L1;

  TBranch *b_vtx_MC;
  TBranch *b_vtx_L1;

  vtx_MC = 0;
  vtx_L1 = 0;

  tree->SetBranchAddress("pv_MC", &vtx_MC, &b_vtx_MC);
  tree->SetBranchAddress("pv_L1", &vtx_L1, &b_vtx_L1);
  

  // ----------------------------------------------------------------------------------------------------------------
  // histograms

  TH1F* h_zres = new TH1F("z0_res_jet", "; z_{L1} - z_{MC} [cm]; jets / 0.02 cm", 200,-2,2);

  TH1F* htmp_zres_vsz[50];  //bins from -25 to +25

  TString hname[50] = {"1","2","3","4","5","6","7","8","9","10",
		       "11","12","13","14","15","16","17","18","19","20",
		       "21","22","23","24","25","26","27","28","29","30",
		       "31","32","33","34","35","36","37","38","39","40",
		       "41","42","43","44","45","46","47","48","49","50"};

  for (int ih=0; ih<50; ih++) {
    htmp_zres_vsz[ih] = new TH1F("htmp_zres_vsz_"+hname[ih], "; z_{L1} - z_{MC} [cm]; jets / 0.02 cm", 200,-2,2);
  }

    
  TH1F* h_all_vsz = new TH1F("all_vsz", ";gen z [cm]; Efficiency",50,-25.,25.);
  TH1F* h_pv1_vsz = new TH1F("pv1_vsz", ";gen z [cm]; Efficiency",50,-25.,25.);
  TH1F* h_pv5_vsz = new TH1F("pv5_vsz", ";gen z [cm]; Efficiency",50,-25.,25.);



  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------
  
  // general variables 
  int nevt = tree->GetEntries();
  
  cout << "number of events = " << nevt << endl;
  

  // ----------------------------------------------------------------------------------------------------------------
  // event loop
  for (int i=0; i<nevt; i++) {      

    tree->GetEntry(i,0);

    float zmc = vtx_MC->at(0);
    float zpv = vtx_L1->at(0);

    // vertex resolution
    h_zres->Fill(zpv - zmc);

    for (int ih=0; ih<50; ih++) {
      if (zmc > (-25+ih*1.0) && zmc < (-24+ih*1.0)) htmp_zres_vsz[ih]->Fill(zpv - zmc);
    }

    h_all_vsz->Fill(zmc);
    if (fabs(zpv-zmc) < 0.1) h_pv1_vsz->Fill(zmc);
    if (fabs(zpv-zmc) < 0.5) h_pv5_vsz->Fill(zmc);


  } // end of event loop
  // ----------------------------------------------------------------------------------------------------------------


  //
  // Output file
  // 
  TFile* fout = new TFile("output_L1PV_"+inputFile+"_"+fitter+".root","recreate");

  // ----------------------------------------------------------------------------------------------------------------
  // write/plot histograms

  h_zres->Sumw2();
  h_zres->Scale(1.0/h_zres->GetSum());
  
  
  TF1* fit = new TF1("fit", "gaus", -0.5,0.5);
  h_zres->Fit("fit","R");
    
  TCanvas cz;
  h_zres->SetAxisRange(-1,1,"X");
  h_zres->Draw("hist");
    
  float mean  = fit->GetParameter(1);
  float sigma = fit->GetParameter(2);
  float emean  = fit->GetParError(1);
  float esigma = fit->GetParError(2);
  char mtxt[500];
  char stxt[500];
  char ttt[500];
  sprintf(mtxt,"mean = %.3f #pm %.3f cm",mean,emean);	
  sprintf(stxt,"sigma = %.3f #pm %.3f cm",sigma,esigma);	
  sprintf(ttt,"RMS = %.3f",h_zres->GetRMS());	
  mySmallText(0.6,0.90,1,mtxt);
  mySmallText(0.6,0.85,1,stxt);
  mySmallText(0.6,0.80,1,ttt);
    
  cz.SaveAs("PV_res.png");
  cz.SaveAs("PV_res.eps");


  TH1F* h_zres_vsz = new TH1F("zres_vsz", ";gen z [cm]; Resolution [mm]",50,-25.,25.);

  TF1* fitz[50];
  for (int ih=0; ih<50; ih++) {
    fitz[ih] = new TF1("fitz_"+hname[ih], "gaus", -1,1);
    htmp_zres_vsz[ih]->Fit("fitz_"+hname[ih],"R");
    sigma = fitz[ih]->GetParameter(2);
    esigma = fitz[ih]->GetParError(2);
    h_zres_vsz->SetBinContent(ih+1,sigma*10);
    h_zres_vsz->SetBinError(ih+1,esigma*10);
  }
  h_zres_vsz->SetAxisRange(0,1,"Y");
  h_zres_vsz->Draw();
  h_zres_vsz->Write();
  
  cz.SaveAs("PV_res_vsz.png");
  cz.SaveAs("PV_res_vsz.eps");


  h_pv1_vsz->Sumw2();
  h_pv5_vsz->Sumw2();
  h_all_vsz->Sumw2();

  TH1F* h_eff1 = (TH1F*) h_pv1_vsz->Clone();
  h_eff1->SetName("eff1");
  h_eff1->Divide(h_pv1_vsz, h_all_vsz, 1.0, 1.0, "B");

  TH1F* h_eff5 = (TH1F*) h_pv5_vsz->Clone();
  h_eff5->SetName("eff5");
  h_eff5->Divide(h_pv5_vsz, h_all_vsz, 1.0, 1.0, "B");

  h_eff1->SetAxisRange(0.5,1.01,"Y");

  TCanvas c;
  h_eff1->Draw("ep");
  h_eff1->Write();
  h_eff5->SetMarkerStyle(24);
  h_eff5->SetMarkerColor(2);
  h_eff5->SetLineColor(2);
  h_eff5->Draw("ep,same");
  h_eff5->Write();

  TLegend* l = new TLegend(0.2,0.2,0.5,0.4);
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_eff1, "#Deltaz(L1,PV) < 1mm","lep");
  l->AddEntry(h_eff5, "#Deltaz(L1,PV) < 5mm","lep");
  l->SetTextFont(42);
  l->Draw();	

  c.SaveAs("PV_eff.png");
  c.SaveAs("PV_eff.eps");

  fout->Close();

}


void SetPlotStyle() {

  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42,"x");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}

void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
  Double_t tsize=0.044;
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}

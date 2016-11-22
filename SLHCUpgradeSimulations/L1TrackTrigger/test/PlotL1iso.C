#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TGraph.h"
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


void PlotL1iso(TString name) {
 
  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;
  
  SetPlotStyle();
  
  
  // ----------------------------------------------------------------------------------------------------------------
  // read ntuples
  TChain* tree = new TChain("L1TrackNtuple/eventTree");
  tree->Add(name+".root");
  
  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }
  

  // ----------------------------------------------------------------------------------------------------------------
  // define leafs & branches

  // tracking particles
  vector<float>* tp_pt;
  vector<float>* tp_eta;
  vector<float>* tp_phi;
  vector<float>* tp_z0;
  vector<int>*   tp_pdgid;
  vector<int>*   tp_momid;
  vector<int>*   tp_nmatch;

  
  // *L1 track* properties, for tracking particles matched to a L1 track
  vector<float>* matchtrk_pt;
  vector<float>* matchtrk_eta;
  vector<float>* matchtrk_phi;
  vector<float>* matchtrk_z0;

  vector<float>* reliso;
  vector<float>* absiso;

  TBranch* b_tp_pt;
  TBranch* b_tp_eta;
  TBranch* b_tp_phi;
  TBranch* b_tp_z0;
  TBranch* b_tp_pdgid;
  TBranch* b_tp_momid;
  TBranch* b_tp_nmatch;

  TBranch* b_matchtrk_pt;
  TBranch* b_matchtrk_eta;
  TBranch* b_matchtrk_phi;
  TBranch* b_matchtrk_z0;

  TBranch* b_reliso; 
  TBranch* b_absiso; 

  tp_pt  = 0;
  tp_eta = 0;
  tp_phi = 0;
  tp_z0  = 0;
  tp_pdgid = 0;
  tp_momid = 0;
  tp_nmatch = 0;

  matchtrk_pt  = 0;
  matchtrk_eta = 0;
  matchtrk_phi = 0;
  matchtrk_z0  = 0;

  reliso = 0; 
  absiso = 0; 


  tree->SetBranchAddress("tp_pt",     &tp_pt,     &b_tp_pt);
  tree->SetBranchAddress("tp_eta",    &tp_eta,    &b_tp_eta);
  tree->SetBranchAddress("tp_phi",    &tp_phi,    &b_tp_phi);
  tree->SetBranchAddress("tp_z0",     &tp_z0,     &b_tp_z0);
  tree->SetBranchAddress("tp_pdgid",  &tp_pdgid,  &b_tp_pdgid);
  tree->SetBranchAddress("tp_momid",  &tp_momid,  &b_tp_momid);
  tree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);

  tree->SetBranchAddress("matchtrk_pt",    &matchtrk_pt,    &b_matchtrk_pt);
  tree->SetBranchAddress("matchtrk_eta",   &matchtrk_eta,   &b_matchtrk_eta);
  tree->SetBranchAddress("matchtrk_phi",   &matchtrk_phi,   &b_matchtrk_phi);
  tree->SetBranchAddress("matchtrk_z0",    &matchtrk_z0,    &b_matchtrk_z0);
  
  tree->SetBranchAddress("reliso", &reliso, &b_reliso);
  tree->SetBranchAddress("absiso", &absiso, &b_absiso);


  // ----------------------------------------------------------------------------------------------------------------
  // histograms
  // ----------------------------------------------------------------------------------------------------------------

  TH1F* h_matchtrk_pt_prompt     = new TH1F("matchtrk_pt_prompt",    ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_matchtrk_pt_nonprompt  = new TH1F("matchtrk_pt_nonprompt", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_matchtrk_eta_prompt    = new TH1F("matchtrk_eta_prompt",   ";Tracking particle #eta; Tracking particles / 0.1",             10, -2.5,   2.5);
  TH1F* h_matchtrk_eta_nonprompt = new TH1F("matchtrk_eta_nonprompt",";Tracking particle #eta; Tracking particles / 0.1",             10, -2.5,   2.5);

  TH1F* h_isotrk_pt_prompt[3];
  TH1F* h_isotrk_pt_nonprompt[3];
  TH1F* h_isotrk_eta_prompt[3];
  TH1F* h_isotrk_eta_nonprompt[3];
  TString nisocut[3] = {"1","15","3"};
  float isocut[3] = {0.1,0.15,0.2};

  for (int ih=0; ih<3; ih++) {
    h_isotrk_pt_prompt[ih]     = new TH1F("isotrk"+nisocut[ih]+"_pt_prompt",     ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100,  0,   100.0);
    h_isotrk_pt_nonprompt[ih]  = new TH1F("isotrk"+nisocut[ih]+"_pt_nonprompt",  ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100,  0,   100.0);
    h_isotrk_eta_prompt[ih]    = new TH1F("isotrk"+nisocut[ih]+"_eta_prompt",    ";Tracking particle #eta; Tracking particles / 1.0 GeV",         10, -2.5,   2.5);
    h_isotrk_eta_nonprompt[ih] = new TH1F("isotrk"+nisocut[ih]+"_eta_nonprompt", ";Tracking particle #eta; Tracking particles / 1.0 GeV",         10, -2.5,   2.5);
  }

  // for ROC curves
  float passIso_prompt[100] = {0};
  float passIso_nonprompt[100] = {0};
  float all_prompt[100] = {0};
  float all_nonprompt[100] = {0};


  
  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------
  
  int nevt = tree->GetEntries();
  cout << "number of events = " << nevt << endl;


  // ----------------------------------------------------------------------------------------------------------------
  // event loop
  for (int i=0; i<nevt; i++) {

    tree->GetEntry(i,0);
    
    // ----------------------------------------------------------------------------------------------------------------
    // tracking particle loop
    for (int it=0; it<(int)tp_pt->size(); it++) {
      
      if (tp_pt->at(it) < 20.0) continue;
      if (fabs(tp_eta->at(it)) > 2.4) continue;
      
      bool prompt = false;
      if (abs(tp_pdgid->at(it)) == 13 && abs(tp_momid->at(it)) == 13) prompt = true;
      else if (abs(tp_pdgid->at(it)) != 13) continue; 
      
      // was the tracking particle matched to a L1 track?
      if (tp_nmatch->at(it) < 1) continue;

      // was the isolation variable calculated for this track?
      if (reliso->at(it) < -990) continue;

      
      // fill matched track histograms
      if (prompt) {
	h_matchtrk_pt_prompt->Fill(tp_pt->at(it));
	h_matchtrk_eta_prompt->Fill(tp_eta->at(it));
      }
      else {
	h_matchtrk_pt_nonprompt->Fill(tp_pt->at(it));
	h_matchtrk_eta_nonprompt->Fill(tp_eta->at(it));
      }

      for (int ih=0; ih<3; ih++) {
	if (reliso->at(it) < isocut[ih]) {
	  if (prompt) {
	    h_isotrk_pt_prompt[ih]->Fill(tp_pt->at(it));
	    h_isotrk_eta_prompt[ih]->Fill(tp_eta->at(it));
	  }
	  else {
	    h_isotrk_pt_nonprompt[ih]->Fill(tp_pt->at(it));
	    h_isotrk_eta_nonprompt[ih]->Fill(tp_eta->at(it));
	  }
	}
      }

      for (int ih=0; ih<100; ih++) {
	if (prompt) all_prompt[ih]++;
	else all_nonprompt[ih]++;
	
	if (reliso->at(it) < 0.01*ih) {
	  if (prompt) passIso_prompt[ih]++;
	  else passIso_nonprompt[ih]++;
	}
      }

    } // end of TP loop
    
  } // end of event loop
  // ----------------------------------------------------------------------------------------------------------------
  


  // ----------------------------------------------------------------------------------------------------------------
  // ROC curves
  // ----------------------------------------------------------------------------------------------------------------

  float eff_prompt[3] = {0};
  float eff_nonprompt[3] = {0};
  int ieff = 0;

  for (int ih=0; ih<100; ih++) {
    if (all_prompt[ih] > 0) passIso_prompt[ih] = passIso_prompt[ih]/all_prompt[ih];
    else passIso_prompt[ih] = 0;
    if (all_nonprompt[ih] > 0) passIso_nonprompt[ih] = passIso_nonprompt[ih]/all_nonprompt[ih];
    else passIso_nonprompt[ih] = 0;

    if ((0.01*ih == 0.1 || 0.01*ih == 0.15 || 0.01*ih == 0.2) && ieff<3) {
      eff_prompt[ieff] = passIso_prompt[ih];
      eff_nonprompt[ieff] = passIso_nonprompt[ih];
      ieff++;
    }

  }

  //
  // Output file
  // 
  TFile* fout = new TFile("output_L1Iso_"+inputFile+"_"+fitter+".root","recreate");

  TGraph* g_eff = new TGraph(100,passIso_nonprompt,passIso_prompt);
  TH2F* h_dummy = new TH2F("dummy", "; efficiency (non-prompt #mu); efficiency (prompt #mu)",90,0.1,1.0,10,0.9,1.0);
  g_eff->SetMarkerStyle(8);
  TGraph* g_eff2 = new TGraph(3,eff_nonprompt,eff_prompt);
  g_eff2->SetMarkerStyle(22);
  g_eff2->SetMarkerColor(2);

  TCanvas cc;
  //gPad->SetGridx();  
  //gPad->SetGridy();
  h_dummy->Draw();
  g_eff->Draw("same,ep");
  g_eff2->Draw("same,ep");

  g_eff->SetName("ROC");
  g_eff2->SetName("ROC_Special");
  g_eff->Write();
  g_eff2->Write();

  mySmallText(0.4,0.35,1,"Relative track isolation efficiency for");
  mySmallText(0.4,0.30,1,"different isolation cuts (reliso < 0.0-1.0)");
  mySmallText(0.4,0.25,2,"Cut values = 0.1, 0.15, 0.2");

  cc.SaveAs("isoeff_promptVSnonprompt.png");
  cc.SaveAs("isoeff_promptVSnonprompt.eps");

  fout->Close();


  // ----------------------------------------------------------------------------------------------------------------
  // efficiency plots  
  // ----------------------------------------------------------------------------------------------------------------

  h_matchtrk_pt_prompt->Sumw2();
  h_matchtrk_pt_nonprompt->Sumw2();
  h_matchtrk_eta_prompt->Sumw2();
  h_matchtrk_eta_nonprompt->Sumw2();

  for (int ih=0; ih<3; ih++) {
    h_isotrk_pt_prompt[ih]->Sumw2();
    h_isotrk_pt_nonprompt[ih]->Sumw2();
    h_isotrk_eta_prompt[ih]->Sumw2();
    h_isotrk_eta_nonprompt[ih]->Sumw2();
  }

  TH1F* h_effiso_pt_prompt[3];
  TH1F* h_effiso_pt_nonprompt[3];
  TH1F* h_effiso_eta_prompt[3];
  TH1F* h_effiso_eta_nonprompt[3];

  for (int ih=0; ih<3; ih++) {
    h_effiso_pt_prompt[ih] = (TH1F*) h_isotrk_pt_prompt[ih]->Clone();
    h_effiso_pt_prompt[ih]->SetName("eff_iso"+nisocut[ih]+"_pt_prompt");
    h_effiso_pt_prompt[ih]->GetYaxis()->SetTitle("Efficiency");
    h_effiso_pt_prompt[ih]->Divide(h_isotrk_pt_prompt[ih], h_matchtrk_pt_prompt, 1.0, 1.0, "B");

    h_effiso_eta_prompt[ih] = (TH1F*) h_isotrk_eta_prompt[ih]->Clone();
    h_effiso_eta_prompt[ih]->SetName("eff_iso"+nisocut[ih]+"_eta_prompt");
    h_effiso_eta_prompt[ih]->GetYaxis()->SetTitle("Efficiency");
    h_effiso_eta_prompt[ih]->Divide(h_isotrk_eta_prompt[ih], h_matchtrk_eta_prompt, 1.0, 1.0, "B");

    h_effiso_pt_nonprompt[ih] = (TH1F*) h_isotrk_pt_nonprompt[ih]->Clone();
    h_effiso_pt_nonprompt[ih]->SetName("eff_iso"+nisocut[ih]+"_pt_nonprompt");
    h_effiso_pt_nonprompt[ih]->GetYaxis()->SetTitle("Efficiency");
    h_effiso_pt_nonprompt[ih]->Divide(h_isotrk_pt_nonprompt[ih], h_matchtrk_pt_nonprompt, 1.0, 1.0, "B");

    h_effiso_eta_nonprompt[ih] = (TH1F*) h_isotrk_eta_nonprompt[ih]->Clone();
    h_effiso_eta_nonprompt[ih]->SetName("eff_iso"+nisocut[ih]+"_eta_nonprompt");
    h_effiso_eta_nonprompt[ih]->GetYaxis()->SetTitle("Efficiency");
    h_effiso_eta_nonprompt[ih]->Divide(h_isotrk_eta_nonprompt[ih], h_matchtrk_eta_nonprompt, 1.0, 1.0, "B");

    h_effiso_pt_prompt[ih]->SetAxisRange(0.0,1.01,"Y");
    h_effiso_pt_nonprompt[ih]->SetAxisRange(0.0,1.01,"Y");
    h_effiso_eta_prompt[ih]->SetAxisRange(0.0,1.01,"Y");
    h_effiso_eta_nonprompt[ih]->SetAxisRange(0.0,1.01,"Y");
  }


  TCanvas c;

  gPad->SetGridx();
  gPad->SetGridy();


  // draw and save plots
  h_effiso_pt_prompt[0]->Draw("ep");
  h_effiso_pt_prompt[1]->SetMarkerStyle(22);
  h_effiso_pt_prompt[1]->SetMarkerColor(2);
  h_effiso_pt_prompt[1]->SetLineColor(2);
  h_effiso_pt_prompt[1]->Draw("same,ep");
  h_effiso_pt_prompt[2]->SetMarkerStyle(24);
  h_effiso_pt_prompt[2]->SetMarkerColor(4);
  h_effiso_pt_prompt[2]->SetLineColor(4);
  h_effiso_pt_prompt[2]->Draw("same,ep");

  TLegend* l = new TLegend(0.2,0.2,0.5,0.4);
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_effiso_pt_prompt[0], "reliso < 0.10","lep");
  l->AddEntry(h_effiso_pt_prompt[1],"reliso < 0.15","lep");
  l->AddEntry(h_effiso_pt_prompt[2], "reliso < 0.20","lep");
  l->SetTextFont(42);
  l->Draw();	

  c.SaveAs("isoeff_pt_prompt.eps");
  c.SaveAs("isoeff_pt_prompt.png");


  h_effiso_pt_nonprompt[0]->Draw("ep");
  h_effiso_pt_nonprompt[1]->SetMarkerStyle(22);
  h_effiso_pt_nonprompt[1]->SetMarkerColor(2);
  h_effiso_pt_nonprompt[1]->SetLineColor(2);
  h_effiso_pt_nonprompt[1]->Draw("same,ep");
  h_effiso_pt_nonprompt[2]->SetMarkerStyle(24);
  h_effiso_pt_nonprompt[2]->SetMarkerColor(4);
  h_effiso_pt_nonprompt[2]->SetLineColor(4);
  h_effiso_pt_nonprompt[2]->Draw("same,ep");
  l->Draw();	
  c.SaveAs("isoeff_pt_nonprompt.eps");
  c.SaveAs("isoeff_pt_nonprompt.png");


  h_effiso_eta_prompt[0]->Draw("ep");
  h_effiso_eta_prompt[1]->SetMarkerStyle(22);
  h_effiso_eta_prompt[1]->SetMarkerColor(2);
  h_effiso_eta_prompt[1]->SetLineColor(2);
  h_effiso_eta_prompt[1]->Draw("same,ep");
  h_effiso_eta_prompt[2]->SetMarkerStyle(24);
  h_effiso_eta_prompt[2]->SetMarkerColor(4);
  h_effiso_eta_prompt[2]->SetLineColor(4);
  h_effiso_eta_prompt[2]->Draw("same,ep");
  l->Draw();	
  c.SaveAs("isoeff_eta_prompt.eps");
  c.SaveAs("isoeff_eta_prompt.png");

  h_effiso_eta_nonprompt[0]->Draw("ep");
  h_effiso_eta_nonprompt[1]->SetMarkerStyle(22);
  h_effiso_eta_nonprompt[1]->SetMarkerColor(2);
  h_effiso_eta_nonprompt[1]->SetLineColor(2);
  h_effiso_eta_nonprompt[1]->Draw("same,ep");
  h_effiso_eta_nonprompt[2]->SetMarkerStyle(24);
  h_effiso_eta_nonprompt[2]->SetMarkerColor(4);
  h_effiso_eta_nonprompt[2]->SetLineColor(4);
  h_effiso_eta_nonprompt[2]->Draw("same,ep");
  l->Draw();	
  c.SaveAs("isoeff_eta_nonprompt.eps");
  c.SaveAs("isoeff_eta_nonprompt.png");


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


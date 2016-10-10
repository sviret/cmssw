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

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text); 


// ----------------------------------------------------------------------------------------------------------------
// Main script
void RunFake(TString type) {

  SetPlotStyle();


  // ----------------------------------------------------------------------------------------------------------------
  // read ntuples
  TChain* tree = new TChain("L1TrackNtuple/eventTree");
  tree->Add(type+".root");
  

  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }


  // ----------------------------------------------------------------------------------------------------------------
  // define leafs & branches

  vector<float>* tp_pt;
  vector<float>* tp_eta;
  vector<float>* tp_phi;
  vector<float>* tp_dxy;
  vector<float>* tp_z0;
  vector<float>* tp_d0;
  vector<int>*   tp_pdgid;
  vector<int>*   tp_nmatch;
  vector<int>*   tp_nstub;
  vector<int>*   tp_nstublayer;
  vector<int>*   tp_eventid;

  vector<float>* trk_pt;
  vector<float>* trk_eta;
  vector<float>* trk_phi;
  vector<float>* trk_chi2;
  vector<float>* trk_consistency;
  vector<int>*   trk_nstub;
  vector<int>*   trk_fake;

  vector<int>*   trk_matchtp_pdgid;
  vector<float>* trk_matchtp_pt;
  vector<float>* trk_matchtp_eta;
  vector<float>* trk_matchtp_phi;
  vector<float>* trk_matchtp_z0;
  vector<float>* trk_matchtp_dxy;

  TBranch* b_tp_pt;
  TBranch* b_tp_eta;
  TBranch* b_tp_phi;
  TBranch* b_tp_dxy;
  TBranch* b_tp_z0;
  TBranch* b_tp_d0;
  TBranch* b_tp_pdgid;
  TBranch* b_tp_nmatch;
  TBranch* b_tp_nstub;
  TBranch* b_tp_nstublayer;
  TBranch* b_tp_eventid;

  TBranch* b_trk_pt; 
  TBranch* b_trk_eta; 
  TBranch* b_trk_phi; 
  TBranch* b_trk_chi2; 
  TBranch* b_trk_consistency; 
  TBranch* b_trk_nstub; 
  TBranch* b_trk_fake;
  
  TBranch* b_trk_matchtp_pdgid;
  TBranch* b_trk_matchtp_pt;
  TBranch* b_trk_matchtp_eta;
  TBranch* b_trk_matchtp_phi;
  TBranch* b_trk_matchtp_z0;
  TBranch* b_trk_matchtp_dxy;


  tp_pt  = 0;
  tp_eta = 0;
  tp_phi = 0;
  tp_dxy = 0;
  tp_z0  = 0;
  tp_d0  = 0;
  tp_pdgid = 0;
  tp_nmatch = 0;
  tp_nstub = 0;
  tp_nstublayer = 0;
  tp_eventid = 0;

  trk_pt = 0; 
  trk_eta = 0; 
  trk_phi = 0; 
  trk_chi2 = 0; 
  trk_consistency = 0; 
  trk_nstub = 0; 
  trk_fake = 0; 

  trk_matchtp_pdgid = 0;
  trk_matchtp_pt = 0;
  trk_matchtp_eta = 0;
  trk_matchtp_phi = 0;
  trk_matchtp_z0 = 0;
  trk_matchtp_dxy = 0;


  tree->SetBranchAddress("tp_pt",     &tp_pt,     &b_tp_pt);
  tree->SetBranchAddress("tp_eta",    &tp_eta,    &b_tp_eta);
  tree->SetBranchAddress("tp_phi",    &tp_phi,    &b_tp_phi);
  tree->SetBranchAddress("tp_dxy",    &tp_dxy,    &b_tp_dxy);
  tree->SetBranchAddress("tp_z0",     &tp_z0,     &b_tp_z0);
  tree->SetBranchAddress("tp_d0",     &tp_d0,     &b_tp_d0);
  tree->SetBranchAddress("tp_pdgid",  &tp_pdgid,  &b_tp_pdgid);
  tree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
  tree->SetBranchAddress("tp_nstub",  &tp_nstub,  &b_tp_nstub);
  tree->SetBranchAddress("tp_nstublayer", &tp_nstublayer, &b_tp_nstublayer);
  tree->SetBranchAddress("tp_eventid",&tp_eventid,&b_tp_eventid);

  tree->SetBranchAddress("trk_pt",   &trk_pt,   &b_trk_pt);
  tree->SetBranchAddress("trk_eta",  &trk_eta,  &b_trk_eta);
  tree->SetBranchAddress("trk_phi",  &trk_phi,  &b_trk_phi);
  tree->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
  tree->SetBranchAddress("trk_consistency", &trk_consistency, &b_trk_consistency);
  tree->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
  tree->SetBranchAddress("trk_fake",  &trk_fake,  &b_trk_fake);
  
  tree->SetBranchAddress("trk_matchtp_pdgid", &trk_matchtp_pdgid, &b_trk_matchtp_pdgid);
  tree->SetBranchAddress("trk_matchtp_pt",    &trk_matchtp_pt,    &b_trk_matchtp_pt);
  tree->SetBranchAddress("trk_matchtp_eta",   &trk_matchtp_eta,   &b_trk_matchtp_eta);
  tree->SetBranchAddress("trk_matchtp_phi",   &trk_matchtp_phi,   &b_trk_matchtp_phi);
  tree->SetBranchAddress("trk_matchtp_z0",    &trk_matchtp_z0,    &b_trk_matchtp_z0);
  tree->SetBranchAddress("trk_matchtp_dxy",   &trk_matchtp_dxy,   &b_trk_matchtp_dxy);


  // ----------------------------------------------------------------------------------------------------------------
  // histograms
  // ----------------------------------------------------------------------------------------------------------------

  TH1F* h_trk_pt = new TH1F("trk_pt", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trk_4stub_pt = new TH1F("trk_4stub_pt", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trk_4stubchi2_pt = new TH1F("trk_4stubchi2_pt", ";Track p_{T} [GeV]; ",40,0,20);

  TH1F* h_trk_pt3 = new TH1F("trk_pt3", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trk_4stub_pt3 = new TH1F("trk_4stub_pt3", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trk_4stubchi2_pt3 = new TH1F("trk_4stubchi2_pt3", ";Track p_{T} [GeV]; ",40,0,20);

  TH1F* h_trkGoodP_pt = new TH1F("trkGoodP_pt", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkGoodNP_pt = new TH1F("trkGoodNP_pt", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkFake_pt = new TH1F("trkFake_pt", ";Track p_{T} [GeV]; ",40,0,20);

  TH1F* h_trk_qc_pt = new TH1F("trk_qc_pt", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkGoodP_qc_pt = new TH1F("trkGoodP_qc_pt", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkGoodNP_qc_pt = new TH1F("trkGoodNP_qc_pt", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkFake_qc_pt = new TH1F("trkFake_qc_pt", ";Track p_{T} [GeV]; ",40,0,20);

  TH1F* h_trk_pt_C = new TH1F("trk_pt_C", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkGoodP_pt_C = new TH1F("trkGoodP_pt_C", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkGoodNP_pt_C = new TH1F("trkGoodNP_pt_C", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkFake_pt_C = new TH1F("trkFake_pt_C", ";Track p_{T} [GeV]; ",40,0,20);

  TH1F* h_trk_qc_pt_C = new TH1F("trk_qc_pt_C", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkGoodP_qc_pt_C = new TH1F("trkGoodP_qc_pt_C", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkGoodNP_qc_pt_C = new TH1F("trkGoodNP_qc_pt_C", ";Track p_{T} [GeV]; ",40,0,20);
  TH1F* h_trkFake_qc_pt_C = new TH1F("trkFake_qc_pt_C", ";Track p_{T} [GeV]; ",40,0,20);

  TH1F* h_trk_eta_pt2 = new TH1F("trk_eta_pt2", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodP_eta_pt2 = new TH1F("trkGoodP_eta_pt2", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodNP_eta_pt2 = new TH1F("trkGoodNP_eta_pt2", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkFake_eta_pt2 = new TH1F("trkFake_eta_pt2", ";Track #eta; ",25,-2.5,2.5);

  TH1F* h_trk_qc_eta_pt2 = new TH1F("trk_qc_eta_pt2", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodP_qc_eta_pt2 = new TH1F("trkGoodP_qc_eta_pt2", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodNP_qc_eta_pt2 = new TH1F("trkGoodNP_qc_eta_pt2", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkFake_qc_eta_pt2 = new TH1F("trkFake_qc_eta_pt2", ";Track #eta; ",25,-2.5,2.5);

  TH1F* h_trk_eta_pt5 = new TH1F("trk_eta_pt5", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodP_eta_pt5 = new TH1F("trkGoodP_eta_pt5", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodNP_eta_pt5 = new TH1F("trkGoodNP_eta_pt5", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkFake_eta_pt5 = new TH1F("trkFake_eta_pt5", ";Track #eta; ",25,-2.5,2.5);

  TH1F* h_trk_qc_eta_pt5 = new TH1F("trk_qc_eta_pt5", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodP_qc_eta_pt5 = new TH1F("trkGoodP_qc_eta_pt5", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodNP_qc_eta_pt5 = new TH1F("trkGoodNP_qc_eta_pt5", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkFake_qc_eta_pt5 = new TH1F("trkFake_qc_eta_pt5", ";Track #eta; ",25,-2.5,2.5);

  TH1F* h_trk_eta_pt10 = new TH1F("trk_eta_pt10", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodP_eta_pt10 = new TH1F("trkGoodP_eta_pt10", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodNP_eta_pt10 = new TH1F("trkGoodNP_eta_pt10", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkFake_eta_pt10 = new TH1F("trkFake_eta_pt10", ";Track #eta; ",25,-2.5,2.5);

  TH1F* h_trk_qc_eta_pt10 = new TH1F("trk_qc_eta_pt10", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodP_qc_eta_pt10 = new TH1F("trkGoodP_qc_eta_pt10", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkGoodNP_qc_eta_pt10 = new TH1F("trkGoodNP_qc_eta_pt10", ";Track #eta; ",25,-2.5,2.5);
  TH1F* h_trkFake_qc_eta_pt10 = new TH1F("trkFake_qc_eta_pt10", ";Track #eta; ",25,-2.5,2.5);

  TH1F* h_tp_pt = new TH1F("tp_pt", ";TP p_{T} [GeV]; ",40,0,20);
  TH1F* h_tp_pt3 = new TH1F("tp_pt3", ";TP p_{T} [GeV]; ",40,0,20);



  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------
  
  int nevt = tree->GetEntries();
  cout << "number of events = " << nevt << endl;

  int nfake2 = 0;
  int nfake2_qc = 0;
  int nall2 = 0;
  int nall2_qc = 0;

  int nfake5 = 0;
  int nfake5_qc = 0;
  int nall5 = 0;
  int nall5_qc = 0;

  int nfake10 = 0;
  int nfake10_qc = 0;
  int nall10 = 0;
  int nall10_qc = 0;


  // ----------------------------------------------------------------------------------------------------------------
  // event loop

  int ntp = 0;
  int ntp_pt3 = 0;
  int ntp0 = 0;
  int ntp0_pt3 = 0;

  int ntp5 = 0;
  int ntp5_pt3 = 0;
  int ntp50 = 0;
  int ntp50_pt3 = 0;

  for (int i=0; i<nevt; i++) {
    
    tree->GetEntry(i,0);

    for (int itp=0; itp<(int)tp_pt->size(); itp++) {

      if (tp_nstub->at(itp) < 4) continue;
      if (tp_nstublayer->at(itp) < 4) continue;

      if (fabs(tp_dxy->at(itp)) > 1) continue;
      if (tp_pt->at(itp) < 2.0) continue;
      if (fabs(tp_eta->at(itp)) > 2.4) continue;

      ntp++;
      h_tp_pt->Fill(tp_pt->at(itp));
      if (tp_pt->at(itp) > 3.0) {
	ntp_pt3++;
	h_tp_pt3->Fill(tp_pt->at(itp));
      }

      if (tp_nstublayer->at(itp) > 4) {
	ntp5++;
	if (tp_pt->at(itp) > 3.0) ntp5_pt3++;
      }

      if (tp_eventid->at(itp) != 0) continue;
      ntp0++;
      if (tp_pt->at(itp) > 3.0) ntp0_pt3++;

      if (tp_nstublayer->at(itp) > 4) {
	ntp50++;
	if (tp_pt->at(itp) > 3.0) ntp50_pt3++;
      }

    }
        
    // ----------------------------------------------------------------------------------------------------------------
    // track loop
    for (int it=0; it<(int)trk_pt->size(); it++) {
      
      if (false) {
	cout << "trk_matchtp_z0->at(it) = " << trk_matchtp_z0->at(it) << endl;
	cout << "trk_matchtp_dxy->at(it) = " << trk_matchtp_dxy->at(it) << endl;
	cout << "trk_fake->at(it) = " << trk_fake->at(it) << endl;
	cout << "trk_pt->at(it) = " << trk_pt->at(it) << endl;
      }


      // ----------------------------------------------------------------------------------------------------------------
      // all tracks

      if (trk_pt->at(it) > 2.0) {
	nall2++;
	if (trk_fake->at(it)==0) nfake2++;
      }
      if (trk_pt->at(it) > 5.0) {
	nall5++;
	if (trk_fake->at(it)==0) nfake5++;
      }
      if (trk_pt->at(it) > 10.0) {
	nall10++;
	if (trk_fake->at(it)==0) nfake10++;
      }

      if (trk_pt->at(it) > 2.0) {
	h_trk_pt->Fill(trk_pt->at(it));
	if (trk_nstub->at(it)>=4) h_trk_4stub_pt->Fill(trk_pt->at(it));
	if (trk_nstub->at(it)>=4 && trk_chi2->at(it)<100.) h_trk_4stubchi2_pt->Fill(trk_pt->at(it));
      }

      if (trk_pt->at(it)>3.0) {
	h_trk_pt3->Fill(trk_pt->at(it));
	if (trk_nstub->at(it)>=4) h_trk_4stub_pt3->Fill(trk_pt->at(it));
	if (trk_nstub->at(it)>=4 && trk_chi2->at(it)<100.) h_trk_4stubchi2_pt3->Fill(trk_pt->at(it));
      }

      if (trk_fake->at(it)==0) h_trkFake_pt->Fill(trk_pt->at(it));
      else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_pt->Fill(trk_pt->at(it));
      else h_trkGoodNP_pt->Fill(trk_pt->at(it));

      if (fabs(trk_eta->at(it)) < 1.0) {
	h_trk_pt_C->Fill(trk_pt->at(it));
	if (trk_fake->at(it)==0) h_trkFake_pt_C->Fill(trk_pt->at(it));
	else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_pt_C->Fill(trk_pt->at(it));
	else h_trkGoodNP_pt_C->Fill(trk_pt->at(it));
      }

      h_trk_eta_pt2->Fill(trk_eta->at(it));
      if (trk_fake->at(it)==0) h_trkFake_eta_pt2->Fill(trk_eta->at(it));
      else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_eta_pt2->Fill(trk_eta->at(it));
      else h_trkGoodNP_eta_pt2->Fill(trk_eta->at(it));

      if (trk_pt->at(it) > 5.0) {
	h_trk_eta_pt5->Fill(trk_eta->at(it));
	if (trk_fake->at(it)==0) h_trkFake_eta_pt5->Fill(trk_eta->at(it));
	else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_eta_pt5->Fill(trk_eta->at(it));
	else h_trkGoodNP_eta_pt5->Fill(trk_eta->at(it));
      }

      if (trk_pt->at(it) > 10.0) {
	h_trk_eta_pt10->Fill(trk_eta->at(it));
	if (trk_fake->at(it)==0) h_trkFake_eta_pt10->Fill(trk_eta->at(it));
	else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_eta_pt10->Fill(trk_eta->at(it));
	else h_trkGoodNP_eta_pt10->Fill(trk_eta->at(it));
      }


      // ----------------------------------------------------------------------------------------------------------------
      // loose quality cuts
      //if (trk_nstub->at(it) > 3 && trk_chi2->at(it) < 100) {
      if (trk_nstub->at(it) > 3) {
      
	if (trk_pt->at(it) > 2.0) {
	  nall2_qc++;
	  if (trk_fake->at(it)==0) nfake2_qc++;
	}
	if (trk_pt->at(it) > 5.0) {
	  nall5_qc++;
	  if (trk_fake->at(it)==0) nfake5_qc++;
	}
	if (trk_pt->at(it) > 10.0) {
	  nall10_qc++;
	  if (trk_fake->at(it)==0) nfake10_qc++;
	}

	h_trk_qc_pt->Fill(trk_pt->at(it));
	if (trk_fake->at(it)==0) h_trkFake_qc_pt->Fill(trk_pt->at(it));
	else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_qc_pt->Fill(trk_pt->at(it));
	else h_trkGoodNP_qc_pt->Fill(trk_pt->at(it));
	
	if (fabs(trk_eta->at(it)) < 1.0) {
	  h_trk_qc_pt_C->Fill(trk_pt->at(it));
	  if (trk_fake->at(it)==0) h_trkFake_qc_pt_C->Fill(trk_pt->at(it));
	  else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_qc_pt_C->Fill(trk_pt->at(it));
	  else h_trkGoodNP_qc_pt_C->Fill(trk_pt->at(it));
	}

	h_trk_qc_eta_pt2->Fill(trk_eta->at(it));
	if (trk_fake->at(it)==0) h_trkFake_qc_eta_pt2->Fill(trk_eta->at(it));
	else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_qc_eta_pt2->Fill(trk_eta->at(it));
	else h_trkGoodNP_qc_eta_pt2->Fill(trk_eta->at(it));

	if (trk_pt->at(it) > 5.0) {
	  h_trk_qc_eta_pt5->Fill(trk_eta->at(it));
	  if (trk_fake->at(it)==0) h_trkFake_qc_eta_pt5->Fill(trk_eta->at(it));
	  else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_qc_eta_pt5->Fill(trk_eta->at(it));
	  else h_trkGoodNP_qc_eta_pt5->Fill(trk_eta->at(it));
	}

	if (trk_pt->at(it) > 10.0) {
	  h_trk_qc_eta_pt10->Fill(trk_eta->at(it));
	  if (trk_fake->at(it)==0) h_trkFake_qc_eta_pt10->Fill(trk_eta->at(it));
	  else if (fabs(trk_matchtp_z0->at(it)) < 30. && fabs(trk_matchtp_dxy->at(it)) < 1.0) h_trkGoodP_qc_eta_pt10->Fill(trk_eta->at(it));
	  else h_trkGoodNP_qc_eta_pt10->Fill(trk_eta->at(it));
	}
	
      }


    } // end of track loop
    // ----------------------------------------------------------------------------------------------------------------
    
  } // end of event loop
  // ----------------------------------------------------------------------------------------------------------------
  
  cout << "# TP/event = " << (float)ntp/nevt << endl;
  cout << "# TP eventID=0/event = " << (float)ntp0/nevt << endl;
  cout << "# TP/event, pt>3 = " << (float)ntp_pt3/nevt << endl;
  cout << "# TP eventID=0/event, pt>3 = " << (float)ntp0_pt3/nevt << endl;
  cout << "# TP(5stub)/event = " << (float)ntp5/nevt << endl;
  cout << "# TP(5stub) eventID=0/event = " << (float)ntp50/nevt << endl;
  cout << "# TP(5stub)/event, pt>3 = " << (float)ntp5_pt3/nevt << endl;
  cout << "# TP(5stub) eventID=0/event, pt>3 = " << (float)ntp50_pt3/nevt << endl;

  cout << "fake rate > 2 GeV  = " << (float)nfake2/nall2 << endl;
  cout << "fake rate > 5 GeV  = " << (float)nfake5/nall5 << endl;
  cout << "fake rate > 10 GeV = " << (float)nfake10/nall10 << endl;

  cout << "fake rate QC > 2 GeV  = " << (float)nfake2_qc/nall2_qc << endl;
  cout << "fake rate QC > 5 GeV  = " << (float)nfake5_qc/nall5_qc << endl;
  cout << "fake rate QC > 10 GeV = " << (float)nfake10_qc/nall10_qc << endl;



  TFile* fout = new TFile("fakes_"+type+".root","recreate");

  h_trk_pt->Write();
  h_trkGoodP_pt->Write();
  h_trkGoodNP_pt->Write();
  h_trkFake_pt->Write();
  h_trk_qc_pt->Write();
  h_trkGoodP_qc_pt->Write();
  h_trkGoodNP_qc_pt->Write();
  h_trkFake_qc_pt->Write();

  h_trk_pt_C->Write();
  h_trkGoodP_pt_C->Write();
  h_trkGoodNP_pt_C->Write();
  h_trkFake_pt_C->Write();
  h_trk_qc_pt_C->Write();
  h_trkGoodP_qc_pt_C->Write();
  h_trkGoodNP_qc_pt_C->Write();
  h_trkFake_qc_pt_C->Write();

  h_trk_eta_pt2->Write();
  h_trkGoodP_eta_pt2->Write();
  h_trkGoodNP_eta_pt2->Write();
  h_trkFake_eta_pt2->Write();
  h_trk_eta_pt5->Write();
  h_trkGoodP_eta_pt5->Write();
  h_trkGoodNP_eta_pt5->Write();
  h_trkFake_eta_pt5->Write();
  h_trk_eta_pt10->Write();
  h_trkGoodP_eta_pt10->Write();
  h_trkGoodNP_eta_pt10->Write();
  h_trkFake_eta_pt10->Write();

  h_trk_qc_eta_pt2->Write();
  h_trkGoodP_qc_eta_pt2->Write();
  h_trkGoodNP_qc_eta_pt2->Write();
  h_trkFake_qc_eta_pt2->Write();
  h_trk_qc_eta_pt5->Write();
  h_trkGoodP_qc_eta_pt5->Write();
  h_trkGoodNP_qc_eta_pt5->Write();
  h_trkFake_qc_eta_pt5->Write();
  h_trk_qc_eta_pt10->Write();
  h_trkGoodP_qc_eta_pt10->Write();
  h_trkGoodNP_qc_eta_pt10->Write();
  h_trkFake_qc_eta_pt10->Write();


  h_trk_pt->Sumw2();
  h_trkGoodP_pt->Sumw2();
  h_trkGoodNP_pt->Sumw2();
  h_trkFake_pt->Sumw2();
  h_trk_qc_pt->Sumw2();
  h_trkGoodP_qc_pt->Sumw2();
  h_trkGoodNP_qc_pt->Sumw2();
  h_trkFake_qc_pt->Sumw2();

  h_trk_pt_C->Sumw2();
  h_trkGoodP_pt_C->Sumw2();
  h_trkGoodNP_pt_C->Sumw2();
  h_trkFake_pt_C->Sumw2();
  h_trk_qc_pt_C->Sumw2();
  h_trkGoodP_qc_pt_C->Sumw2();
  h_trkGoodNP_qc_pt_C->Sumw2();
  h_trkFake_qc_pt_C->Sumw2();


  h_trk_eta_pt2->Sumw2();
  h_trkGoodP_eta_pt2->Sumw2();
  h_trkGoodNP_eta_pt2->Sumw2();
  h_trkFake_eta_pt2->Sumw2();
  h_trk_eta_pt5->Sumw2();
  h_trkGoodP_eta_pt5->Sumw2();
  h_trkGoodNP_eta_pt5->Sumw2();
  h_trkFake_eta_pt5->Sumw2();
  h_trk_eta_pt10->Sumw2();
  h_trkGoodP_eta_pt10->Sumw2();
  h_trkGoodNP_eta_pt10->Sumw2();
  h_trkFake_eta_pt10->Sumw2();

  h_trk_qc_eta_pt2->Sumw2();
  h_trkGoodP_qc_eta_pt2->Sumw2();
  h_trkGoodNP_qc_eta_pt2->Sumw2();
  h_trkFake_qc_eta_pt2->Sumw2();
  h_trk_qc_eta_pt5->Sumw2();
  h_trkGoodP_qc_eta_pt5->Sumw2();
  h_trkGoodNP_qc_eta_pt5->Sumw2();
  h_trkFake_qc_eta_pt5->Sumw2();
  h_trk_qc_eta_pt10->Sumw2();
  h_trkGoodP_qc_eta_pt10->Sumw2();
  h_trkGoodNP_qc_eta_pt10->Sumw2();
  h_trkFake_qc_eta_pt10->Sumw2();

  
  // total rate vs pt
  h_tp_pt->Scale(1.0/nevt);
  h_tp_pt3->Scale(1.0/nevt);
  h_tp_pt->SetLineColor(4);
  h_tp_pt->SetLineStyle(2);
  h_tp_pt3->SetLineColor(4);
  h_tp_pt3->SetLineStyle(2);
  h_tp_pt->GetXaxis()->SetTitle("Track p_{T} [GeV]");
  h_tp_pt->GetYaxis()->SetTitle("Tracks / event");
  h_tp_pt3->GetXaxis()->SetTitle("Track p_{T} [GeV]");
  h_tp_pt3->GetYaxis()->SetTitle("Tracks / event");

  h_tp_pt->Write();
  h_tp_pt3->Write();

  TCanvas ctot;
  TH1F* h_rate_pt = (TH1F*) h_trk_pt->Clone();
  h_rate_pt->SetName("rate_pt");
  h_rate_pt->Scale(1.0/nevt);
  h_rate_pt->Write();
  h_rate_pt->GetYaxis()->SetTitle("Tracks / event");
  float max = h_tp_pt->GetMaximum();
  if (h_rate_pt->GetMaximum() > max) max = h_rate_pt->GetMaximum();
  h_tp_pt->SetAxisRange(0,max*1.05,"Y");  
  h_tp_pt->Draw();
  h_rate_pt->Draw("hist,same");
  char txt[500];
  sprintf(txt,"average # tracks/event = %.1f",h_rate_pt->GetSum());
  mySmallText(0.5,0.85,1,txt);
  char txt3[500];
  sprintf(txt3,"average # TPs(stubs in #geq 4 layers)/");
  char txt2[500];
  sprintf(txt2,"event = %.1f",h_tp_pt->GetSum());
  mySmallText(0.5,0.79,4,txt3);
  mySmallText(0.5,0.74,4,txt2);
  ctot.SaveAs("FakePlots/tot_vspt_"+type+".png");
  ctot.SaveAs("FakePlots/tot_vspt_"+type+".eps");

  TH1F* h_rate_4stub_pt = (TH1F*) h_trk_4stub_pt->Clone();
  h_rate_4stub_pt->SetName("rate_4stub_pt");
  h_rate_4stub_pt->Scale(1.0/nevt);
  h_rate_4stub_pt->Write();
  h_rate_4stub_pt->GetYaxis()->SetTitle("Tracks / event");
  max = h_tp_pt->GetMaximum();
  if (h_rate_4stub_pt->GetMaximum() > max) max = h_rate_4stub_pt->GetMaximum();
  h_tp_pt->SetAxisRange(0,max*1.05,"Y");  
  h_tp_pt->Draw();
  h_rate_4stub_pt->Draw("same");
  sprintf(txt,"average # tracks/event = %.1f",h_rate_4stub_pt->GetSum());
  mySmallText(0.5,0.85,1,txt);
  mySmallText(0.5,0.79,4,txt3);
  mySmallText(0.5,0.74,4,txt2);
  h_tp_pt->Draw("same");
  ctot.SaveAs("FakePlots/tot_vspt_4stub_"+type+".png");
  ctot.SaveAs("FakePlots/tot_vspt_4stub_"+type+".eps");

  TH1F* h_rate_4stubchi2_pt = (TH1F*) h_trk_4stubchi2_pt->Clone();
  h_rate_4stubchi2_pt->SetName("rate_4stubchi2_pt");
  h_rate_4stubchi2_pt->Scale(1.0/nevt);
  h_rate_4stubchi2_pt->Write();
  h_rate_4stubchi2_pt->GetYaxis()->SetTitle("Tracks / event");
  max = h_tp_pt->GetMaximum();
  if (h_rate_4stubchi2_pt->GetMaximum() > max) max = h_rate_4stubchi2_pt->GetMaximum();
  h_tp_pt->SetAxisRange(0,max*1.05,"Y");  
  h_tp_pt->Draw();
  h_rate_4stubchi2_pt->Draw("same");
  sprintf(txt,"average # tracks/event = %.1f",h_rate_4stubchi2_pt->GetSum());
  mySmallText(0.5,0.85,1,txt);
  mySmallText(0.5,0.79,4,txt3);
  mySmallText(0.5,0.74,4,txt2);
  h_tp_pt->Draw("same");
  ctot.SaveAs("FakePlots/tot_vspt_4stubchi2_"+type+".png");
  ctot.SaveAs("FakePlots/tot_vspt_4stubchi2_"+type+".eps");

  // 3 GeV

  TH1F* h_rate_pt3 = (TH1F*) h_trk_pt3->Clone();
  h_rate_pt3->SetName("rate_pt3");
  h_rate_pt3->Scale(1.0/nevt);
  h_rate_pt3->Write();
  h_rate_pt3->GetYaxis()->SetTitle("Tracks / event");
  max = h_tp_pt3->GetMaximum();
  if (h_rate_pt3->GetMaximum() > max) max = h_rate_pt3->GetMaximum();
  h_tp_pt3->SetAxisRange(0,max*1.05,"Y");  
  h_tp_pt3->Draw();
  h_rate_pt3->Draw("same");
  sprintf(txt,"average # tracks/event = %.1f",h_rate_pt3->GetSum());
  mySmallText(0.5,0.85,1,txt);
  sprintf(txt2,"event = %.1f",h_tp_pt3->GetSum());
  mySmallText(0.5,0.79,4,txt3);
  mySmallText(0.5,0.74,4,txt2);
  ctot.SaveAs("FakePlots/tot_vspt3_"+type+".png");
  ctot.SaveAs("FakePlots/tot_vspt3_"+type+".eps");

  TH1F* h_rate_4stub_pt3 = (TH1F*) h_trk_4stub_pt3->Clone();
  h_rate_4stub_pt3->SetName("rate_4stub_pt3");
  h_rate_4stub_pt3->Scale(1.0/nevt);
  h_rate_4stub_pt3->Write();
  h_rate_4stub_pt3->GetYaxis()->SetTitle("Tracks / event");
  max = h_tp_pt3->GetMaximum();
  if (h_rate_4stub_pt3->GetMaximum() > max) max = h_rate_4stub_pt3->GetMaximum();
  h_tp_pt3->SetAxisRange(0,max*1.05,"Y");  
  h_tp_pt3->Draw();
  h_rate_4stub_pt3->Draw("same");
  sprintf(txt,"average # tracks/event = %.1f",h_rate_4stub_pt3->GetSum());
  mySmallText(0.5,0.85,1,txt);
  mySmallText(0.5,0.79,4,txt3);
  mySmallText(0.5,0.74,4,txt2);
  h_tp_pt3->Draw("same");
  ctot.SaveAs("FakePlots/tot_vspt3_4stub_"+type+".png");
  ctot.SaveAs("FakePlots/tot_vspt3_4stub_"+type+".eps");

  TH1F* h_rate_4stubchi2_pt3 = (TH1F*) h_trk_4stubchi2_pt3->Clone();
  h_rate_4stubchi2_pt3->SetName("rate_4stubchi2_pt3");
  h_rate_4stubchi2_pt3->Scale(1.0/nevt);
  h_rate_4stubchi2_pt3->Write();
  h_rate_4stubchi2_pt3->GetYaxis()->SetTitle("Tracks / event");
  max = h_tp_pt3->GetMaximum();
  if (h_rate_4stubchi2_pt3->GetMaximum() > max) max = h_rate_4stubchi2_pt3->GetMaximum();
  h_tp_pt3->SetAxisRange(0,max*1.05,"Y");  
  h_tp_pt3->Draw();
  h_rate_4stubchi2_pt3->Draw("same");
  sprintf(txt,"average # tracks/event = %.1f",h_rate_4stubchi2_pt3->GetSum());
  mySmallText(0.5,0.85,1,txt);
  mySmallText(0.5,0.79,4,txt3);
  mySmallText(0.5,0.74,4,txt2);
  h_tp_pt3->Draw("same");
  ctot.SaveAs("FakePlots/tot_vspt3_4stubchi2_"+type+".png");
  ctot.SaveAs("FakePlots/tot_vspt3_4stubchi2_"+type+".eps");

  //overlays

  h_rate_pt->Draw("hist");
  h_rate_4stub_pt->SetLineColor(2);
  h_rate_4stub_pt->Draw("hist,same");
  h_rate_4stubchi2_pt->SetLineColor(4);
  h_rate_4stubchi2_pt->Draw("hist,same");
  ctot.SaveAs("FakePlots/tot_vspt_overlay_"+type+".png");
  ctot.SaveAs("FakePlots/tot_vspt_overlay_"+type+".eps");

  h_rate_pt3->Draw("hist");
  h_rate_4stub_pt3->SetLineColor(2);
  h_rate_4stub_pt3->Draw("hist,same");
  h_rate_4stubchi2_pt3->SetLineColor(4);
  h_rate_4stubchi2_pt3->Draw("hist,same");
  ctot.SaveAs("FakePlots/tot_vspt3_overlay_"+type+".png");
  ctot.SaveAs("FakePlots/tot_vspt3_overlay_"+type+".eps");



  // actual fakes
  TH1F* h_eff_pt = (TH1F*) h_trk_pt->Clone();
  h_eff_pt->SetName("eff_pt");
  h_eff_pt->GetYaxis()->SetTitle("Fake rate");
  h_eff_pt->Divide(h_trkFake_pt, h_trk_pt, 1.0, 1.0, "B");
  
  TH1F* h_eff_qc_pt = (TH1F*) h_trk_qc_pt->Clone();
  h_eff_qc_pt->SetName("eff_qc_pt");
  h_eff_qc_pt->GetYaxis()->SetTitle("Fake rate");
  h_eff_qc_pt->Divide(h_trkFake_qc_pt, h_trk_qc_pt, 1.0, 1.0, "B");

  h_eff_pt->Write();
  h_eff_qc_pt->Write();

  TH1F* h_eff_pt_C = (TH1F*) h_trk_pt_C->Clone();
  h_eff_pt_C->SetName("eff_pt_C");
  h_eff_pt_C->GetYaxis()->SetTitle("Fake rate");
  h_eff_pt_C->Divide(h_trkFake_pt_C, h_trk_pt_C, 1.0, 1.0, "B");
  
  TH1F* h_eff_qc_pt_C = (TH1F*) h_trk_qc_pt_C->Clone();
  h_eff_qc_pt_C->SetName("eff_qc_pt_C");
  h_eff_qc_pt_C->GetYaxis()->SetTitle("Fake rate");
  h_eff_qc_pt_C->Divide(h_trkFake_qc_pt_C, h_trk_qc_pt_C, 1.0, 1.0, "B");

  h_eff_pt_C->Write();
  h_eff_qc_pt_C->Write();


  TCanvas c;
  h_eff_pt->SetAxisRange(0.0,1.0,"Y");
  h_eff_pt->Draw();
  h_eff_qc_pt->SetLineColor(4);
  h_eff_qc_pt->SetMarkerColor(4);
  h_eff_qc_pt->SetMarkerStyle(24);
  h_eff_qc_pt->Draw("same");

  TLegend* ll = new TLegend(0.22,0.76,0.38,0.9);
  ll->SetFillStyle(0);
  ll->SetBorderSize(0);
  ll->SetTextSize(0.04);
  ll->AddEntry(h_eff_pt," No quality cuts","lep");
  ll->AddEntry(h_eff_qc_pt," #geq 4 stubs, #chi^{2} < 100","lep");
  ll->SetTextFont(42);
  ll->Draw();	

  c.SaveAs("FakePlots/fake_vspt_"+type+".png");
  c.SaveAs("FakePlots/fake_vspt_"+type+".eps");



  TH1F* h_eff_eta_pt2 = (TH1F*) h_trk_eta_pt2->Clone();
  h_eff_eta_pt2->SetName("eff_eta_pt2");
  h_eff_eta_pt2->GetYaxis()->SetTitle("Fake rate");
  h_eff_eta_pt2->Divide(h_trkFake_eta_pt2, h_trk_eta_pt2, 1.0, 1.0, "B");
  
  TH1F* h_eff_qc_eta_pt2 = (TH1F*) h_trk_qc_eta_pt2->Clone();
  h_eff_qc_eta_pt2->SetName("eff_qc_eta_pt2");
  h_eff_qc_eta_pt2->GetYaxis()->SetTitle("Fake rate");
  h_eff_qc_eta_pt2->Divide(h_trkFake_qc_eta_pt2, h_trk_qc_eta_pt2, 1.0, 1.0, "B");

  TH1F* h_eff_eta_pt5 = (TH1F*) h_trk_eta_pt5->Clone();
  h_eff_eta_pt5->SetName("eff_eta_pt5");
  h_eff_eta_pt5->GetYaxis()->SetTitle("Fake rate");
  h_eff_eta_pt5->Divide(h_trkFake_eta_pt5, h_trk_eta_pt5, 1.0, 1.0, "B");
  
  TH1F* h_eff_qc_eta_pt5 = (TH1F*) h_trk_qc_eta_pt5->Clone();
  h_eff_qc_eta_pt5->SetName("eff_qc_eta_pt5");
  h_eff_qc_eta_pt5->GetYaxis()->SetTitle("Fake rate");
  h_eff_qc_eta_pt5->Divide(h_trkFake_qc_eta_pt5, h_trk_qc_eta_pt5, 1.0, 1.0, "B");

  TH1F* h_eff_eta_pt10 = (TH1F*) h_trk_eta_pt10->Clone();
  h_eff_eta_pt10->SetName("eff_eta_pt10");
  h_eff_eta_pt10->GetYaxis()->SetTitle("Fake rate");
  h_eff_eta_pt10->Divide(h_trkFake_eta_pt10, h_trk_eta_pt10, 1.0, 1.0, "B");
  
  TH1F* h_eff_qc_eta_pt10 = (TH1F*) h_trk_qc_eta_pt10->Clone();
  h_eff_qc_eta_pt10->SetName("eff_qc_eta_pt10");
  h_eff_qc_eta_pt10->GetYaxis()->SetTitle("Fake rate");
  h_eff_qc_eta_pt10->Divide(h_trkFake_qc_eta_pt10, h_trk_qc_eta_pt10, 1.0, 1.0, "B");

  h_eff_eta_pt2->Write();
  h_eff_eta_pt5->Write();
  h_eff_eta_pt10->Write();

  h_eff_qc_eta_pt2->Write();
  h_eff_qc_eta_pt5->Write();
  h_eff_qc_eta_pt10->Write();



  h_eff_eta_pt2->SetAxisRange(0.0,1.3,"Y");
  h_eff_eta_pt2->Draw();
  h_eff_eta_pt5->SetLineColor(4);
  h_eff_eta_pt5->SetMarkerColor(4);
  h_eff_eta_pt5->SetMarkerStyle(8);
  h_eff_eta_pt5->Draw("same");
  h_eff_eta_pt10->SetLineColor(2);
  h_eff_eta_pt10->SetMarkerColor(2);
  h_eff_eta_pt10->SetMarkerStyle(24);
  h_eff_eta_pt10->Draw("same");

  TLegend* l2 = new TLegend(0.22,0.76,0.38,0.9);
  l2->SetFillStyle(0);
  l2->SetBorderSize(0);
  l2->SetTextSize(0.04);
  l2->AddEntry(h_eff_eta_pt2," pt > 2 GeV","lep");
  l2->AddEntry(h_eff_eta_pt5," pt > 5 GeV","lep");
  l2->AddEntry(h_eff_eta_pt10," pt > 10 GeV","lep");
  l2->SetTextFont(42);
  l2->Draw();	

  mySmallText(0.5,0.85,1,"w/o quality cuts");

  c.SaveAs("FakePlots/fake_vseta_"+type+".png");
  c.SaveAs("FakePlots/fake_vseta_"+type+".eps");


  h_eff_qc_eta_pt2->SetAxisRange(0,1.3,"Y");
  h_eff_qc_eta_pt2->Draw();
  h_eff_qc_eta_pt5->SetLineColor(4);
  h_eff_qc_eta_pt5->SetMarkerColor(4);
  h_eff_qc_eta_pt5->SetMarkerStyle(8);
  h_eff_qc_eta_pt5->Draw("same");
  h_eff_qc_eta_pt10->SetLineColor(2);
  h_eff_qc_eta_pt10->SetMarkerColor(2);
  h_eff_qc_eta_pt10->SetMarkerStyle(24);
  h_eff_qc_eta_pt10->Draw("same");

  l2->Draw();	

  mySmallText(0.5,0.85,1,"# stubs #geq 4, #chi^{2} < 100");

  c.SaveAs("FakePlots/fake_vseta_qc_"+type+".png");
  c.SaveAs("FakePlots/fake_vseta_qc_"+type+".eps");




  // non-prompt fraction
  TH1F* h_eff2_pt = (TH1F*) h_trk_pt->Clone();
  h_eff2_pt->SetName("effNP_pt");
  h_eff2_pt->GetYaxis()->SetTitle("Fake rate");
  h_eff2_pt->Divide(h_trkGoodNP_pt, h_trk_pt, 1.0, 1.0, "B");
  
  TH1F* h_eff2_qc_pt = (TH1F*) h_trk_qc_pt->Clone();
  h_eff2_qc_pt->SetName("effNP_qc_pt");
  h_eff2_qc_pt->GetYaxis()->SetTitle("Fake rate");
  h_eff2_qc_pt->Divide(h_trkGoodNP_qc_pt, h_trk_qc_pt, 1.0, 1.0, "B");

  h_eff2_pt->Write();
  h_eff2_qc_pt->Write();


  TH1F* h_eff2_eta_pt2 = (TH1F*) h_trk_eta_pt2->Clone();
  h_eff2_eta_pt2->SetName("effNP_eta_pt2");
  h_eff2_eta_pt2->GetYaxis()->SetTitle("GoodNP rate");
  h_eff2_eta_pt2->Divide(h_trkGoodNP_eta_pt2, h_trk_eta_pt2, 1.0, 1.0, "B");
  
  TH1F* h_eff2_qc_eta_pt2 = (TH1F*) h_trk_qc_eta_pt2->Clone();
  h_eff2_qc_eta_pt2->SetName("effNP_qc_eta_pt2");
  h_eff2_qc_eta_pt2->GetYaxis()->SetTitle("GoodNP rate");
  h_eff2_qc_eta_pt2->Divide(h_trkGoodNP_qc_eta_pt2, h_trk_qc_eta_pt2, 1.0, 1.0, "B");

  TH1F* h_eff2_eta_pt5 = (TH1F*) h_trk_eta_pt5->Clone();
  h_eff2_eta_pt5->SetName("effNP_eta_pt5");
  h_eff2_eta_pt5->GetYaxis()->SetTitle("GoodNP rate");
  h_eff2_eta_pt5->Divide(h_trkGoodNP_eta_pt5, h_trk_eta_pt5, 1.0, 1.0, "B");
  
  TH1F* h_eff2_qc_eta_pt5 = (TH1F*) h_trk_qc_eta_pt5->Clone();
  h_eff2_qc_eta_pt5->SetName("effNP_qc_eta_pt5");
  h_eff2_qc_eta_pt5->GetYaxis()->SetTitle("GoodNP rate");
  h_eff2_qc_eta_pt5->Divide(h_trkGoodNP_qc_eta_pt5, h_trk_qc_eta_pt5, 1.0, 1.0, "B");

  TH1F* h_eff2_eta_pt10 = (TH1F*) h_trk_eta_pt10->Clone();
  h_eff2_eta_pt10->SetName("effNP_eta_pt10");
  h_eff2_eta_pt10->GetYaxis()->SetTitle("GoodNP rate");
  h_eff2_eta_pt10->Divide(h_trkGoodNP_eta_pt10, h_trk_eta_pt10, 1.0, 1.0, "B");
  
  TH1F* h_eff2_qc_eta_pt10 = (TH1F*) h_trk_qc_eta_pt10->Clone();
  h_eff2_qc_eta_pt10->SetName("effNP_qc_eta_pt10");
  h_eff2_qc_eta_pt10->GetYaxis()->SetTitle("GoodNP rate");
  h_eff2_qc_eta_pt10->Divide(h_trkGoodNP_qc_eta_pt10, h_trk_qc_eta_pt10, 1.0, 1.0, "B");

  h_eff2_eta_pt2->Write();
  h_eff2_eta_pt5->Write();
  h_eff2_eta_pt10->Write();

  h_eff2_qc_eta_pt2->Write();
  h_eff2_qc_eta_pt5->Write();
  h_eff2_qc_eta_pt10->Write();



  // prompt fraction
  TH1F* h_eff3_pt = (TH1F*) h_trk_pt->Clone();
  h_eff3_pt->SetName("effP_pt");
  h_eff3_pt->GetYaxis()->SetTitle("Fake rate");
  h_eff3_pt->Divide(h_trkGoodP_pt, h_trk_pt, 1.0, 1.0, "B");
  
  TH1F* h_eff3_qc_pt = (TH1F*) h_trk_qc_pt->Clone();
  h_eff3_qc_pt->SetName("effP_qc_pt");
  h_eff3_qc_pt->GetYaxis()->SetTitle("Fake rate");
  h_eff3_qc_pt->Divide(h_trkGoodP_qc_pt, h_trk_qc_pt, 1.0, 1.0, "B");

  h_eff3_pt->Write();
  h_eff3_qc_pt->Write();


  TH1F* h_eff3_eta_pt2 = (TH1F*) h_trk_eta_pt2->Clone();
  h_eff3_eta_pt2->SetName("effP_eta_pt2");
  h_eff3_eta_pt2->GetYaxis()->SetTitle("GoodP rate");
  h_eff3_eta_pt2->Divide(h_trkGoodP_eta_pt2, h_trk_eta_pt2, 1.0, 1.0, "B");
  
  TH1F* h_eff3_qc_eta_pt2 = (TH1F*) h_trk_qc_eta_pt2->Clone();
  h_eff3_qc_eta_pt2->SetName("effP_qc_eta_pt2");
  h_eff3_qc_eta_pt2->GetYaxis()->SetTitle("GoodP rate");
  h_eff3_qc_eta_pt2->Divide(h_trkGoodP_qc_eta_pt2, h_trk_qc_eta_pt2, 1.0, 1.0, "B");

  TH1F* h_eff3_eta_pt5 = (TH1F*) h_trk_eta_pt5->Clone();
  h_eff3_eta_pt5->SetName("effP_eta_pt5");
  h_eff3_eta_pt5->GetYaxis()->SetTitle("GoodP rate");
  h_eff3_eta_pt5->Divide(h_trkGoodP_eta_pt5, h_trk_eta_pt5, 1.0, 1.0, "B");
  
  TH1F* h_eff3_qc_eta_pt5 = (TH1F*) h_trk_qc_eta_pt5->Clone();
  h_eff3_qc_eta_pt5->SetName("effP_qc_eta_pt5");
  h_eff3_qc_eta_pt5->GetYaxis()->SetTitle("GoodP rate");
  h_eff3_qc_eta_pt5->Divide(h_trkGoodP_qc_eta_pt5, h_trk_qc_eta_pt5, 1.0, 1.0, "B");

  TH1F* h_eff3_eta_pt10 = (TH1F*) h_trk_eta_pt10->Clone();
  h_eff3_eta_pt10->SetName("effP_eta_pt10");
  h_eff3_eta_pt10->GetYaxis()->SetTitle("GoodP rate");
  h_eff3_eta_pt10->Divide(h_trkGoodP_eta_pt10, h_trk_eta_pt10, 1.0, 1.0, "B");
  
  TH1F* h_eff3_qc_eta_pt10 = (TH1F*) h_trk_qc_eta_pt10->Clone();
  h_eff3_qc_eta_pt10->SetName("effP_qc_eta_pt10");
  h_eff3_qc_eta_pt10->GetYaxis()->SetTitle("GoodP rate");
  h_eff3_qc_eta_pt10->Divide(h_trkGoodP_qc_eta_pt10, h_trk_qc_eta_pt10, 1.0, 1.0, "B");

  h_eff3_eta_pt2->Write();
  h_eff3_eta_pt5->Write();
  h_eff3_eta_pt10->Write();

  h_eff3_qc_eta_pt2->Write();
  h_eff3_qc_eta_pt5->Write();
  h_eff3_qc_eta_pt10->Write();


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



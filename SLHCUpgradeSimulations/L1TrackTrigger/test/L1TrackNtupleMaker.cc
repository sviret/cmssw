//////////////////////////////////////////////////////////////////////
//                                                                  //
//  Analyzer for making mini-ntuple for L1 track performance plots  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"
#include "DataFormats/SiPixelDetId/interface/StackedTrackerDetId.h"

////////////////
// PHYSICS TOOLS
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/LorentzVector.h"


///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

//////////////
// NAMESPACES
using namespace std;
using namespace edm;


//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class L1TrackNtupleMaker : public edm::EDAnalyzer
{
public:

  // Constructor/destructor
  explicit L1TrackNtupleMaker(const edm::ParameterSet& iConfig);
  virtual ~L1TrackNtupleMaker();

  // Mandatory methods
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  
protected:
  
private:
  
  //-----------------------------------------------------------------------------------------------
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config; 
  
  int MyProcess;        // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool Slim;            // only keep the branches that we really need
  bool DebugMode;       // lots of debug printout statements
  bool SaveAllTracks;   // store in ntuples not only truth-matched tracks but ALL tracks
  bool SaveStubs;       // option to save also stubs in the ntuples (makes them large...)
  int L1Tk_nPar;        // use 4 or 5 parameter track fit? 
  int TP_minNStub;      // require TPs to have >= minNStub (defining efficiency denominator) (==0 means to only require >= 1 cluster)
  int TP_minNStubLayer; // require TPs to have stubs in >= minNStubLayer layers/disks (defining efficiency denominator)
  double TP_minPt;      // save TPs with pt > minPt 
  double TP_maxEta;     // save TPs with |eta| < maxEta 
  double TP_maxZ0;      // save TPs with |z0| < maxZ0 
  int L1Tk_minNStub;    // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)
  
  edm::InputTag L1TrackInputTag;        // L1 track collection
  edm::InputTag MCTruthTrackInputTag;   // MC truth collection


  bool TrackIsolation; // do track isolation?
  float PTmin; 
  float ETAmax; 
  float TrackPTmin;
  float TrackETAmax;
  float TrackChi2max;
  int TrackNStubmin;

  float IsoTrackZmax;
  float IsoTrackChi2max;
  float IsoTrackNStubmin;
  float IsoTrackPTmin;
  float IsoDRmin;
  float IsoDRmax;
  float IsoDZmax; 

  bool TrackingInJets; // do tracking in jets?

  bool PrimaryVertex; // save primary vertex information?


  //-----------------------------------------------------------------------------------------------
  // tree & branches for mini-ntuple

  TTree* eventTree;

  // all L1 tracks
  std::vector<float>* m_trk_pt;
  std::vector<float>* m_trk_eta;
  std::vector<float>* m_trk_phi;
  std::vector<float>* m_trk_d0;    //(filled if L1Tk_nPar==5, else 999)
  std::vector<float>* m_trk_z0;
  std::vector<float>* m_trk_chi2; 
  std::vector<float>* m_trk_consistency; 
  std::vector<int>*   m_trk_nstub;
  std::vector<int>*   m_trk_genuine;
  std::vector<int>*   m_trk_loose;
  std::vector<int>*   m_trk_unknown;
  std::vector<int>*   m_trk_combinatoric;
  std::vector<int>*   m_trk_fake;  //0 fake, 1 track from primary interaction, 2 secondary track
  std::vector<int>*   m_trk_matchtp_pdgid;
  std::vector<float>* m_trk_matchtp_pt;
  std::vector<float>* m_trk_matchtp_eta;
  std::vector<float>* m_trk_matchtp_phi;
  std::vector<float>* m_trk_matchtp_z0;
  std::vector<float>* m_trk_matchtp_dxy;
  std::vector<int>*   m_trk_injet;        //is the track within dR<0.4 of a genjet with pt > 30 GeV?
  std::vector<int>*   m_trk_injet_highpt; //is the track within dR<0.4 of a genjet with pt > 100 GeV?

  // all tracking particles
  std::vector<float>* m_tp_pt;
  std::vector<float>* m_tp_eta;
  std::vector<float>* m_tp_phi;
  std::vector<float>* m_tp_dxy;
  std::vector<float>* m_tp_d0;
  std::vector<float>* m_tp_z0;
  std::vector<float>* m_tp_d0_prod;
  std::vector<float>* m_tp_z0_prod;
  std::vector<int>*   m_tp_pdgid; 
  std::vector<int>*   m_tp_momid;          //PDG ID of the mother of the genparticle corresponding to the TP
  std::vector<int>*   m_tp_nmatch;
  std::vector<int>*   m_tp_nloosematch;
  std::vector<int>*   m_tp_nstub;          //total number of stubs associated with the TP
  std::vector<int>*   m_tp_nstublayer;     //number of layers/disks with at least one stub associated to the TP
  std::vector<int>*   m_tp_ngenstublayer;  //number of layers/disks with at least one GENUINE stub associated to the TP
  std::vector<int>*   m_tp_eventid;
  std::vector<int>*   m_tp_injet;
  std::vector<int>*   m_tp_injet_highpt;

  // *L1 track* properties if m_tp_nmatch > 0
  std::vector<float>* m_matchtrk_pt;
  std::vector<float>* m_matchtrk_eta;
  std::vector<float>* m_matchtrk_phi;
  std::vector<float>* m_matchtrk_d0; //this variable is only filled if L1Tk_nPar==5
  std::vector<float>* m_matchtrk_z0;
  std::vector<float>* m_matchtrk_chi2; 
  std::vector<float>* m_matchtrk_consistency; 
  std::vector<int>*   m_matchtrk_nstub;
  std::vector<int>*   m_matchtrk_injet;
  std::vector<int>*   m_matchtrk_injet_highpt;

  // *L1 track* properties if m_tp_nloosematch > 0
  std::vector<float>* m_loosematchtrk_pt;
  std::vector<float>* m_loosematchtrk_eta;
  std::vector<float>* m_loosematchtrk_phi;
  std::vector<float>* m_loosematchtrk_d0; //this variable is only filled if L1Tk_nPar==5
  std::vector<float>* m_loosematchtrk_z0;
  std::vector<float>* m_loosematchtrk_chi2; 
  std::vector<float>* m_loosematchtrk_consistency; 
  std::vector<int>*   m_loosematchtrk_nstub;
  std::vector<int>*   m_loosematchtrk_injet;
  std::vector<int>*   m_loosematchtrk_injet_highpt;

  // ALL stubs
  std::vector<float>* m_allstub_x;
  std::vector<float>* m_allstub_y;
  std::vector<float>* m_allstub_z;
  std::vector<float>* m_allstub_pt;
  std::vector<float>* m_allstub_ptsign;
  std::vector<int>*   m_allstub_isBarrel;
  std::vector<int>*   m_allstub_layer;
  std::vector<int>*   m_allstub_isPS;

  // isolation variable
  std::vector<float>* m_reliso;
  std::vector<float>* m_absiso;

  // primary vertex
  std::vector<float>* m_pv_L1;
  std::vector<float>* m_pv_MC;

  // "track jet variables" (for each gen jet, store the sum of pt of TPs / tracks inside jet cone)
  std::vector<float>* m_jet_eta;
  std::vector<float>* m_jet_pt;
  std::vector<float>* m_jet_tp_sumpt;
  std::vector<float>* m_jet_matchtrk_sumpt;
  std::vector<float>* m_jet_loosematchtrk_sumpt;
  std::vector<float>* m_jet_trk_sumpt;
  

};


//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
L1TrackNtupleMaker::L1TrackNtupleMaker(edm::ParameterSet const& iConfig) : 
  config(iConfig)
{

  MyProcess        = iConfig.getParameter< int >("MyProcess");
  Slim             = iConfig.getParameter< bool >("Slim");
  DebugMode        = iConfig.getParameter< bool >("DebugMode");
  SaveAllTracks    = iConfig.getParameter< bool >("SaveAllTracks");
  SaveStubs        = iConfig.getParameter< bool >("SaveStubs");
  L1Tk_nPar        = iConfig.getParameter< int >("L1Tk_nPar");
  TP_minNStub      = iConfig.getParameter< int >("TP_minNStub");
  TP_minNStubLayer = iConfig.getParameter< int >("TP_minNStubLayer");
  TP_minPt         = iConfig.getParameter< double >("TP_minPt");
  TP_maxEta        = iConfig.getParameter< double >("TP_maxEta");
  TP_maxZ0         = iConfig.getParameter< double >("TP_maxZ0");
  L1TrackInputTag      = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
  MCTruthTrackInputTag = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
  L1Tk_minNStub    = iConfig.getParameter< int >("L1Tk_minNStub");

  TrackIsolation = iConfig.getParameter< bool >("TrackIsolation");

  PTmin = (float)iConfig.getParameter<double>("PTmin");
  ETAmax = (float)iConfig.getParameter<double>("ETAmax");
  TrackPTmin = (float)iConfig.getParameter<double>("TrackPTmin");
  TrackETAmax = (float)iConfig.getParameter<double>("TrackETAmax");
  TrackChi2max = (float)iConfig.getParameter<double>("TrackChi2max");
  TrackNStubmin = (int)iConfig.getParameter<int>("TrackNStubmin");

  IsoTrackZmax = (float)iConfig.getParameter<double>("IsoTrackZmax");
  IsoTrackChi2max = (float)iConfig.getParameter<double>("IsoTrackChi2max");
  IsoTrackPTmin = (float)iConfig.getParameter<double>("IsoTrackPTmin");
  IsoTrackNStubmin = (int)iConfig.getParameter<int>("IsoTrackNStubmin");
  IsoDRmin = (float)iConfig.getParameter<double>("IsoDRmin");
  IsoDRmax = (float)iConfig.getParameter<double>("IsoDRmax");
  IsoDZmax = (float)iConfig.getParameter<double>("IsoDZmax");

  TrackingInJets = iConfig.getParameter< bool >("TrackingInJets");

  PrimaryVertex = iConfig.getParameter< bool >("PrimaryVertex");

}

/////////////
// DESTRUCTOR
L1TrackNtupleMaker::~L1TrackNtupleMaker()
{
}  

//////////
// END JOB
void L1TrackNtupleMaker::endJob()
{
  // things to be done at the exit of the event Loop
  cerr << "L1TrackNtupleMaker::endJob" << endl;

}

////////////
// BEGIN JOB
void L1TrackNtupleMaker::beginJob()
{

  // things to be done before entering the event Loop
  cerr << "L1TrackNtupleMaker::beginJob" << endl;


  //-----------------------------------------------------------------------------------------------
  // book histograms / make ntuple
  edm::Service<TFileService> fs;


  // initilize
  m_trk_pt    = new std::vector<float>;
  m_trk_eta   = new std::vector<float>;
  m_trk_phi   = new std::vector<float>;
  m_trk_z0    = new std::vector<float>;
  m_trk_d0    = new std::vector<float>;
  m_trk_chi2  = new std::vector<float>;
  m_trk_nstub = new std::vector<int>;
  m_trk_consistency  = new std::vector<float>;
  m_trk_genuine      = new std::vector<int>;
  m_trk_loose        = new std::vector<int>;
  m_trk_unknown      = new std::vector<int>;
  m_trk_combinatoric = new std::vector<int>;
  m_trk_fake = new std::vector<int>;
  m_trk_matchtp_pdgid = new std::vector<int>;
  m_trk_matchtp_pt = new std::vector<float>;
  m_trk_matchtp_eta = new std::vector<float>;
  m_trk_matchtp_phi = new std::vector<float>;
  m_trk_matchtp_z0 = new std::vector<float>;
  m_trk_matchtp_dxy = new std::vector<float>;
  m_trk_injet = new std::vector<int>;
  m_trk_injet_highpt = new std::vector<int>;

  m_tp_pt     = new std::vector<float>;
  m_tp_eta    = new std::vector<float>;
  m_tp_phi    = new std::vector<float>;
  m_tp_dxy    = new std::vector<float>;
  m_tp_d0     = new std::vector<float>;
  m_tp_z0     = new std::vector<float>;
  m_tp_d0_prod = new std::vector<float>;
  m_tp_z0_prod = new std::vector<float>;
  m_tp_pdgid  = new std::vector<int>;
  m_tp_momid  = new std::vector<int>;
  m_tp_nmatch = new std::vector<int>;
  m_tp_nloosematch = new std::vector<int>;
  m_tp_nstub  = new std::vector<int>;
  m_tp_nstublayer = new std::vector<int>;
  m_tp_ngenstublayer = new std::vector<int>;
  m_tp_eventid = new std::vector<int>;
  m_tp_injet = new std::vector<int>;
  m_tp_injet_highpt = new std::vector<int>;

  m_matchtrk_pt    = new std::vector<float>;
  m_matchtrk_eta   = new std::vector<float>;
  m_matchtrk_phi   = new std::vector<float>;
  m_matchtrk_z0    = new std::vector<float>;
  m_matchtrk_d0    = new std::vector<float>;
  m_matchtrk_chi2  = new std::vector<float>;
  m_matchtrk_nstub = new std::vector<int>;
  m_matchtrk_consistency = new std::vector<float>;
  m_matchtrk_injet = new std::vector<int>;
  m_matchtrk_injet_highpt = new std::vector<int>;
  
  m_loosematchtrk_pt    = new std::vector<float>;
  m_loosematchtrk_eta   = new std::vector<float>;
  m_loosematchtrk_phi   = new std::vector<float>;
  m_loosematchtrk_z0    = new std::vector<float>;
  m_loosematchtrk_d0    = new std::vector<float>;
  m_loosematchtrk_chi2  = new std::vector<float>;
  m_loosematchtrk_nstub = new std::vector<int>;
  m_loosematchtrk_consistency = new std::vector<float>;
  m_loosematchtrk_injet = new std::vector<int>;
  m_loosematchtrk_injet_highpt = new std::vector<int>;

  m_allstub_x = new std::vector<float>;
  m_allstub_y = new std::vector<float>;
  m_allstub_z = new std::vector<float>;
  m_allstub_pt     = new std::vector<float>;
  m_allstub_ptsign = new std::vector<float>;
  m_allstub_isBarrel = new std::vector<int>;
  m_allstub_layer    = new std::vector<int>;
  m_allstub_isPS     = new std::vector<int>;

  m_reliso = new std::vector<float>;
  m_absiso = new std::vector<float>;

  m_pv_L1 = new std::vector<float>;
  m_pv_MC = new std::vector<float>;

  m_jet_eta = new std::vector<float>;
  m_jet_pt = new std::vector<float>;
  m_jet_tp_sumpt = new std::vector<float>;
  m_jet_matchtrk_sumpt = new std::vector<float>;
  m_jet_loosematchtrk_sumpt = new std::vector<float>;
  m_jet_trk_sumpt = new std::vector<float>;



  // ntuple
  eventTree = fs->make<TTree>("eventTree", "Event tree");

  if (SaveAllTracks) {
    eventTree->Branch("trk_pt",    &m_trk_pt);
    eventTree->Branch("trk_eta",   &m_trk_eta);
    eventTree->Branch("trk_phi",   &m_trk_phi);
    eventTree->Branch("trk_chi2",  &m_trk_chi2);
    eventTree->Branch("trk_nstub", &m_trk_nstub);
    if (!Slim) {
      eventTree->Branch("trk_d0",    &m_trk_d0);
      eventTree->Branch("trk_z0",    &m_trk_z0);
      eventTree->Branch("trk_consistency",  &m_trk_consistency);
      eventTree->Branch("trk_genuine",      &m_trk_genuine);
      eventTree->Branch("trk_loose",        &m_trk_loose);
      eventTree->Branch("trk_unknown",      &m_trk_unknown);
      eventTree->Branch("trk_combinatoric", &m_trk_combinatoric);
      eventTree->Branch("trk_fake", &m_trk_fake);
      eventTree->Branch("trk_matchtp_pdgid",&m_trk_matchtp_pdgid);
      eventTree->Branch("trk_matchtp_pt",   &m_trk_matchtp_pt);
      eventTree->Branch("trk_matchtp_eta",  &m_trk_matchtp_eta);
      eventTree->Branch("trk_matchtp_phi",  &m_trk_matchtp_phi);
      eventTree->Branch("trk_matchtp_z0",   &m_trk_matchtp_z0);
      eventTree->Branch("trk_matchtp_dxy",  &m_trk_matchtp_dxy);
    }
    if (TrackingInJets) {
      eventTree->Branch("trk_injet", &m_trk_injet);
      eventTree->Branch("trk_injet_highpt", &m_trk_injet_highpt);
    }
  }

  eventTree->Branch("tp_pt",     &m_tp_pt);
  eventTree->Branch("tp_eta",    &m_tp_eta);
  eventTree->Branch("tp_phi",    &m_tp_phi);
  eventTree->Branch("tp_dxy",    &m_tp_dxy);
  eventTree->Branch("tp_d0",     &m_tp_d0);
  eventTree->Branch("tp_z0",     &m_tp_z0);
  eventTree->Branch("tp_pdgid",  &m_tp_pdgid);
  eventTree->Branch("tp_momid",  &m_tp_momid);
  eventTree->Branch("tp_nmatch", &m_tp_nmatch);
  eventTree->Branch("tp_nloosematch", &m_tp_nloosematch);
  eventTree->Branch("tp_nstub",  &m_tp_nstub);
  eventTree->Branch("tp_eventid",&m_tp_eventid);
  if (!Slim) {
    eventTree->Branch("tp_d0_prod",&m_tp_d0_prod);
    eventTree->Branch("tp_z0_prod",&m_tp_z0_prod);
    eventTree->Branch("tp_ngenstublayer", &m_tp_ngenstublayer);
    eventTree->Branch("tp_nstublayer",    &m_tp_nstublayer);
  }
  if (TrackingInJets) {
    eventTree->Branch("tp_injet",     &m_tp_injet);
    eventTree->Branch("tp_injet_highpt",     &m_tp_injet_highpt);
  }

  eventTree->Branch("matchtrk_pt",      &m_matchtrk_pt);
  eventTree->Branch("matchtrk_eta",     &m_matchtrk_eta);
  eventTree->Branch("matchtrk_phi",     &m_matchtrk_phi);
  eventTree->Branch("matchtrk_z0",      &m_matchtrk_z0);
  eventTree->Branch("matchtrk_d0",      &m_matchtrk_d0);
  eventTree->Branch("matchtrk_chi2",    &m_matchtrk_chi2);
  eventTree->Branch("matchtrk_nstub",   &m_matchtrk_nstub);
  if (!Slim) eventTree->Branch("matchtrk_consistency", &m_matchtrk_consistency);
  if (TrackingInJets) {
    eventTree->Branch("matchtrk_injet",    &m_matchtrk_injet);
    eventTree->Branch("matchtrk_injet_highpt",    &m_matchtrk_injet_highpt);
  }

  eventTree->Branch("loosematchtrk_pt",      &m_loosematchtrk_pt);
  eventTree->Branch("loosematchtrk_eta",     &m_loosematchtrk_eta);
  eventTree->Branch("loosematchtrk_phi",     &m_loosematchtrk_phi);
  eventTree->Branch("loosematchtrk_z0",      &m_loosematchtrk_z0);
  eventTree->Branch("loosematchtrk_d0",      &m_loosematchtrk_d0);
  eventTree->Branch("loosematchtrk_chi2",    &m_loosematchtrk_chi2);
  eventTree->Branch("loosematchtrk_nstub",   &m_loosematchtrk_nstub);
  if (!Slim) eventTree->Branch("loosematchtrk_consistency", &m_loosematchtrk_consistency);
  if (TrackingInJets) {
    eventTree->Branch("loosematchtrk_injet",   &m_loosematchtrk_injet);
    eventTree->Branch("loosematchtrk_injet_highpt",   &m_loosematchtrk_injet_highpt);
  }

  if (SaveStubs) {
    eventTree->Branch("allstub_x", &m_allstub_x);
    eventTree->Branch("allstub_y", &m_allstub_y);
    eventTree->Branch("allstub_z", &m_allstub_z);
    eventTree->Branch("allstub_pt",    &m_allstub_pt);
    eventTree->Branch("allstub_ptsign",&m_allstub_ptsign);
    eventTree->Branch("allstub_isBarrel", &m_allstub_isBarrel);
    eventTree->Branch("allstub_layer",    &m_allstub_layer);
    eventTree->Branch("allstub_isPS",     &m_allstub_isPS);
  }

  if (TrackIsolation) {
    eventTree->Branch("reliso", &m_reliso);
    eventTree->Branch("absiso", &m_absiso);
  }

  if (PrimaryVertex) {
    eventTree->Branch("pv_L1", &m_pv_L1);
    eventTree->Branch("pv_MC", &m_pv_MC);
  }

  if (TrackingInJets) {
    eventTree->Branch("jet_eta", &m_jet_eta);
    eventTree->Branch("jet_pt", &m_jet_pt);
    eventTree->Branch("jet_tp_sumpt", &m_jet_tp_sumpt);
    eventTree->Branch("jet_matchtrk_sumpt", &m_jet_matchtrk_sumpt);
    if (!Slim) eventTree->Branch("jet_loosematchtrk_sumpt", &m_jet_loosematchtrk_sumpt);
    eventTree->Branch("jet_trk_sumpt", &m_jet_trk_sumpt);
  }

}


//////////
// ANALYZE
void L1TrackNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if (DebugMode) cout << "start L1TrackNtupleMaker::analyze" << endl;

  if (!(MyProcess==13 || MyProcess==11 || MyProcess==211 || MyProcess==6 || MyProcess==15 || MyProcess==1)) {
    cout << "The specified MyProcess is invalid! Exiting..." << endl;
    return;
  }

  if ( !(L1Tk_nPar==4 || L1Tk_nPar==5) ) {
    cout << "Invalid number of track parameters, specified L1Tk_nPar == " << L1Tk_nPar << " but only 4/5 are valid options! Exiting..." << endl;
    return;
  }


  // clear variables
  m_trk_pt->clear();
  m_trk_eta->clear();
  m_trk_phi->clear();
  m_trk_d0->clear();
  m_trk_z0->clear();
  m_trk_chi2->clear();
  m_trk_nstub->clear();
  m_trk_consistency->clear();
  m_trk_genuine->clear();
  m_trk_loose->clear();
  m_trk_unknown->clear();
  m_trk_combinatoric->clear();
  m_trk_fake->clear();
  m_trk_matchtp_pdgid->clear();
  m_trk_matchtp_pt->clear();
  m_trk_matchtp_eta->clear();
  m_trk_matchtp_phi->clear();
  m_trk_matchtp_z0->clear();
  m_trk_matchtp_dxy->clear();
  m_trk_injet->clear();
  m_trk_injet_highpt->clear();
  
  m_tp_pt->clear();
  m_tp_eta->clear();
  m_tp_phi->clear();
  m_tp_dxy->clear();
  m_tp_d0->clear();
  m_tp_z0->clear();
  m_tp_d0_prod->clear();
  m_tp_z0_prod->clear();
  m_tp_pdgid->clear();
  m_tp_momid->clear();
  m_tp_nmatch->clear();
  m_tp_nloosematch->clear();
  m_tp_nstub->clear();
  m_tp_nstublayer->clear();
  m_tp_ngenstublayer->clear();
  m_tp_eventid->clear();
  m_tp_injet->clear();
  m_tp_injet_highpt->clear();

  m_matchtrk_pt->clear();
  m_matchtrk_eta->clear();
  m_matchtrk_phi->clear();
  m_matchtrk_z0->clear();
  m_matchtrk_d0->clear();
  m_matchtrk_chi2->clear();
  m_matchtrk_consistency->clear();
  m_matchtrk_nstub->clear();
  m_matchtrk_injet->clear();
  m_matchtrk_injet_highpt->clear();
  
  m_loosematchtrk_pt->clear();
  m_loosematchtrk_eta->clear();
  m_loosematchtrk_phi->clear();
  m_loosematchtrk_z0->clear();
  m_loosematchtrk_d0->clear();
  m_loosematchtrk_chi2->clear();
  m_loosematchtrk_consistency->clear();
  m_loosematchtrk_nstub->clear();
  m_loosematchtrk_injet->clear();
  m_loosematchtrk_injet_highpt->clear();
  
  m_allstub_x->clear();
  m_allstub_y->clear();
  m_allstub_z->clear();
  m_allstub_pt->clear();
  m_allstub_ptsign->clear();
  m_allstub_isBarrel->clear();
  m_allstub_layer->clear();
  m_allstub_isPS->clear();

  m_reliso->clear();
  m_absiso->clear();
  
  m_pv_L1->clear();
  m_pv_MC->clear();
  
  m_jet_eta->clear();
  m_jet_pt->clear();
  m_jet_tp_sumpt->clear();
  m_jet_matchtrk_sumpt->clear();
  m_jet_loosematchtrk_sumpt->clear();
  m_jet_trk_sumpt->clear();


  //-----------------------------------------------------------------------------------------------
  // retrieve various containers
  //-----------------------------------------------------------------------------------------------

  // L1 tracks
  if (DebugMode) cout << "get TTTracks" << endl;
  edm::Handle< std::vector< TTTrack< Ref_PixelDigi_ > > > TTTrackHandle;
  iEvent.getByLabel(L1TrackInputTag, TTTrackHandle);
  
  // L1 stubs
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > TTStubHandle;
  if (SaveStubs) {
    if (DebugMode) cout << "get TTStubs" << endl;
    iEvent.getByLabel("TTStubsFromPixelDigis", "StubAccepted", TTStubHandle);
  }

  // MC truth association maps
  if (DebugMode) cout << "get MC truth association maps" << endl;
  edm::Handle< TTTrackAssociationMap< Ref_PixelDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByLabel(MCTruthTrackInputTag, MCTruthTTTrackHandle);
  edm::Handle< TTClusterAssociationMap< Ref_PixelDigi_ > > MCTruthTTClusterHandle;
  iEvent.getByLabel("TTClusterAssociatorFromPixelDigis", "ClusterAccepted", MCTruthTTClusterHandle);
  edm::Handle< TTStubAssociationMap< Ref_PixelDigi_ > > MCTruthTTStubHandle;
  iEvent.getByLabel("TTStubAssociatorFromPixelDigis", "StubAccepted", MCTruthTTStubHandle);

  // tracking particles
  if (DebugMode) cout << "get tracking particles" << endl;
  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  edm::Handle< std::vector< TrackingVertex > > TrackingVertexHandle;
  iEvent.getByLabel("mix", "MergedTrackTruth", TrackingParticleHandle);
  iEvent.getByLabel("mix", "MergedTrackTruth", TrackingVertexHandle);



  // ----------------------------------------------------------------------------------------------
  // save primary vertex?
  // ----------------------------------------------------------------------------------------------

  if (PrimaryVertex) {

    if (DebugMode) cout << "getting L1/MC primary vertex" << endl;
    edm::Handle< std::vector<L1TkPrimaryVertex> > L1PVHandle;
    iEvent.getByLabel("L1TkPrimaryVertex","", L1PVHandle);
    
    edm::Handle< std::vector<L1TkPrimaryVertex> > MCPVHandle;
    iEvent.getByLabel("L1TkPrimaryVertexMC","", MCPVHandle);
    
    // by convention, the first vertex in the collection is the one that should be used by default
    
    float zvtx_L1 = -999;
    if (L1PVHandle.isValid()) {
      std::vector<L1TkPrimaryVertex>::const_iterator vtxIter = L1PVHandle->begin();
      zvtx_L1 = vtxIter->getZvertex();
    }
    float zvtx_MC = -999;
    if (MCPVHandle.isValid()) {
      std::vector<L1TkPrimaryVertex>::const_iterator vtxIterMC = MCPVHandle->begin();
      zvtx_MC = vtxIterMC->getZvertex();
    }

    m_pv_L1->push_back(zvtx_L1);
    m_pv_MC->push_back(zvtx_MC);
  }


  // ----------------------------------------------------------------------------------------------
  // do tracking in jets?
  // ----------------------------------------------------------------------------------------------

  std::vector<math::XYZTLorentzVector> v_jets;
  std::vector<int> v_jets_highpt;
  //std::vector<int> v_jets_tt;

  if (TrackingInJets) {

    // gen jets
    if (DebugMode) cout << "get genjets" << endl;
    edm::Handle< std::vector<reco::GenJet> > GenJetHandle;
    iEvent.getByLabel("ak4GenJets", "", GenJetHandle);

    
    if (GenJetHandle.isValid()) {
      
      if (DebugMode) cout << "loop over genjets" << endl;
      std::vector<reco::GenJet>::const_iterator iterGenJet;
      for ( iterGenJet = GenJetHandle->begin(); iterGenJet != GenJetHandle->end(); ++iterGenJet ) {
	
	reco::GenJet myJet = reco::GenJet(*iterGenJet);
	
	if (myJet.pt() < 30.0) continue;
	if (fabs(myJet.eta()) > 2.5) continue;
	
	if (DebugMode) cout << "genjet pt = " << myJet.pt() << ", eta = " << myJet.eta() << endl;

	bool ishighpt = false;
	if (myJet.pt() > 100.0)	ishighpt = true;

	/*
	bool istt = false;
	// loop over the constituents of the genjets
	if (DebugMode) cout << "loop over genjet constituents" << endl;

	std::vector<const reco::GenParticle*> jetConstituents = myJet.getGenConstituents();

	for (std::vector<const reco::GenParticle*>::const_iterator genCand = jetConstituents.begin(); genCand != jetConstituents.end(); ++genCand) {
	  
	  if (DebugMode) cout << "constituent pdgid = " << (*genCand)->pdgId() << endl;
	  
	  int nmoms = (int) (*genCand)->numberOfMothers();
	  if (nmoms > 0) {
	    const reco::Candidate* dau = (*genCand);
	    const reco::Candidate* mom = (*genCand)->mother();
	    while ( dau->numberOfMothers() == 1) {
	      if (DebugMode) cout << "constituent's mom's pdgid = " << mom->pdgId() << endl;
	      if (fabs(mom->pdgId()) == 6 || fabs(mom->pdgId()) == 24 || fabs(mom->pdgId()) == 5) istt = true;
	      dau = mom;
	      mom = dau->mother();
	    } 
	  }// end loop over jet gen constituents
	}
	*/

	math::XYZTLorentzVector jetP4 = myJet.p4();
	v_jets.push_back(jetP4);
	if (ishighpt) v_jets_highpt.push_back(1);
	else v_jets_highpt.push_back(0);
	//if (istt) v_jets_tt.push_back(1);
	//else v_jets_tt.push_back(0);
	
      }// end loop over genjets
    }// end isValid

  }// end TrackingInJets

  const int NJETS = 10;
  float jets_tp_sumpt[NJETS] = {0};       //sum pt of TPs with dR<0.4 of jet
  float jets_matchtrk_sumpt[NJETS] = {0}; //sum pt of tracks matched to TP with dR<0.4 of jet
  float jets_loosematchtrk_sumpt[NJETS] = {0}; //sum pt of tracks matched to TP with dR<0.4 of jet
  float jets_trk_sumpt[NJETS] = {0};      //sum pt of all tracks with dR<0.4 of jet


  // ----------------------------------------------------------------------------------------------
  // if saving stubs, need geometry handles etc.
  // ----------------------------------------------------------------------------------------------

  // geomtry 
  edm::ESHandle< StackedTrackerGeometry >         StackedGeometryHandle;
  const StackedTrackerGeometry*                   theStackedGeometry = 0;
  if (SaveStubs) {  
    iSetup.get< StackedTrackerGeometryRecord >().get(StackedGeometryHandle);
    theStackedGeometry = StackedGeometryHandle.product();
  }

  // magnetic field
  edm::ESHandle< MagneticField > magneticFieldHandle;
  iSetup.get< IdealMagneticFieldRecord >().get(magneticFieldHandle);
  double mMagneticFieldStrength = 0;
  if (SaveStubs) {  
    const MagneticField* theMagneticField = magneticFieldHandle.product();
    mMagneticFieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();
  }


  // ----------------------------------------------------------------------------------------------
  // loop over ALL L1 stubs
  // ----------------------------------------------------------------------------------------------

  if (SaveStubs) {
  
    if (DebugMode) cout << endl << "Loop over all stubs... " << endl;

    edmNew::DetSetVector<TTStub<Ref_PixelDigi_> >::const_iterator iterStubDet;
    for ( iterStubDet = TTStubHandle->begin(); iterStubDet != TTStubHandle->end(); ++iterStubDet ) {
      
      edmNew::DetSet<TTStub<Ref_PixelDigi_> >::const_iterator iterTTStub;
      for ( iterTTStub = iterStubDet->begin(); iterTTStub != iterStubDet->end(); ++iterTTStub ) {
	
	const TTStub<Ref_PixelDigi_>* stub=iterTTStub;
	      
	StackedTrackerDetId thisDetId = stub->getDetId();
	
	bool isPS = theStackedGeometry->isPSModule(thisDetId);
	bool isBarrel = thisDetId.isBarrel();

	int layer = -1;
	if (isBarrel) layer = thisDetId.iLayer();
	else layer = thisDetId.iDisk();
	  
	GlobalPoint posStub = theStackedGeometry->findGlobalPosition(stub);

	float tmp_stub_x = posStub.x();
	float tmp_stub_y = posStub.y();
	float tmp_stub_z = posStub.z();
    
	float tmp_stub_pt = theStackedGeometry->findRoughPt(mMagneticFieldStrength,stub);
	float tmp_stub_sign = 1.0;

	float trigBend = stub->getTriggerBend();
	if (trigBend<0 && tmp_stub_z<125.0) 
	  tmp_stub_sign = (-1)*tmp_stub_sign;
	else if (trigBend>0 && tmp_stub_z>125.0)
	  tmp_stub_sign = (-1)*tmp_stub_sign;

	m_allstub_x->push_back(tmp_stub_x);
	m_allstub_y->push_back(tmp_stub_y);
	m_allstub_z->push_back(tmp_stub_z);
	m_allstub_pt->push_back(tmp_stub_pt);
	m_allstub_ptsign->push_back(tmp_stub_sign);
	m_allstub_isPS->push_back(isPS);
	m_allstub_isBarrel->push_back(isBarrel);
	m_allstub_layer->push_back(layer);
	
      }
    }

  }// end if save stubs	


  // ----------------------------------------------------------------------------------------------
  // loop over L1 tracks
  // ----------------------------------------------------------------------------------------------

  if (SaveAllTracks) {
    
    if (DebugMode) {
      cout << endl << "Loop over L1 tracks!" << endl;
      cout << endl << "Looking at " << L1Tk_nPar << "-parameter tracks!" << endl;
    }
    
    int this_l1track = 0;
    std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterL1Track;
    for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
      
      edm::Ptr< TTTrack< Ref_PixelDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
      this_l1track++;
      
      float tmp_trk_pt   = iterL1Track->getMomentum(L1Tk_nPar).perp();
      float tmp_trk_eta  = iterL1Track->getMomentum(L1Tk_nPar).eta();
      float tmp_trk_phi  = iterL1Track->getMomentum(L1Tk_nPar).phi();
      float tmp_trk_z0   = iterL1Track->getPOCA(L1Tk_nPar).z(); //cm

      float tmp_trk_d0 = -999;
      if (L1Tk_nPar == 5) {
	float tmp_trk_x0   = iterL1Track->getPOCA(L1Tk_nPar).x();
	float tmp_trk_y0   = iterL1Track->getPOCA(L1Tk_nPar).y();	
	tmp_trk_d0 = -tmp_trk_x0*sin(tmp_trk_phi) + tmp_trk_y0*cos(tmp_trk_phi);
      }

      float tmp_trk_chi2 = iterL1Track->getChi2(L1Tk_nPar);
      float tmp_trk_consistency = iterL1Track->getStubPtConsistency(L1Tk_nPar);
      int tmp_trk_nstub  = (int) iterL1Track->getStubRefs().size();
            
      int tmp_trk_genuine = 0;
      int tmp_trk_loose = 0;
      int tmp_trk_unknown = 0;
      int tmp_trk_combinatoric = 0;
      if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) tmp_trk_loose = 1;
      if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
      if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) tmp_trk_unknown = 1;
      if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) tmp_trk_combinatoric = 1;
      
      if (DebugMode) {
	cout << "L1 track, pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi 
	     << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub;
	if (tmp_trk_genuine) cout << " (is genuine)" << endl; 
	if (tmp_trk_unknown) cout << " (is unknown)" << endl; 
	if (tmp_trk_combinatoric) cout << " (is combinatoric)" << endl; 
      }
      
      m_trk_pt ->push_back(tmp_trk_pt);
      m_trk_eta->push_back(tmp_trk_eta);
      m_trk_phi->push_back(tmp_trk_phi);
      m_trk_z0 ->push_back(tmp_trk_z0);
      if (L1Tk_nPar==5) m_trk_d0->push_back(tmp_trk_d0);
      else m_trk_d0->push_back(999.);
      m_trk_chi2 ->push_back(tmp_trk_chi2);
      m_trk_consistency->push_back(tmp_trk_consistency);
      m_trk_nstub->push_back(tmp_trk_nstub);
      m_trk_genuine->push_back(tmp_trk_genuine);
      m_trk_loose->push_back(tmp_trk_loose);
      m_trk_unknown->push_back(tmp_trk_unknown);
      m_trk_combinatoric->push_back(tmp_trk_combinatoric);
      

      // ----------------------------------------------------------------------------------------------
      // for studying the fake rate
      // ----------------------------------------------------------------------------------------------

      edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

      int myFake = 0;

      int myTP_pdgid = -999;
      float myTP_pt = -999;
      float myTP_eta = -999;
      float myTP_phi = -999;
      float myTP_z0 = -999;
      float myTP_dxy = -999;

      if (my_tp.isNull()) myFake = 0;
      else {
	int tmp_eventid = my_tp->eventId().event();

	if (tmp_eventid > 0) myFake = 2;
	else myFake = 1;

	myTP_pdgid = my_tp->pdgId();
	myTP_pt = my_tp->p4().pt();
	myTP_eta = my_tp->p4().eta();
	myTP_phi = my_tp->p4().phi();
	myTP_z0 = my_tp->vertex().z();
	
	float myTP_x0 = my_tp->vertex().x();
	float myTP_y0 = my_tp->vertex().y();
	myTP_dxy = sqrt(myTP_x0*myTP_x0 + myTP_y0*myTP_y0);

	if (DebugMode) {
	  cout << "TP matched to track has pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta() 
	       << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z()
	       << " pdgid = " <<  my_tp->pdgId() << " dxy = " << myTP_dxy << endl;
	}
      }

      m_trk_fake->push_back(myFake);

      m_trk_matchtp_pdgid->push_back(myTP_pdgid);
      m_trk_matchtp_pt->push_back(myTP_pt);
      m_trk_matchtp_eta->push_back(myTP_eta);
      m_trk_matchtp_phi->push_back(myTP_phi);
      m_trk_matchtp_z0->push_back(myTP_z0);
      m_trk_matchtp_dxy->push_back(myTP_dxy);


      // ----------------------------------------------------------------------------------------------
      // for tracking in jets
      // ----------------------------------------------------------------------------------------------

      if (TrackingInJets) {

	if (DebugMode) cout << "doing tracking in jets now" << endl;

	int InJet = 0;
	int InJetHighpt = 0;
	//int InJetTT = 0;

	for (int ij=0; ij<(int)v_jets.size(); ij++) {
	  float deta = tmp_trk_eta - (v_jets.at(ij)).eta();
	  float dphi = tmp_trk_phi - (v_jets.at(ij)).phi();
	  while (dphi > 3.14159) dphi = fabs(2*3.14159 - dphi);
	  float dR = sqrt(deta*deta + dphi*dphi);

	  if (dR < 0.4) {
	    InJet = 1;
	    //if (v_jets_tt.at(ij) == 1) InJetTT = 1;
	    if (v_jets_highpt.at(ij) == 1) InJetHighpt = 1;
	    if (ij<NJETS) jets_trk_sumpt[ij] += tmp_trk_pt;
	  }
	}

	m_trk_injet->push_back(InJet);
	m_trk_injet_highpt->push_back(InJetHighpt);

      }//end tracking in jets

    }//end track loop

  }//end if SaveAllTracks



  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------

  if (DebugMode) cout << endl << "Loop over tracking particles!" << endl;

  int this_tp = 0;
  std::vector< TrackingParticle >::const_iterator iterTP;
  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
 
    edm::Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, this_tp);
    this_tp++;

    int tmp_eventid = iterTP->eventId().event();
    if (MyProcess != 1 && tmp_eventid > 0) continue; //only care about tracking particles from the primary interaction (except for MyProcess==1, i.e. looking at all TPs)

    float tmp_tp_pt  = iterTP->pt();
    float tmp_tp_eta = iterTP->eta();
    float tmp_tp_phi = iterTP->phi(); 
    float tmp_tp_vz  = iterTP->vz();
    float tmp_tp_vx  = iterTP->vx();
    float tmp_tp_vy  = iterTP->vy();
    int tmp_tp_pdgid = iterTP->pdgId();
    float tmp_tp_z0_prod = tmp_tp_vz;
    float tmp_tp_d0_prod = -tmp_tp_vx*sin(tmp_tp_phi) + tmp_tp_vy*cos(tmp_tp_phi);

    if (MyProcess==13 && abs(tmp_tp_pdgid) != 13) continue;
    if (MyProcess==11 && abs(tmp_tp_pdgid) != 11) continue;
    if ((MyProcess==6 || MyProcess==15 || MyProcess==211) && abs(tmp_tp_pdgid) != 211) continue;

    if (tmp_tp_pt < TP_minPt) continue;
    if (fabs(tmp_tp_eta) > TP_maxEta) continue;


    // ----------------------------------------------------------------------------------------------
    // get the mother of the TP

    if (DebugMode) cout << "get mother of the TP" << endl;

    int tmp_tp_momid = -999;

    const reco::GenParticleRefVector & genp = iterTP->genParticles();
    for (reco::GenParticleRefVector::const_iterator ig=genp.begin(); ig!=genp.end(); ++ig) {

      if ((*ig)->numberOfMothers() > 0) {
	const reco::Candidate* mom = (*ig)->mother();
	if (tmp_tp_momid == -999) tmp_tp_momid = mom->pdgId();
	else {
	  cout << "MORE THAN ONE MOTHER !?!?" << endl;
	  break;
	}
      }
    }


    // ----------------------------------------------------------------------------------------------
    // get d0/z0 propagated back to the IP

    float tmp_tp_t = tan(2.0*atan(1.0)-2.0*atan(exp(-tmp_tp_eta)));

    float delx = -tmp_tp_vx;
    float dely = -tmp_tp_vy;
	
    float A = 0.01*0.5696;
    float Kmagnitude = A / tmp_tp_pt;
    
    float tmp_tp_charge = tp_ptr->charge();
    float K = Kmagnitude * tmp_tp_charge;
    float d = 0;
	
    float tmp_tp_x0p = delx - (d + 1./(2. * K)*sin(tmp_tp_phi));
    float tmp_tp_y0p = dely + (d + 1./(2. * K)*cos(tmp_tp_phi));
    float tmp_tp_rp = sqrt(tmp_tp_x0p*tmp_tp_x0p + tmp_tp_y0p*tmp_tp_y0p);
    float tmp_tp_d0 = tmp_tp_charge*tmp_tp_rp - (1. / (2. * K));
	
    tmp_tp_d0 = tmp_tp_d0*(-1); //fix d0 sign

    static double pi = 4.0*atan(1.0);
    float delphi = tmp_tp_phi-atan2(-K*tmp_tp_x0p,K*tmp_tp_y0p);
    if (delphi<-pi) delphi+=2.0*pi;
    if (delphi>pi) delphi-=2.0*pi;
    float tmp_tp_z0 = tmp_tp_vz+tmp_tp_t*delphi/(2.0*K);
    // ----------------------------------------------------------------------------------------------
    

    if (fabs(tmp_tp_z0) > TP_maxZ0) continue;


    // for pions in ttbar, only consider TPs coming from near the IP!
    float dxy = sqrt(tmp_tp_vx*tmp_tp_vx + tmp_tp_vy*tmp_tp_vy);
    float tmp_tp_dxy = dxy;
    if (MyProcess==6 && (dxy > 1.0)) continue;


    if (DebugMode) cout << "Tracking particle, pt: " << tmp_tp_pt << " eta: " << tmp_tp_eta << " phi: " << tmp_tp_phi 
			<< " z0: " << tmp_tp_z0 << " d0: " << tmp_tp_d0 
			<< " z_prod: " << tmp_tp_z0_prod << " d_prod: " << tmp_tp_d0_prod 
			<< " pdgid: " << tmp_tp_pdgid << " mother pdgid: " << tmp_tp_momid << " eventID: " << iterTP->eventId().event()
			<< " ttclusters " << MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size() 
			<< " ttstubs " << MCTruthTTStubHandle->findTTStubRefs(tp_ptr).size()
			<< " tttracks " << MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr).size() << endl;

    
    
    // ----------------------------------------------------------------------------------------------
    // only consider TPs associated with >= 1 cluster, or >= X stubs, or have stubs in >= X layers (configurable options)
    
    if (MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size() < 1) {
      if (DebugMode) cout << "No matching TTClusters for TP, continuing..." << endl;
      continue;
    }
    
    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
    int nStubTP = (int) theStubRefs.size(); 

    // how many layers/disks have stubs?
    int hasStubInLayer[11] = {0};
    for (unsigned int is=0; is<theStubRefs.size(); is++) {
      
      StackedTrackerDetId thisDetId( theStubRefs.at(is)->getDetId() );
      //GlobalPoint posStub = theStackedGeometry->findGlobalPosition( &(*theStubRefs.at(is)) );
      //bool isPS = theStackedGeometry->isPSModule(thisDetId);

      bool isBarrel = thisDetId.isBarrel();
      int layer = -1;
      if (isBarrel) layer = thisDetId.iLayer()-1; //fill in array as entries 0-5
      else layer = thisDetId.iDisk()+5;           //fill in array as entries 6-10

       //treat genuine stubs separately (==2 is genuine, ==1 is not)
      if (MCTruthTTStubHandle->findTrackingParticlePtr(theStubRefs.at(is)).isNull() && hasStubInLayer[layer]<2)
	hasStubInLayer[layer] = 1;
      else 
	hasStubInLayer[layer] = 2;
    }


    int nStubLayerTP = 0;
    int nStubLayerTP_g = 0;
    for (int isum=0; isum<11; isum++) {
      if ( hasStubInLayer[isum] >= 1) nStubLayerTP   += 1;
      if ( hasStubInLayer[isum] == 2) nStubLayerTP_g += 1;
    }


    if (DebugMode) cout << "TP is associated with " << nStubTP << " stubs, and has stubs in " << nStubLayerTP << " different layers/disks, and has GENUINE stubs in "
			<< nStubLayerTP_g << " layers " << endl;

    
    if (TP_minNStub > 0) {
      if (DebugMode) cout << "Only consider TPs with >= " << TP_minNStub << " stubs" << endl;
      if (nStubTP < TP_minNStub) {
	if (DebugMode) cout << "TP fails minimum nbr stubs requirement! Continuing..." << endl;
	continue;
      }
    }
    if (TP_minNStubLayer > 0) {
      if (DebugMode) cout << "Only consider TPs with stubs in >= " << TP_minNStubLayer << " layers/disks" << endl;
      if (nStubLayerTP < TP_minNStubLayer) {
	if (DebugMode) cout << "TP fails stubs in minimum nbr of layers/disks requirement! Continuing..." << endl;
	continue;
      }
    }

    

    // ----------------------------------------------------------------------------------------------
    // look for L1 tracks matched to the tracking particle

    std::vector< edm::Ptr< TTTrack< Ref_PixelDigi_ > > > matchedTracks = MCTruthTTTrackHandle->findTTTrackPtrs(tp_ptr);
    
    int nMatch = 0;
    int nLooseMatch = 0;
    int i_track = -1;
    int i_loosetrack = -1;
    float i_chi2dof = 99999;
    float i_loosechi2dof = 99999;

    if (matchedTracks.size() > 0) { 
    
      if (DebugMode && (matchedTracks.size()>1)) cout << "TrackingParticle has more than one matched L1 track!" << endl;


      // ----------------------------------------------------------------------------------------------
      // loop over matched L1 tracks
      // here, "match" means tracks that can be associated to a TrackingParticle with at least one hit of at least one of its clusters 
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SLHCTrackerTriggerSWTools#MC_truth_for_TTTrack

      for (int it=0; it<(int)matchedTracks.size(); it++) {

	bool tmp_trk_genuine = false;
	bool tmp_trk_loose = false;
	if (MCTruthTTTrackHandle->isGenuine(matchedTracks.at(it))) tmp_trk_genuine = true;
	if (MCTruthTTTrackHandle->isLooselyGenuine(matchedTracks.at(it))) tmp_trk_loose = true;
	if (!tmp_trk_loose) continue;


	if (DebugMode) {
	  if (MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it)).isNull()) {
	    cout << "track matched to TP is NOT uniquely matched to a TP" << endl;
	  }
	  else {
	    edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
	    cout << "TP matched to track matched to TP ... tp pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta() 
		 << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z() << endl;
	  }
	  cout << "   ... matched L1 track has pt = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).perp() 
	       << " eta = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).eta()
	       << " phi = " << matchedTracks.at(it)->getMomentum(L1Tk_nPar).phi()
	       << " chi2 = " << matchedTracks.at(it)->getChi2(L1Tk_nPar) 
	       << " consistency = " << matchedTracks.at(it)->getStubPtConsistency(L1Tk_nPar) 
	       << " z0 = " << matchedTracks.at(it)->getPOCA(L1Tk_nPar).z() 
	       << " nstub = " << matchedTracks.at(it)->getStubRefs().size();
	  if (tmp_trk_genuine) cout << " (genuine!) " << endl;
	  else if (tmp_trk_loose) cout << " (loosely genuine!) " << endl;
	  else cout << " (NOT loose/genuine) !!!!!" << endl;
	}


	// ----------------------------------------------------------------------------------------------
	// further require L1 track to be (loosely) genuine, that there is only one TP matched to the track
	// + have >= L1Tk_minNStub stubs for it to be a valid match (only relevant is your track collection
	// e.g. stores 3-stub tracks but at plot level you require >= 4 stubs (--> tracklet case)

	int tmp_trk_nstub = matchedTracks.at(it)->getStubRefs().size();

	if (tmp_trk_nstub < L1Tk_minNStub) continue;
	
	float dmatch_pt  = 999;
	float dmatch_eta = 999;
	float dmatch_phi = 999;
	int match_id = 999;

	edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(matchedTracks.at(it));
	dmatch_pt  = fabs(my_tp->p4().pt() - tmp_tp_pt);
	dmatch_eta = fabs(my_tp->p4().eta() - tmp_tp_eta);
	dmatch_phi = fabs(my_tp->p4().phi() - tmp_tp_phi);
	match_id = my_tp->pdgId();

	float tmp_trk_chi2dof = (matchedTracks.at(it)->getChi2(L1Tk_nPar)) / (2*tmp_trk_nstub - L1Tk_nPar);
	
	// ensure that track is uniquely matched to the TP we are looking at!
	if (dmatch_pt<0.1 && dmatch_eta<0.1 && dmatch_phi<0.1 && tmp_tp_pdgid==match_id) { 
	  nLooseMatch++;
	  if (i_loosetrack < 0 || tmp_trk_chi2dof < i_loosechi2dof) {
	    i_loosetrack = it;
	    i_loosechi2dof = tmp_trk_chi2dof;
	  }

	  if (tmp_trk_genuine) {
	    nMatch++;
	    if (i_track < 0 || tmp_trk_chi2dof < i_chi2dof) {
	      i_track = it;
	      i_chi2dof = tmp_trk_chi2dof;
	    }
	  }
	}


      }// end loop over matched L1 tracks

    }// end has at least 1 matched L1 track
    // ----------------------------------------------------------------------------------------------
        

    float tmp_matchtrk_pt   = -999;
    float tmp_matchtrk_eta  = -999;
    float tmp_matchtrk_phi  = -999;
    float tmp_matchtrk_z0   = -999;
    float tmp_matchtrk_d0   = -999;
    float tmp_matchtrk_chi2 = -999;
    float tmp_matchtrk_consistency = -999;
    int tmp_matchtrk_nstub  = -999;

    float tmp_loosematchtrk_pt   = -999;
    float tmp_loosematchtrk_eta  = -999;
    float tmp_loosematchtrk_phi  = -999;
    float tmp_loosematchtrk_z0   = -999;
    float tmp_loosematchtrk_d0   = -999;
    float tmp_loosematchtrk_chi2 = -999;
    float tmp_loosematchtrk_consistency = -999;
    int tmp_loosematchtrk_nstub  = -999;


    if (nMatch > 1 && DebugMode) cout << "WARNING *** 2 or more matches to genuine L1 tracks ***" << endl;
    if (nLooseMatch > 1 && DebugMode) cout << "WARNING *** 2 or more matches to loosely genuine L1 tracks ***" << endl;


    if (nMatch > 0) {
      tmp_matchtrk_pt   = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).perp();
      tmp_matchtrk_eta  = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).eta();
      tmp_matchtrk_phi  = matchedTracks.at(i_track)->getMomentum(L1Tk_nPar).phi();
      tmp_matchtrk_z0   = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).z();

      if (L1Tk_nPar == 5) {
	float tmp_matchtrk_x0 = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).x();
	float tmp_matchtrk_y0 = matchedTracks.at(i_track)->getPOCA(L1Tk_nPar).y();
	tmp_matchtrk_d0 = -tmp_matchtrk_x0*sin(tmp_matchtrk_phi) + tmp_matchtrk_y0*cos(tmp_matchtrk_phi);
      }

      tmp_matchtrk_chi2 = matchedTracks.at(i_track)->getChi2(L1Tk_nPar);
      tmp_matchtrk_consistency = matchedTracks.at(i_track)->getStubPtConsistency(L1Tk_nPar);
      tmp_matchtrk_nstub  = (int) matchedTracks.at(i_track)->getStubRefs().size();
    }

    if (nLooseMatch > 0) {
      tmp_loosematchtrk_pt   = matchedTracks.at(i_loosetrack)->getMomentum(L1Tk_nPar).perp();
      tmp_loosematchtrk_eta  = matchedTracks.at(i_loosetrack)->getMomentum(L1Tk_nPar).eta();
      tmp_loosematchtrk_phi  = matchedTracks.at(i_loosetrack)->getMomentum(L1Tk_nPar).phi();
      tmp_loosematchtrk_z0   = matchedTracks.at(i_loosetrack)->getPOCA(L1Tk_nPar).z();

      if (L1Tk_nPar == 5) {
	float tmp_loosematchtrk_x0 = matchedTracks.at(i_loosetrack)->getPOCA(L1Tk_nPar).x();
	float tmp_loosematchtrk_y0 = matchedTracks.at(i_loosetrack)->getPOCA(L1Tk_nPar).y();
	tmp_loosematchtrk_d0 = -tmp_loosematchtrk_x0*sin(tmp_loosematchtrk_phi) + tmp_loosematchtrk_y0*cos(tmp_loosematchtrk_phi);
      }

      tmp_loosematchtrk_chi2 = matchedTracks.at(i_loosetrack)->getChi2(L1Tk_nPar);
      tmp_loosematchtrk_consistency = matchedTracks.at(i_loosetrack)->getStubPtConsistency(L1Tk_nPar);
      tmp_loosematchtrk_nstub  = (int) matchedTracks.at(i_loosetrack)->getStubRefs().size();
    }


    m_tp_pt->push_back(tmp_tp_pt);
    m_tp_eta->push_back(tmp_tp_eta);
    m_tp_phi->push_back(tmp_tp_phi);
    m_tp_dxy->push_back(tmp_tp_dxy);
    m_tp_z0->push_back(tmp_tp_z0);
    m_tp_d0->push_back(tmp_tp_d0);
    m_tp_z0_prod->push_back(tmp_tp_z0_prod);
    m_tp_d0_prod->push_back(tmp_tp_d0_prod);
    m_tp_pdgid->push_back(tmp_tp_pdgid);
    m_tp_momid->push_back(tmp_tp_momid);
    m_tp_nmatch->push_back(nMatch);
    m_tp_nloosematch->push_back(nLooseMatch);
    m_tp_nstub->push_back(nStubTP);
    m_tp_nstublayer->push_back(nStubLayerTP);
    m_tp_ngenstublayer->push_back(nStubLayerTP_g);
    m_tp_eventid->push_back(tmp_eventid);

    m_matchtrk_pt ->push_back(tmp_matchtrk_pt);
    m_matchtrk_eta->push_back(tmp_matchtrk_eta);
    m_matchtrk_phi->push_back(tmp_matchtrk_phi);
    m_matchtrk_z0 ->push_back(tmp_matchtrk_z0);
    m_matchtrk_d0 ->push_back(tmp_matchtrk_d0);
    m_matchtrk_chi2 ->push_back(tmp_matchtrk_chi2);
    m_matchtrk_consistency->push_back(tmp_matchtrk_consistency);
    m_matchtrk_nstub->push_back(tmp_matchtrk_nstub);

    m_loosematchtrk_pt ->push_back(tmp_loosematchtrk_pt);
    m_loosematchtrk_eta->push_back(tmp_loosematchtrk_eta);
    m_loosematchtrk_phi->push_back(tmp_loosematchtrk_phi);
    m_loosematchtrk_z0 ->push_back(tmp_loosematchtrk_z0);
    m_loosematchtrk_d0 ->push_back(tmp_loosematchtrk_d0);
    m_loosematchtrk_chi2 ->push_back(tmp_loosematchtrk_chi2);
    m_loosematchtrk_consistency->push_back(tmp_loosematchtrk_consistency);
    m_loosematchtrk_nstub->push_back(tmp_loosematchtrk_nstub);


    // ----------------------------------------------------------------------------------------------
    // for tracking in jets 
    // ----------------------------------------------------------------------------------------------

    if (TrackingInJets) {

      if (DebugMode) cout << "check if TP/matched track is within jet" << endl;

      int tp_InJet = 0;
      int matchtrk_InJet = 0;
      int loosematchtrk_InJet = 0;
      int tp_InJetHighpt = 0;
      int matchtrk_InJetHighpt = 0;
      int loosematchtrk_InJetHighpt = 0;
      //int tp_InJetTT = 0;
      //int matchtrk_InJetTT = 0;
      //int loosematchtrk_InJetTT = 0;
      
      for (int ij=0; ij<(int)v_jets.size(); ij++) {
	float deta = tmp_tp_eta - (v_jets.at(ij)).eta();
	float dphi = tmp_tp_phi - (v_jets.at(ij)).phi();
	while (dphi > 3.14159) dphi = fabs(2*3.14159 - dphi);
	float dR = sqrt(deta*deta + dphi*dphi);
	if (dR < 0.4) {
	  tp_InJet = 1;
	  //if (v_jets_tt.at(ij) == 1) tp_InJetTT = 1;
	  if (v_jets_highpt.at(ij) == 1) tp_InJetHighpt = 1;
	  if (ij<NJETS) jets_tp_sumpt[ij] += tmp_tp_pt;
	}

	if (nMatch > 0) {
	  deta = tmp_matchtrk_eta - (v_jets.at(ij)).eta();
	  dphi = tmp_matchtrk_phi - (v_jets.at(ij)).phi();
	  while (dphi > 3.14159) dphi = fabs(2*3.14159 - dphi);
	  dR = sqrt(deta*deta + dphi*dphi); 
	  if (dR < 0.4) {
	    matchtrk_InJet = 1;
	    //if (v_jets_tt.at(ij) == 1) matchtrk_InJetTT = 1;
	    if (v_jets_highpt.at(ij) == 1) matchtrk_InJetHighpt = 1;
	    if (ij<NJETS) jets_matchtrk_sumpt[ij] += tmp_matchtrk_pt;
	  }
	}

	if (nLooseMatch > 0) {
	  deta = tmp_loosematchtrk_eta - (v_jets.at(ij)).eta();
	  dphi = tmp_loosematchtrk_phi - (v_jets.at(ij)).phi();
	  while (dphi > 3.14159) dphi = fabs(2*3.14159 - dphi);
	  dR = sqrt(deta*deta + dphi*dphi); 
	  if (dR < 0.4) {
	    loosematchtrk_InJet = 1;
	    //if (v_jets_tt.at(ij) == 1) loosematchtrk_InJetTT = 1;
	    if (v_jets_highpt.at(ij) == 1) loosematchtrk_InJetHighpt = 1;
	    if (ij<NJETS) jets_loosematchtrk_sumpt[ij] += tmp_loosematchtrk_pt;
	  }
	}
      }
      
      m_tp_injet->push_back(tp_InJet);
      m_tp_injet_highpt->push_back(tp_InJetHighpt);
      m_matchtrk_injet->push_back(matchtrk_InJet);
      m_matchtrk_injet_highpt->push_back(matchtrk_InJetHighpt);
      m_loosematchtrk_injet->push_back(loosematchtrk_InJet);
      m_loosematchtrk_injet_highpt->push_back(loosematchtrk_InJetHighpt);

    }//end TrackingInJets


    // ----------------------------------------------------------------------------------------------
    // track isolation 
    // ----------------------------------------------------------------------------------------------

    float absiso = -999;
    float reliso = -999;

    if (TrackIsolation) {
      
      if (DebugMode) cout << "calculate track isolation" << endl;

      bool TrackForIsolation = true;

      if (tmp_tp_pt < PTmin) TrackForIsolation = false;
      if (fabs(tmp_tp_eta) > ETAmax) TrackForIsolation = false;  
      if (nMatch == 0) TrackForIsolation = false;
      
      if ((tmp_matchtrk_pt < TrackPTmin) || (fabs(tmp_matchtrk_eta) > TrackETAmax) || (tmp_matchtrk_chi2 > TrackChi2max) || (tmp_matchtrk_nstub < TrackNStubmin) ) 
	TrackForIsolation = false;
      
      
      // ----------------------------------------------------------------------------------------------
      // now continue with calculating isolation variable
      // ----------------------------------------------------------------------------------------------

      if (TrackForIsolation) {
      
	float sumPt = 0;
	
	// loop over L1 tracks
	int ncentral = 0;
	int this_l1track = 0;
	std::vector< TTTrack< Ref_PixelDigi_ > >::const_iterator iterL1Track;
	for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
	  
	  edm::Ptr< TTTrack< Ref_PixelDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
	  this_l1track++;
	  
	  float iso_trk_pt   = iterL1Track->getMomentum(L1Tk_nPar).perp();
	  float iso_trk_eta  = iterL1Track->getMomentum(L1Tk_nPar).eta();
	  float iso_trk_phi  = iterL1Track->getMomentum(L1Tk_nPar).phi();
	  float iso_trk_z0   = iterL1Track->getPOCA(L1Tk_nPar).z(); //cm
	  float iso_trk_chi2 = iterL1Track->getChi2(L1Tk_nPar);
	  int iso_trk_nstub  = (int) iterL1Track->getStubRefs().size();

	  // avoid the central track
	  if ( fabs(iso_trk_pt-tmp_matchtrk_pt)<0.01 && fabs(iso_trk_eta-tmp_matchtrk_eta)<0.01 && fabs(iso_trk_phi-tmp_matchtrk_phi)<0.01) {
	    ncentral++;
	    if (ncentral == 1) continue; //if there are duplicate tracks, we don't want to ignore those
	  }


	  if (iso_trk_z0 > IsoTrackZmax) continue;
	  if (iso_trk_chi2 > IsoTrackChi2max) continue;
	  if (iso_trk_nstub < IsoTrackNStubmin) continue;
	  if (iso_trk_pt < IsoTrackPTmin) continue;
	  
	  float dz = fabs(iso_trk_z0 - tmp_matchtrk_z0);

	  if (dz > IsoDZmax) continue;
	  
	  float deta = iso_trk_eta - tmp_matchtrk_eta;
	  float dphi = iso_trk_phi - tmp_matchtrk_phi;
	  while (dphi > 3.14159) dphi = fabs(2*3.14159 - dphi);
	  float dR = sqrt(deta*deta + dphi*dphi);
	  
	  if (dR < IsoDRmin || dR > IsoDRmax) continue;
	  
	  sumPt += iso_trk_pt;
	  
	}// end track loop

	absiso = sumPt;
	reliso = sumPt / tmp_matchtrk_pt;

	if (DebugMode) cout << "rel iso = " << reliso << " absiso = " << absiso << endl;

      }
    
      m_absiso->push_back(absiso);
      m_reliso->push_back(reliso);

    }//end TrackIsolation
    
    
  } //end loop tracking particles
    
  
  if (TrackingInJets) {
    for (int ij=0; ij<(int)v_jets.size(); ij++) {
      if (ij<NJETS) {
	m_jet_eta->push_back((v_jets.at(ij)).eta());
	m_jet_pt->push_back((v_jets.at(ij)).pt());
	m_jet_tp_sumpt->push_back(jets_tp_sumpt[ij]);
	m_jet_matchtrk_sumpt->push_back(jets_matchtrk_sumpt[ij]);
	m_jet_loosematchtrk_sumpt->push_back(jets_loosematchtrk_sumpt[ij]);
	m_jet_trk_sumpt->push_back(jets_trk_sumpt[ij]);
      }
    }
  }


  eventTree->Fill();


} // end of analyze()


///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1TrackNtupleMaker);

#ifndef _PATTERNGENERATOR_H_
#define _PATTERNGENERATOR_H_

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <TChain.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TFile.h>
#include <TSystemDirectory.h>
#include "SectorTree.h"

using namespace std;

/**
   \brief Creates a pattern bank from root muons simulation files.
   Each event contains only one track. To be used in the pattern creation, it must have at least one stub per used layer (except in the case of fake stubs). 
**/
class PatternGenerator{
 private:
  map<int,int> variableRes; // number of DC bits used for each layer
  bool variableRes_state_cache;
  bool cache_is_uptodate;
  float ptMin;
  float ptMax;
  float etaMin;
  float etaMax;
  int nbMaxFakeSuperstrips;
  vector<int> tracker_layers;
  vector<int> inactive_layers;
  string particuleDirName;

  // Containers to load the TTree branches
  int m_stub;
  vector<int>                   m_stub_modid;  // layer/ladder/module of the stub
  vector<int>                   m_stub_pdg; // PDG ID of the particle 
  vector<float>                 m_stub_strip;  // Strip du cluster interne du stub
  vector<float>                 m_stub_ptGEN;  // PT de la particule initiale
  vector<float>                 m_stub_etaGEN;  // Eta de la particule initiale

  vector<int>                   *p_m_stub_modid;
  vector<int>                   *p_m_stub_pdg;
  vector<float>                 *p_m_stub_strip;
  vector<float>                 *p_m_stub_ptGEN;
  vector<float>                 *p_m_stub_etaGEN;

  TChain* createTChain(string directoryName, string tchainName);
  /**
     If coverageEstimation!=NULL we do not create patterns, we just test the coverage of the existing bank with new tracks
   **/
  int generate(TChain* TT, int* evtIndex, int evtNumber, int* nbTrack, SectorTree* sectors, map<int,pair<float,float> > eta_limits, int* coverageEstimation=NULL);
  
 public:
 /**
     \brief Constructor
  **/
  PatternGenerator();

  /**
     \brief Change the minimum PT accepted to create a pattern (default is 2)
     \param minp The minimum PT in GeV/c
  **/
  void setMinPT(float minp);
  /**
     \brief Change the maximum PT accepted to create a pattern (default is 100)
     \param maxp The maximum PT in GeV/c
  **/
  void setMaxPT(float maxp);
 /**
     \brief Change the minimum Eta accepted to create a pattern
     \param mine The minimum Eta
  **/
  void setMinEta(float mine);
  /**
     \brief Change the maximum Eta accepted to create a pattern
     \param maxe The maximum Eta
  **/
  void setMaxEta(float maxe);
  /**
     \brief Change the maximum number of fake superstrips that can be used in a pattern
     \param mf The maximum number of fake superstrips in a pattern
  **/
  void setMaxFakeSuperstrips(int mf);
 /**
     \brief Sets the layers used
     \param l A vector of int, each integer is a layer number (ex : 8 9 10 for the last 3 layers)
  **/
  void setLayers(vector<int> l);
 /**
     \brief Sets the inactive layers (will have only fake stubs but will be present in the bank)
     \param l A vector of int, each integer is a layer number (ex : 8 9 10 for the last 3 layers)
  **/
  void setInactiveLayers(vector<int> l);
  /**
     \brief Set the name of the directory containing the root files with muons informations
     \param f The name of the directory
  **/
  void setParticuleDirName(string f);
  /**
     \brief Set the number of DC bits used (0 to 3). If 0, adaptative patterns are not used
     \param nb The number of DC bits to use
     \param l The ID of the layer (0 by default)
   **/
  void setVariableResolution(int nb, int l=0);
  /**
     \brief Is the variable resolutions pattern system activated
     \return True if activated
  **/
  bool getVariableResolutionState();
  /**
     \brief Generates the patterns
     \param sectors This structure contains the interesting sectors, patterns will be added to this structure
     \param step Number of events used at each step
     \param threshold The process will stop when the size of the bank does not grow above this threshold (%) beteen 2 steps (or when we don't have any event).
     \param eta_limits For each sector the minimum eta and the maximum eta on which the sector does exist. If a sector is missing in the list, no limits are taken into account.
  **/
  void generate(SectorTree* sectors, int step, float threshold, map<int,pair<float,float> > eta_limits);
};
#endif

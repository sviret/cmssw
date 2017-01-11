#ifndef _PATTERNTRUNK_H_
#define _PATTERNTRUNK_H_

#include <map>
#include <vector>
#include <cmath>
#include <cstring>
#ifdef __APPLE__
#include <boost/mpl/and.hpp> // needed for boost v1.35
#endif
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include "GradedPattern.h"

#ifdef IPNL_USE_CUDA
#include "gpu_struct.h"
#endif

using namespace std;

/**
   \brief A PatternTrunk can contain one low definition pattern and all the associated full definitions patterns (all high definition patterns contnained have the same low resolution version). Used to store patterns and compute variable resolution patterns : we keep the low resolution versions and compute the DC bits values from the high resolution patterns.
**/

class PatternTrunk{
 public:
  /**
     \brief Constructor
     \param p The low definition pattern (will be copied)
  **/
  PatternTrunk(Pattern* p);

  PatternTrunk(GradedPattern* p);

  /**
     \brief Default constructor
  **/
  PatternTrunk();
  ~PatternTrunk();
  /**
     \brief Add a new full definition pattern to the structure
     Increment the low definition pattern grade. If the FD pattern is already contained, increments the grade.
  **/
  void addFDPattern(Pattern* p);

  /**
     \brief Add a new full definition pattern to the structure and update the average pt of this pattern
     Increment the low definition pattern grade. If the FD pattern is already contained, increments the grade.
  **/
  void addFDPattern(Pattern* p, float pt, int pdg);


  /**
     \brief Get the list of all the full definition patterns
     \return A vector of pointers on the Patterns (each Pattern is a copy)
  **/
  vector<GradedPattern*> getFDPatterns();
 /**
     \brief Get the low definition pattern
     \return A pointer on a copy of the Pattern
  **/
  GradedPattern* getLDPattern();
 /**
     \brief Get the average PT of the tracks that have created the low definition pattern
     \return The PT in GeV/c
  **/
  float getLDPatternPT();
 /**
     \brief Get the number of tracks that have created the low definition pattern
     \return The number of tracks leading to this pattern
  **/
  int getLDPatternGrade() const;
  /**
     \brief Get the number of FD patterns
     \return The number of FD patterns stored in the PatternTrunk
  **/
  int getFDPatternNumber();

  /**
     \brief Set the DC bits of the LD patterns. All FD patterns are removed.
     \param r The number of DC bits used between FD and LD for each layer
  **/
  void computeAdaptativePattern(vector<int> r);

  /**
     \brief Link the low definition patterns to the detector structure
     \param d The detector
  **/  
  void link(Detector& d);

#ifdef IPNL_USE_CUDA
 /**
     \brief Link the low definition patterns to the detector structure
     \param p the pattern bank structure on the device
     \param d The detector structure on the device
     \param pattern_index The index of the pattern in the bank
     \param sec The ladders in the sector
     \param modules The modules in the sector (one vector per ladder)
     \param layers List of layers IDs
  **/  
  void linkCuda(patternBank* p, deviceDetector* d, int pattern_index, const vector< vector<int> >& sec, const vector<map<int, vector<int> > >& modules, vector<int> layers, unsigned int* cache);
#endif

  /**
     \brief Returns a copy of the active pattern
     \param active_threshold The minimum number of active super strips to activate the pattern
     \return A pointer on the copy
  **/
  GradedPattern* getActivePattern(int active_threshold);
  /**
     \brief Returns a copy of the active pattern
     \param max_nb_missing_hit The maximum number of non active layers to activate the pattern 
     \param active_threshold The minimum number of active super strips to activate the pattern
     \return A pointer on the copy
  **/
  GradedPattern* getActivePatternUsingMissingHit(int max_nb_missing_hit, int active_threshold);

  /**
     \brief Check if the high resolution pattern is already in the bank when DC bits are activated
     \param hp The attern to check
     \result True if the pattern is already in tha bank, false otherwise
   **/
  bool checkPattern(Pattern* hp);

  /**
     \brief Set the orderInChip information to the LD pattern
     \param i The theoretical order in the chip
   **/
  void setOrderInChip(int i);
  /**
     \brief Get the orderInChip information of the LD pattern
     \return The theoretical order in the chip
   **/
  int getOrderInChip() const;
  /**
     \brief Update the data and the DC bits of the LDP according to the parameter
     \param gp The new pattern
   **/
  void updateWithPattern(GradedPattern* gp);

 private:
  GradedPattern* lowDefPattern;
  map<string, GradedPattern*> fullDefPatterns;

  void deleteFDPatterns();
  /**
     \brief Change the DC bits of the LDP to take the parameter's DC bits into account (used while merging banks)
     \param p A new pattern
  **/
  void updateDCBits(GradedPattern* p);
  /**
     \brief Update the data of the LDP with the parameters of the new pattern
     \param p A new pattern
  **/
  void updatePatternData(const GradedPattern& p);

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const{
    ar << lowDefPattern;
    ar << fullDefPatterns;
  }
  
  template<class Archive> void load(Archive & ar, const unsigned int version){
    if(lowDefPattern!=NULL)
      delete lowDefPattern;
    ar >> lowDefPattern;
    ar >> fullDefPatterns;
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

};
#endif

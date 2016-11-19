//This class implementes the duplicate removal
#ifndef FPGAPURGEDUPLICATE_H
#define FPGAPURGEDUPLICATE_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAPurgeDuplicate:public FPGAProcessBase{

public:

  FPGAPurgeDuplicate(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
           << " to output "<<output<<endl;
    }
    if (output=="trackout1"||
        output=="trackout2"||
        output=="trackout3"||
        output=="trackout4"||
        output=="trackout5"||
        output=="trackout6"||
        output=="trackout7"||
        output=="trackout8"||
        output=="trackout9"||
        output=="trackout10"){
    FPGACleanTrack* tmp=dynamic_cast<FPGACleanTrack*>(memory);
    assert(tmp!=0);
    outputtracklets_.push_back(tmp);
    return;
    }
    cout << "Did not find output : "<<output<<endl;
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
           << " to input "<<input<<endl;
    }
    if (input=="trackin1"||
        input=="trackin2"||
        input=="trackin3"||
        input=="trackin4"||
        input=="trackin5"||
        input=="trackin6"||
        input=="trackin7"||
        input=="trackin8"||
        input=="trackin9"||
        input=="trackin10"){
        FPGATrackFit* tmp=dynamic_cast<FPGATrackFit*>(memory);
        assert(tmp!=0);
        inputtracklets_.push_back(tmp);
        return;  
    }
    cout << "Did not find input : "<<input<<endl;
    assert(0);
  }

  void execute(std::vector<FPGATrack*>& outputtracks_) {

    inputtracks_.clear();
    
    //cout << "In FPGAPurgeDuplicate::execute() " << iSector_ << endl;
    //cout << "Number of input track lists  : " << inputtracklets_.size() << endl;
    for (unsigned int i=0;i<inputtracklets_.size();i++) {
      if(inputtracklets_[i]->nTracks()==0) continue;
      //  cout << "Input tracklist " << i << " " << inputtracklets_[i]->getName()
      //<< " with " << inputtracklets_[i]->nTracks() << endl;
      for(unsigned int j=0;j<inputtracklets_[i]->nTracks();j++){
        FPGATrack* aTrack=inputtracklets_[i]->getTrack(j)->getTrack();
        inputtracks_.push_back(aTrack);
    //    cout << "FPGATrack with duplicate = "<<aTrack->duplicate()<<endl;
      }
    }
    if(inputtracks_.size()==0) return;
    
    unsigned int numTrk = inputtracks_.size();
    //cout << "Sector = " << iSector_ << " Number of tracks = " << numTrk << endl;

    // Set the sector for FPGATrack, enabling the ability to do adjacent sector removal
    for(unsigned int strk=0; strk<numTrk; strk++) {
      inputtracks_[strk]->setSector(iSector_);

    }
    for(unsigned int itrk=0; itrk<numTrk-1; itrk++) { // numTrk-1 since last track has no other to compare to

      // If primary track is a duplicate, it cannot veto any...move on
      if(inputtracks_[itrk]->duplicate()==1) continue;

      int nStubP = 0;
      int nStubS[300] = {0};
      int nShare[300] = {0};

      // Get and count primary stubs
      std::map<int, int> stubsTrk1 = inputtracks_[itrk]->stubID();
      nStubP = stubsTrk1.size();

      for(unsigned int jtrk=itrk+1; jtrk<numTrk; jtrk++) {
        // Skip duplicate tracks
        if(inputtracks_[jtrk]->duplicate()==1) continue;

        // Get and count secondary stubs
        std::map<int, int> stubsTrk2 = inputtracks_[jtrk]->stubID();
        nStubS[jtrk] = stubsTrk2.size();

        // Count shared stubs
        for(std::map<int, int>::iterator  st=stubsTrk1.begin(); st!=stubsTrk1.end(); st++) {
          if(stubsTrk2.find(st->first) != stubsTrk2.end()) {
            if(st->second == stubsTrk2[st->first] && st->second != 63) nShare[jtrk]++;
          }
        }
      }

      // Tag duplicates
      for(unsigned int jtrk=itrk+1; jtrk<numTrk; jtrk++) {
        // Skip duplicate tracks
        if(inputtracks_[jtrk]->duplicate()==1) continue;
        if((nStubP-nShare[jtrk] < minIndStubs) && (nStubP <  nStubS[jtrk])) {
          inputtracks_[itrk]->setDuplicate(true);
        }
        if((nStubS[jtrk]-nShare[jtrk] < minIndStubs) && (nStubS[jtrk] <= nStubP)) {
          inputtracks_[jtrk]->setDuplicate(true);
        }
      }

    } //loop over primary track

    //Add tracks to output
    for(unsigned int i=0;i<inputtracklets_.size();i++) {
      for(unsigned int j=0;j<inputtracklets_[i]->nTracks();j++) {
        //cout << inputtracklets_[i]->getTrack(j)->getTrack()->duplicate() << " " << inputtracklets_[i]->getTrack(j)->trackfitstr() << endl;
        if(inputtracklets_[i]->getTrack(j)->getTrack()->duplicate()==0) {
          outputtracklets_[i]->addTrack(inputtracklets_[i]->getTrack(j));
        }
        //outputtracklets_[i]->addTrack(inputtracklets_[i]->getTrack(j));
        //For root file:
        outputtracks_.push_back(inputtracklets_[i]->getTrack(j)->getTrack());
      }
    }
  }

private:

  std::vector<FPGATrack*> inputtracks_;
  std::vector<FPGATrackFit*> inputtracklets_;
  std::vector<FPGACleanTrack*> outputtracklets_;
  
};


#endif

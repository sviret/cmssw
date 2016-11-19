//This class implementes the tracklet engine
#ifndef FPGAMATCHENGINE_H
#define FPGAMATCHENGINE_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAMatchEngine:public FPGAProcessBase{

public:

  FPGAMatchEngine(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    layer_=0;
    disk_=0;
    string subname=name.substr(8,2);
    if (subname=="L1") layer_=1;
    if (subname=="L2") layer_=2;
    if (subname=="L3") layer_=3;
    if (subname=="L4") layer_=4;
    if (subname=="L5") layer_=5;
    if (subname=="L6") layer_=6;
    if (subname=="F1") disk_=1;
    if (subname=="F2") disk_=2;
    if (subname=="F3") disk_=3;
    if (subname=="F4") disk_=4;
    if (subname=="F5") disk_=5;
    if (subname=="B1") disk_=-1;
    if (subname=="B2") disk_=-2;
    if (subname=="B3") disk_=-3;
    if (subname=="B4") disk_=-4;
    if (subname=="B5") disk_=-5;
    if (layer_==0&&disk_==0) {
      cout << name<<" subname = "<<subname<<" "<<layer_<<" "<<disk_<<endl;
    }
    assert((layer_!=0)||(disk_!=0));

  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="matchout") {
      FPGACandidateMatch* tmp=dynamic_cast<FPGACandidateMatch*>(memory);
      assert(tmp!=0);
      candmatches_=tmp;
      return;
    }
    assert(0);

  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="vmstubin") {
      FPGAVMStubs* tmp=dynamic_cast<FPGAVMStubs*>(memory);
      assert(tmp!=0);
      vmstubs_=tmp;
      return;
    }
    if (input=="vmprojin") {
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      vmprojs_=tmp;
      return;
    }
    cout << "Could not find input : "<<input<<endl;
    assert(0);
  }

  void execute() {

    unsigned int countall=0;
    unsigned int countpass=0;

    
    for(unsigned int j=0;j<vmprojs_->nTracklets();j++){

      FPGATracklet* proj=vmprojs_->getFPGATracklet(j);
      for(unsigned int i=0;i<vmstubs_->nStubs();i++){
	std::pair<FPGAStub*,L1TStub*> stub=vmstubs_->getStub(i);
	countall++;
	if (layer_>0){
	  //cout << "Stub: layer = "<<layer_<<" "<<stub.first->phivm().value()<<" "
	  //	     <<stub.first->zvm().value()<<endl;
	}
	//cout << "Adding match in "<<getName()<<endl;
	if (layer_>0) {
	  //cout << "  Proj: "<<proj->phiprojvm(layer_)
	  //     <<" "<<proj->zprojvm(layer_)<<" "<<proj->zproj(layer_)<<endl;
	  int iphivmdiff=stub.first->phivm().value()-proj->phiprojvm(layer_);
	  int izvmdiff=stub.first->zvm().value()-proj->zprojvm(layer_);
	  //cout << "ME diff : "<<layer_<<" "<<iphivmdiff<<" "<<izvmdiff<<endl;
	  if (doMEMatch){
	    if (abs(iphivmdiff)>1) continue;
	    if (abs(izvmdiff)>1000) continue; //dummy for now
	  }
	}
	if (disk_>0) {
	  int iphivmdiff=stub.first->phivm().value()-proj->phiprojvmdisk(disk_);
	  int irvmdiff=stub.first->rvm().value()-proj->rprojvmdisk(disk_);
	  //cout << getName()<<" "<<disk_<<" "<<iphivmdiff<<" "<<irvmdiff<<endl;
	  if (doMEMatch){
	    if (abs(iphivmdiff)>1) continue;
	    if (abs(irvmdiff)>1000) continue; //dummy for now
	  }
	  //cout << "FPGAMatchEngine::execute found candidate match in disk :"
	  //     <<disk_<<endl;
	}
	//cout << "FPGAMatchEngine "<<getName()<<" stub.r "
	//     <<stub.second->r()<<" "<<vmstubs_->getName()<<endl;
	candmatches_->addMatch(proj,stub);
	countpass++;
	if (countall>=MAXME) break;
      }
      if (countall>=MAXME) break;
    }

    if (writeME) {
      static ofstream out("matchengine.txt");
      out << getName()<<" "<<countall<<" "<<countpass<<endl;
    }

  }


private:

  FPGAVMStubs* vmstubs_;
  FPGAVMProjections* vmprojs_;

  FPGACandidateMatch* candmatches_;

  int layer_;
  int disk_;
 
};

#endif

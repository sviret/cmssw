//This class implementes the projection tranceiver
#ifndef FPGAMATCHTRANSCEIVER_H
#define FPGAMATCHTRANSCEIVER_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAMatchTransceiver:public FPGAProcessBase{

public:

  FPGAMatchTransceiver(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="matchout1"||
	output=="matchout2"||
	output=="matchout3"||
	output=="matchout4"||
	output=="matchout5"||
	output=="matchout6"||
	output=="matchout7"||
	output=="matchout8"||
	output=="matchout9"||
	output=="matchout10"||
	output=="matchout11"||
	output=="matchout12"||
	output=="matchout13"||
	output=="matchout14"||
	output=="matchout15"||
	output=="matchout16"||
	output=="matchout17"||
	output=="matchout18"||
	output=="matchout19"||
	output=="matchout20"||
	output=="matchout21"
	){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      outputmatches_.push_back(tmp);
      return;
    }
    cout << "In FPGAMatchTransceiver: Did not find output  = "<<output<<endl;
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="matchin1"||
	input=="matchin2"||
	input=="matchin3"||
	input=="matchin4"||
	input=="matchin5"||
	input=="matchin6"||
	input=="matchin7"||
	input=="matchin8"||
	input=="matchin9"||
	input=="matchin10"||
	input=="matchin11"||
	input=="matchin12"||
	input=="matchin13"||
	input=="matchin14"||
	input=="matchin15"||
	input=="matchin16"||
	input=="matchin17"||
	input=="matchin18"||
	input=="matchin19"||
	input=="matchin20"||
	input=="matchin21"||
	input=="matchin22"||
	input=="matchin23"||
	input=="matchin24"||
	input=="matchin25"||
	input=="matchin26"||
	input=="matchin27"||
	input=="matchin28"||
	input=="matchin29"||
	input=="matchin30"||
	input=="matchin31"||
	input=="matchin32"||
	input=="matchin33"||
	input=="matchin34"||
	input=="matchin35"||
	input=="matchin36"||
	input=="matchin37"||
	input=="matchin38"||
	input=="matchin39"||
	input=="matchin40"||
	input=="matchin41"||
	input=="matchin42"||
	input=="matchin43"||
	input=="matchin44"||
	input=="matchin45"||
	input=="matchin46"||
	input=="matchin47"||
	input=="matchin48"||
	input=="matchin49"||
	input=="matchin50"
	){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      inputmatches_.push_back(tmp);
      return;
    }
    
    cout << "In FPGAMatchTransceiver: Did not find input  = "<<input<<endl;
    assert(0);
  }

  //this->inputmatches_ to other->outputmatches_ 
  void execute(FPGAMatchTransceiver* other){
    bool write=false; 
    if (write) cout << "MT name = "<<getName()<<endl;
    int count=0;
    for(unsigned int i=0;i<inputmatches_.size();i++){
      string basename=inputmatches_[i]->getName().substr(0,10);
      if (write) cout << "input : "<<inputmatches_[i]->getName()
      		      << " "<<basename<<" #matches = "
		      <<inputmatches_[i]->nMatches()<<endl;
      bool wrote=false;
      for(unsigned int j=0;j<other->outputmatches_.size();j++){
	if (write) cout << "output  : "<<other->outputmatches_[j]->getName()<<endl;
	std::size_t found = other->outputmatches_[j]->getName().find(basename);
	if (found!=std::string::npos){
	  if (write) cout << "Will write to "
			  <<other->outputmatches_[j]->getName()<<endl;
	  wrote=true;
	  for(unsigned int l=0;l<inputmatches_[i]->nMatches();l++){
	    //cout << getName()<<" from = "<<inputmatches_[i]->getName()
	    //	 << " to = "<<other->outputmatches_[j]->getName()
	    // 	 << " "<<iSector_
	    // 	 << " "<<inputmatches_[i]->getMatch(l).first<<endl;
	    other->outputmatches_[j]->addMatch(inputmatches_[i]->getMatch(l));
	  }
	  count+=inputmatches_[i]->nMatches();
	  //continue;
	}	
      }
      assert(wrote);
    }
    if (writeMatchTransceiver) {
      static ofstream out("matchtransceiver.txt");
      out << getName() << " " 
	  << count << endl;
    }
  }
  

private:

  vector<FPGAFullMatch*> inputmatches_;

  vector<FPGAFullMatch*> outputmatches_;

};

#endif

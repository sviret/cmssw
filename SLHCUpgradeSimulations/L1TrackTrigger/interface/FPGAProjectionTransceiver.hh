//This class implementes the projection tranceiver
#ifndef FPGAPROJECTIONTRANSCEIVER_H
#define FPGAPROJECTIONTRANSCEIVER_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAProjectionTransceiver:public FPGAProcessBase{

public:

  FPGAProjectionTransceiver(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="projout"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojections_.push_back(tmp);
      return;
    }
    /*
    if (output=="output_L1L2_2"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory
      assert(tmp!=0);
      outputprojectionsL1L2_2_=tmp;
      return;
    }
    if (output=="output_L1L2_3"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL1L2_3_=tmp;
      return;
    }
    if (output=="output_L1L2_4"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL1L2_4_=tmp;
      return;
    }
    if (output=="output_L3L4_1"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL3L4_1_=tmp;
      return;
    }
    if (output=="output_L3L4_2"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL3L4_2_=tmp;
      return;
    }
    if (output=="output_L3L4_3"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL3L4_3_=tmp;
      return;
    }
    if (output=="output_L3L4_4"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL3L4_4_=tmp;
      return;
    }
    if (output=="output_L5L6_1"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL5L6_1_=tmp;
      return;
    }
    if (output=="output_L5L6_2"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL5L6_2_=tmp;
      return;
    }
    if (output=="output_L5L6_3"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL5L6_3_=tmp;
      return;
    }
    if (output=="output_L5L6_4"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      outputprojectionsL5L6_4_=tmp;
      return;
    }
    */
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }

    /*
    if (input=="input_L1L2_1"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL1L2_1_=tmp;
      return;
    }
    if (input=="input_L1L2_2"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL1L2_2_=tmp;
      return;
    }
    if (input=="input_L1L2_3"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL1L2_3_=tmp;
      return;
    }
    if (input=="input_L1L2_4"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL1L2_4_=tmp;
      return;
    }
    if (input=="input_L3L4_1"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL3L4_1_=tmp;
      return;
    }
    if (input=="input_L3L4_2"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL3L4_2_=tmp;
      return;
    }
    if (input=="input_L3L4_3"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL3L4_3_=tmp;
      return;
    }
    if (input=="input_L3L4_4"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL3L4_4_=tmp;
      return;
    }
    if (input=="input_L5L6_1"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL5L6_1_=tmp;
      return;
    }
    if (input=="input_L5L6_2"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL5L6_2_=tmp;
      return;
    }
    if (input=="input_L5L6_3"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL5L6_3_=tmp;
      return;
    }
    if (input=="input_L5L6_4"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojectionsL5L6_4_=tmp;
      return;
    }
    */
    if (input=="projin"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputprojections_.push_back(tmp);
      return;
    }

    assert(0);
  }

  //Copy otherSector->inputprojections_ to this->outputprojections_ 
  void execute(FPGAProjectionTransceiver* otherSector){

    if (!doProjections) return;

    unsigned int count=0;
    //cout << "in FPGAProjectionTransceiver "<<otherSector->inputprojections_.size()<<endl;
    for(unsigned int i=0;i<otherSector->inputprojections_.size();i++){
      FPGATrackletProjections* otherProj=otherSector->inputprojections_[i];
      for (unsigned int l=0;l<otherProj->nTracklets();l++){
	count++;
	FPGATracklet* tracklet=otherProj->getFPGATracklet(l);
	FPGAWord fpgaphi;
        string seedpair=otherProj->getSeedPair();
        string target=otherProj->getTarget();
	//cout << "FPGAProjectionTransceiver "<<getName()<<" "
	//     << otherProj->getName()<<"  = "<<otherProj->getSeedPair()
	//     << " : "<<otherProj->getTarget()
	//     << " tracklet layer disk : "<<tracklet->layer()
	//     << " "<<tracklet->disk()<<endl;

	if (seedpair=="F1L1"&&target=="F5") seedpair="F1F2";
	if (seedpair=="F1L1"&&target=="F4") seedpair="F1F2";
	if (seedpair=="F1L1"&&target=="F3") seedpair="F1F2";
	if (seedpair=="F1L1"&&target=="F2") seedpair="F3F4";

	if (seedpair=="B1L1"&&target=="B5") seedpair="B1B2";
	if (seedpair=="B1L1"&&target=="B4") seedpair="B1B2";
	if (seedpair=="B1L1"&&target=="B3") seedpair="B1B2";
	if (seedpair=="B1L1"&&target=="B2") seedpair="B3B4";

	if (seedpair=="L1L2"&&target=="F5") seedpair="F3F4";
	if (seedpair=="L1L2"&&target=="F4") seedpair="F1F2";
	if (seedpair=="L1L2"&&target=="F3") seedpair="F1F2";
	if (seedpair=="L1L2"&&target=="F2") seedpair="F3F4";
	if (seedpair=="L3L4"&&target=="F2") seedpair="F3F4";
	if (seedpair=="L1L2"&&target=="F1") seedpair="F3F4";
	if (seedpair=="L3L4"&&target=="F1") seedpair="F3F4";

	if (seedpair=="L1L2"&&target=="B5") seedpair="B3B4";
	if (seedpair=="L1L2"&&target=="B4") seedpair="B1B2";
	if (seedpair=="L1L2"&&target=="B3") seedpair="B1B2";
	if (seedpair=="L1L2"&&target=="B2") seedpair="B3B4";
	if (seedpair=="L3L4"&&target=="B2") seedpair="B3B4";
	if (seedpair=="L1L2"&&target=="B1") seedpair="B3B4";
	if (seedpair=="L3L4"&&target=="B1") seedpair="B3B4";

	if (seedpair=="F1F2"&&target=="L1") seedpair="L5L6";
	if (seedpair=="F1F2"&&target=="L2") seedpair="L5L6";
	if (seedpair=="B1B2"&&target=="L1") seedpair="L5L6";
	if (seedpair=="B1B2"&&target=="L2") seedpair="L5L6";

	string region="";

	if (target[0]=='F'||target[0]=='B') {
	  int disk=target[1]-'0';
	  if (target[0]=='B') disk=-disk;
	  FPGAWord rproj=tracklet->fpgarprojdisk(disk);
	  //FIXME should not use the floats...
	  int ir=2*(rproj.value()*krprojshiftdisk)/(rmaxdisk-rmindisk)+1;
	  if (target[0]=='F'&&ir==1) region="D5";
	  if (target[0]=='F'&&ir==2) region="D6";
	  if (target[0]=='B'&&ir==1) region="D7";
	  if (target[0]=='B'&&ir==2) region="D8";

	  //cout << "disk = "<<disk<<" "<<ir<<" "<<region<<endl;
	} 


	if (target[0]=='L') {
	  int layer=target[1]-'0';
	  FPGAWord zproj=tracklet->fpgazproj(layer);
	  int iz=4+(zproj.value()>>(zproj.nbits()-3));
	  iz=iz/2+1;
	  if (iz==1) region="D1";
	  if (iz==2) region="D2";
	  if (iz==3) region="D3";
	  if (iz==4) region="D4";

	  //cout << "layer = "<<layer<<" "<<iz<<" "<<region<<endl;
	} 


	assert(region!="");

	int foundTarget=0;
	
	for(unsigned int ii=0;ii<outputprojections_.size();ii++){
	  FPGATrackletProjections* projs=outputprojections_[ii];
	  string name=projs->getName();
	  string seeding=name.substr(name.size()-4,4);
	  string toTarget=name.substr(name.size()-9,2);
	  string toRegion=name.substr(name.size()-7,2);
	  //cout << "  target "<<name<<" "<<seeding<<" "<<toTarget
	  //     << " "<<toRegion<<endl;
	  if (seeding!=seedpair) continue;
	  if (target!=toTarget) continue;
	  if (region!=toRegion) continue;
	  count++;
          if (count>MAXPROJECTIONTRANSCEIVER) continue;
	  projs->addTracklet(tracklet);
	  foundTarget++;
	}

	if (foundTarget==0) {
	  cout << "ERROR in FPGAProjectionTransceiver "<<getName()
	       <<" no target for :"<< otherProj->getName() 
	       <<" seedpair "<<seedpair<<endl;
	}

	assert(foundTarget<2);

	/*
	int disk1=disk_;
	if (tracklet->t()<0.0) disk1=-disk_;
	//Handle PT that handles both disk and layer
	if (layer_!=0&&disk1!=0) {
	  if (abs(tracklet->disk())){
	    if (layer_<3){
	      fpgaphi=tracklet->fpgaphiproj(layer_);
	      layer=true;
	    }
	    else {
	      fpgaphi=tracklet->fpgaphiprojdisk(disk1);
	      disk=true;
	    }
	  } else {
	    //cout << "FPGAProjectionTransceiver "<<tracklet<<" layer_ = "<<layer_<<endl;
	    if (!tracklet->fpgazproj(layer_).atExtreme()){
	      fpgaphi=tracklet->fpgaphiproj(layer_);
	      layer=true;
	    } else {
	      fpgaphi=tracklet->fpgaphiprojdisk(disk1);
	      if (!tracklet->fpgarprojdisk(disk1).atExtreme()){
		disk=true;
	      }
	    }
	  }
	  if (fpgaphi.atExtreme()) {
	    cout << "FPGAProjectionTransceiver: Warning skipping projection"<<endl;
	    continue;
	  }
	}
	//Handle PT to only a layer
	if (layer_!=0&&disk1==0) {
	  fpgaphi=tracklet->fpgaphiproj(layer_);
	  layer=true;
	}
	//Handle PT to only a disk
	if (disk1!=0&&layer_==0) {
	  fpgaphi=tracklet->fpgaphiprojdisk(disk1);
	  disk=true;
	}
	
	assert(disk1!=0||layer_!=0);
      

	assert(disk||layer);

	//cout << "(bool) layer disk : "<<layer<<" "<<disk<<endl;

	*/

	/*	

	int iphivmRaw=fpgaphi.value()>>(fpgaphi.nbits()-5);
    
	//cout << "FPGAProjectionTransceiver iphivmRaw "<<iphivmRaw << " "
	//     <<((fpgaphi.value()+1)>>(fpgaphi.nbits()-5))<<endl;
	if (iphivmRaw<4||iphivmRaw>27) {
	  cout << "FPGAProjectionTransceiver "<<getName()<<" will skip projection"<<endl;
	  continue;
	}
	assert(iphivmRaw>=4);
	assert(iphivmRaw<=27);

	int iphi=(iphivmRaw-4)>>3;

	//cout << "FPGAProjectionTranceiver "<<getName()<<" layer fpgaphi iphivmRaw iphi : "<<layer_<<" "<<fpgaphi.value()<<" "<<iphivmRaw<<" "<<iphi<<endl;

    
	assert(iphi>=0);
	assert(iphi<=2);

	if (iphi==0) {
	  if (layer) {
	    assert(outputprojLPHI1!=0);
	    //cout << "Adding tracklet "<<otherProj->getFPGATracklet(l)<<" to "<<outputprojLPHI1->getName()<<endl;
	    outputprojLPHI1->addTracklet(otherProj->getFPGATracklet(l));
	  }
	  if (disk) {
	    if (outputprojDPHI1==0) {
	      cout << "FPGAProjectionTransceiver in : "<<getName()<< " outputprojDPHI1 is zero"<<endl;
	    }
	    assert(outputprojDPHI1!=0);
	    //cout << "Adding tracklet "<<otherProj->getFPGATracklet(l)<<" to "<<outputprojDPHI1->getName()<<endl;
	    //cout << "FPGAProjectionTransceiver add projection to : "<<outputprojDPHI1->getName()<<endl;
	    outputprojDPHI1->addTracklet(otherProj->getFPGATracklet(l));
	  }
	}

	if (iphi==1) {
	  if (layer) {
	    if (outputprojLPHI2==0) {
	      cout << "FPGAProjectionTransceiver in : "<<getName()<< " outputprojLPHI2 is zero"<<endl;
	    }
	    assert(outputprojLPHI2!=0);
	    //cout << "Adding tracklet "<<otherProj->getFPGATracklet(l)<<" to "<<outputprojLPHI2->getName()<<endl;
	    outputprojLPHI2->addTracklet(otherProj->getFPGATracklet(l));
	  }
	  if (disk) {
	    if (outputprojDPHI2==0) {
	      cout << "FPGAProjectionTransceiver in : "<<getName()<< " outputprojDPHI2 is zero"<<endl;
	    }
	    assert(outputprojDPHI2!=0);
	    //cout << "Adding tracklet "<<otherProj->getFPGATracklet(l)<<" to "<<outputprojDPHI2->getName()<<endl;
	    outputprojDPHI2->addTracklet(otherProj->getFPGATracklet(l));
	  }
	}
	
	if (iphi==2) {
	  if (layer) {
	    //cout << "In getName = "<<getName()<<endl;
	    assert(outputprojLPHI3!=0);
	    //cout << "Adding tracklet "<<otherProj->getFPGATracklet(l)<<" to "<<outputprojLPHI3->getName()<<endl;
	    outputprojLPHI3->addTracklet(otherProj->getFPGATracklet(l));
	  }
	  if (disk) {
	    if (outputprojDPHI3==0) {
	      cout << "FPGAProjectionTransceiver in : "<<getName()<< " outputprojDPHI3 is zero"<<endl;
	    }
	    assert(outputprojDPHI3!=0);
	    //cout << "FPGAProjectionTransceiver add projection to : "<<outputprojDPHI3->getName()<<endl;
	    //cout << "Adding tracklet "<<otherProj->getFPGATracklet(l)<<" to "<<outputprojDPHI3->getName()<<endl;
	    outputprojDPHI3->addTracklet(otherProj->getFPGATracklet(l));
	  }
	}
	*/

      }

    }



    if (writeProjectionTransceiver) {
      static ofstream out("projectiontransceiver.txt");
      out << getName() << " " 
	  << count << endl;
    }

  }


private:

  vector<FPGATrackletProjections*> inputprojections_;

  vector<FPGATrackletProjections*> outputprojections_;

};

#endif

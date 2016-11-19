//This class implementes the projection router
#ifndef FPGAPROJECTIONROUTER_H
#define FPGAPROJECTIONROUTER_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAProjectionRouter:public FPGAProcessBase{

public:

  FPGAProjectionRouter(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    string subname=name.substr(0,5);
    string subname1=name.substr(0,6);
    string subname2=name.substr(0,3);
    //cout << "subname : "<<subname<<endl;
    layer_=0;
    disk_=0;
    vmprojPHI1R1_=0;
    vmprojPHI1R2_=0;
    vmprojPHI2R1_=0;
    vmprojPHI2R2_=0;
    vmprojPHI3R1_=0;
    vmprojPHI3R2_=0;
    vmprojPHI4R1_=0;
    vmprojPHI4R2_=0;

    dct_=name[6]-'0';
    if (subname2=="PRD") {
      dct_=name[7]-'0';
    }

    if (subname=="PR_L1") layer_=1;
    if (subname=="PR_L2") layer_=2;
    if (subname=="PR_L3") layer_=3;
    if (subname=="PR_L4") layer_=4;
    if (subname=="PR_L5") layer_=5;
    if (subname=="PR_L6") layer_=6;
    if (subname1=="PRD_F1") disk_=1;
    if (subname1=="PRD_F2") disk_=2;
    if (subname1=="PRD_F3") disk_=3;
    if (subname1=="PRD_F4") disk_=4;
    if (subname1=="PRD_F5") disk_=5;
    if (subname1=="PRD_B1") disk_=-1;
    if (subname1=="PRD_B2") disk_=-2;
    if (subname1=="PRD_B3") disk_=-3;
    if (subname1=="PRD_B4") disk_=-4;
    if (subname1=="PRD_B5") disk_=-5;
    allproj_=0;
    assert(disk_!=0||layer_!=0);
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="allprojout"){
      FPGAAllProjections* tmp=dynamic_cast<FPGAAllProjections*>(memory);
      assert(tmp!=0);
      allproj_=tmp;
      return;
    }
    if (output=="vmprojoutPHI1X1"){
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmprojPHI1Z1_=tmp;
      }
      if (disk_!=0) {
	vmprojPHI1R1_=tmp;
      }
      return;
    }
    if (output=="vmprojoutPHI1X2"){
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmprojPHI1Z2_=tmp;
      }
      if (disk_!=0) {
	vmprojPHI1R2_=tmp;
      }
      return;
    }
    if (output=="vmprojoutPHI2X1"){
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmprojPHI2Z1_=tmp;
      }
      if (disk_!=0) {
	vmprojPHI2R1_=tmp;
      }
      return;
    }
    if (output=="vmprojoutPHI2X2"){
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmprojPHI2Z2_=tmp;
      }
      if (disk_!=0) {
	vmprojPHI2R2_=tmp;
      }
      return;
    }
    if (output=="vmprojoutPHI3X1"){
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmprojPHI3Z1_=tmp;
      }
      if (disk_!=0) {
	vmprojPHI3R1_=tmp;
      }
      return;
    }
    if (output=="vmprojoutPHI3X2"){
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmprojPHI3Z2_=tmp;
      }
      if (disk_!=0) {
	vmprojPHI3R2_=tmp;
      }
      return;
    }
    if (output=="vmprojoutPHI4X1"){
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmprojPHI4Z1_=tmp;
      }
      if (disk_!=0) {
	vmprojPHI4R1_=tmp;
      }
      return;
    }
    if (output=="vmprojoutPHI4X2"){
      FPGAVMProjections* tmp=dynamic_cast<FPGAVMProjections*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmprojPHI4Z2_=tmp;
      }
      if (disk_!=0) {
	vmprojPHI4R2_=tmp;
      }
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
    if (input=="proj1in"||input=="proj2in"||
	input=="proj3in"||input=="proj4in"||
	input=="proj5in"||input=="proj6in"||
	input=="proj7in"||input=="proj8in"||
	input=="proj9in"||input=="proj10in"||
	input=="proj11in"||input=="proj12in"||
	input=="proj13in"||input=="proj14in"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputproj_.push_back(tmp);
      return;
    }
    if (input=="projplusin"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputplusproj_=tmp;
      return;
    }
    if (input=="projminusin"){
      FPGATrackletProjections* tmp=dynamic_cast<FPGATrackletProjections*>(memory);
      assert(tmp!=0);
      inputminusproj_=tmp;
      return;
    }
    cout << "Could not find input : "<<input<<endl;
    assert(0);
  }

  void execute() {

    //cout << "FPGAProjectionRouter : "<<getName()<<endl;

    unsigned int count=0;

    if (layer_!=0) {
      for (unsigned int j=0;j<inputproj_.size();j++){
	//cout << "Inputproj : "<<inputproj_[j]->getName()<<" "
        //    <<inputproj_[j]->nTracklets()<<endl;
	for (unsigned int i=0;i<inputproj_[j]->nTracklets();i++){
	  count++;
	  if (count>MAXPROJROUTER) continue;
	  //cout << "Doing projection"<<endl;
	  
	  FPGAWord fpgaphi=inputproj_[j]->getFPGATracklet(i)->fpgaphiproj(layer_);
	  FPGAWord fpgaz=inputproj_[j]->getFPGATracklet(i)->fpgazproj(layer_);

	  //if (inputproj_[j]->getFPGATracklet(i)->plusNeighbor(layer_)) {
	  //  cout << "Found plus neighbor in : "<<inputproj_[j]->getName()<<endl;
	  //} 
	  
	  //skip if projection is out of range!
	  if (fpgaz.atExtreme()) continue;
	  if (fpgaphi.atExtreme()) continue;


	  //	  cout<<" string: "<<inputproj_[j]->getFPGATracklet(i)->trackletprojstr(i)<<endl; 

	
	  int iz=4+(fpgaz.value()>>(fpgaz.nbits()-3));
	  int iphitmp=fpgaphi.value();
	  if ((layer_%2)==1) iphitmp-=(1<<(fpgaphi.nbits()-3));  
	  int iphi=iphitmp>>(fpgaphi.nbits()-2);

	  int iztmp=iz/2+1;
	  
	  //cout << "FPGAProjectionRouter : "<<getName()
	  //     <<" iz = "<<iz<<" iztmp = "<<iztmp
	  //     <<" dct = "<<dct_<<endl;
	  
	  iz=iz%2;

	  if (dct_!=iztmp) iz=1-iz;  //test hack
	  
	  //if (layer_==3) {
	  //  cout << "Will add to allproj_ in "<<getName()
	  //	 <<" z = "<<fpgaz.value()*kzproj
	  //		 <<" from :"<<inputproj_[j]->getName()
	  //		 <<endl;
	  //}

	  assert(allproj_!=0);

	  unsigned int index=allproj_->nTracklets();

	  allproj_->addTracklet(inputproj_[j]->getFPGATracklet(i));
	  
	  if (iphi==0&&iz==0) {
	    vmprojPHI1Z1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  if (iphi==0&&iz==1) {
	    vmprojPHI1Z2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  
	  if (iphi==1&&iz==0) {
	    vmprojPHI2Z1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  if (iphi==1&&iz==1) {
	    vmprojPHI2Z2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  
	  if (iphi==2&&iz==0) {
	    vmprojPHI3Z1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  if (iphi==2&&iz==1) {
	    vmprojPHI3Z2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	
	  if (iphi==3&&iz==0) {
	    vmprojPHI4Z1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  if (iphi==3&&iz==1) {
	    vmprojPHI4Z2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	}
      }
    } else {
      for (unsigned int j=0;j<inputproj_.size();j++){
	//cout << "Inputproj : "<<inputproj_[j]->getName()<<" "
        //     <<inputproj_[j]->nTracklets()<<endl;
	for (unsigned int i=0;i<inputproj_[j]-> nTracklets();i++){
	  count++;
	  if (count>MAXPROJROUTER) continue;

	  FPGAWord fpgaphi=inputproj_[j]->getFPGATracklet(i)->fpgaphiprojdisk(disk_);
	  FPGAWord fpgar=inputproj_[j]->getFPGATracklet(i)->fpgarprojdisk(disk_);

	  //cout << "fpgaphi fpgar : "<<fpgaphi.nbits()<<" "<<fpgar.nbits()<<endl;
	 
	  // cout<<" string: "<<inputproj_[j]->getFPGATracklet(i)->trackletprojstrD(i)<<endl; 

	  
	  //skip if projection is out of range!
	  if (fpgar.atExtreme()) continue;
	  if (fpgaphi.atExtreme()) continue;
	
	  // Print out Projected Tracklets	  


	  //FIXME this should come from the integer calculations!!!
	  int ir=(1<<2)*(fpgar.value()*krprojshiftdisk)/(rmaxdisk-rmindisk);
          

	  int iphitmp=fpgaphi.value();
	  if ((disk_%2)==0) iphitmp-=(1<<(fpgaphi.nbits()-3));  
	  int iphi=iphitmp>>(fpgaphi.nbits()-2);
	  
	  //cout << "iphi ir : "<<iphi<<" "<<ir<<endl;
	  
	  ir=ir%2;
	  
	  //cout << "Will add to allproj_ in "<<getName()<<endl;
	  assert(allproj_!=0);

	  unsigned int index=allproj_->nTracklets();
	  allproj_->addTracklet(inputproj_[j]->getFPGATracklet(i));
	  
	  if (iphi==0&&ir==0) {
	    //cout << "Adding projection to "<<vmprojPHI1R1_->getName()<<endl;
	    vmprojPHI1R1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  if (iphi==0&&ir==1) {
	    //cout << "Adding projection to "<<vmprojPHI1R2_->getName()<<endl;
	    vmprojPHI1R2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  
	  if (iphi==1&&ir==0) {
	    //cout << "Adding projection to "<<vmprojPHI2R1_->getName()<<endl;
	    vmprojPHI2R1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  if (iphi==1&&ir==1) {
	    //cout << "Adding projection to "<<vmprojPHI2R2_->getName()<<endl;
	    vmprojPHI2R2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  
	  if (iphi==2&&ir==0) {
	    //cout << "Adding projection to "<<vmprojPHI3R1_->getName()<<endl;
	    vmprojPHI3R1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  if (iphi==2&&ir==1) {
	    //cout << "Adding projection to "<<vmprojPHI3R2_->getName()<<endl;
	    vmprojPHI3R2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	
	  if (iphi==3&&ir==0&&vmprojPHI4R1_!=0) {
	    //cout << "Adding projection to "<<vmprojPHI4R1_->getName()<<endl;
	    vmprojPHI4R1_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	  if (iphi==3&&ir==1&&vmprojPHI4R2_!=0) {
	    //cout << "Adding projection to "<<vmprojPHI4R2_->getName()<<endl;
	    vmprojPHI4R2_->addTracklet(inputproj_[j]->getFPGATracklet(i),index);
	  }
	}
      }
    }

     
    if (writeAllProjections) {
      static ofstream out("allprojections.txt");
      out << getName()<<" "<<allproj_->nTracklets()<<endl;
      /*
      if( allproj_->nTracklets() > 0 && disk_ == 0  && allproj_->nTracklets() < 4 ){
	// 	out << getName() << " " << allproj_->nTracklets() << endl;
      	for(uint i  = 0; i  < allproj_->nTracklets(); i++){
	  out <<i<<" "<<getName()<<" "<< allproj_->getFPGATracklet(i)->trackletprojstr(i)<<" "<<layer_<<endl;
	}
      }
      if( allproj_->nTracklets() > 0 && layer_ == 0  && allproj_->nTracklets() < 4 ){
	out << getName() << " " << allproj_->nTracklets() << endl;
      	for(uint i  = 0; i  < allproj_->nTracklets(); i++){
	  out <<i<<" "<< allproj_->getFPGATracklet(i)->trackletprojstrD(i)<<" "<<layer_<<endl;
	}
      } 
      */
    }

    if (layer_!=0) {
      if (writeVMProjections) {
	static ofstream out("vmprojections.txt");

	/*
	if(vmprojPHI1Z1_->nTracklets() != 0  ){
	  out << vmprojPHI1Z1_->getName() << " " << vmprojPHI1Z1_->nTracklets() << endl; 
	  for(uint i = 0; i < vmprojPHI1Z1_->nTracklets(); i++){
	    out<<i<<" "<< vmprojPHI1Z1_->getFPGATracklet(i)->trackletprojstr(i)<<endl;;
	  }
	}
	if(vmprojPHI1Z2_->nTracklets() != 0  ){
	  out << vmprojPHI1Z2_->getName() << " " << vmprojPHI1Z2_->nTracklets() << endl; 
	  for(uint i = 0; i < vmprojPHI1Z2_->nTracklets(); i++){
	    out<<i<<" "<< vmprojPHI1Z2_->getFPGATracklet(i)->trackletprojstr(i)<<endl;;
	  }
	}
	if(vmprojPHI2Z1_->nTracklets() != 0  ){
	  out << vmprojPHI2Z1_->getName() << " " << vmprojPHI2Z1_->nTracklets() << endl; 
	  for(uint i = 0; i < vmprojPHI2Z1_->nTracklets(); i++){
	    out<<i<<" "<< vmprojPHI2Z1_->getFPGATracklet(i)->trackletprojstr(i)<<endl;;
	  }
	}
	if(vmprojPHI2Z2_->nTracklets() != 0  ){
	  out << vmprojPHI2Z2_->getName() << " " << vmprojPHI2Z2_->nTracklets() << endl; 
	  for(uint i = 0; i < vmprojPHI2Z2_->nTracklets(); i++){
	    out<<i<<" "<< vmprojPHI2Z2_->getFPGATracklet(i)->trackletprojstr(i)<<endl;;
	  }
	}

	cout<<"here 2"<<endl;
	if(vmprojPHI3Z1_->nTracklets() != 0 && layer_< 7  ){ 
	  out << vmprojPHI3Z1_->getName() << " " << vmprojPHI3Z1_->nTracklets() << endl; 
	  // for(uint i = 0; i < vmprojPHI3Z1_->nTracklets(); i++){
	    // out<<i<<" "<< vmprojPHI3Z1_->getFPGATracklet(i)->trackletprojstr(i)<<endl;;
	  // }
	}

	cout<<"here 3"<<endl;
	if(vmprojPHI3Z2_->nTracklets() != 0  ){
	  out << vmprojPHI3Z2_->getName() << " " << vmprojPHI3Z2_->nTracklets() << endl; 
	  //	  for(uint i = 0; i < vmprojPHI3Z2_->nTracklets(); i++){
	  //  out<<i<<" "<< vmprojPHI3Z2_->getFPGATracklet(i)->trackletprojstr(i)<<endl;;
	  //}
	}
	cout<<"here! 4"<<endl;

	if (layer_%2==0) {
	  if(vmprojPHI4Z1_->nTracklets() != 0  ) out << vmprojPHI4Z1_->getName() << " " << vmprojPHI4Z1_->nTracklets() << endl; 
	  for(uint i = 0; i < vmprojPHI4Z1_->nTracklets(); i++){
	    //    out<<i<<" "<< vmprojPHI4Z1_->getFPGATracklet(i)->trackletprojstrD(i)<<endl;;
	  }
	  
	  if(vmprojPHI4Z2_->nTracklets() != 0  ) out << vmprojPHI4Z2_->getName() << " " << vmprojPHI4Z2_->nTracklets() << endl; 
	  for(uint i = 0; i < vmprojPHI4Z2_->nTracklets(); i++){
	    // out<<i<<" "<< vmprojPHI4Z2_->getFPGATracklet(i)->trackletprojstrD(i)<<endl;;
	  }
	}
	*/
	out << vmprojPHI1Z1_->getName() << " " << vmprojPHI1Z1_->nTracklets() << endl; 
	out << vmprojPHI1Z2_->getName() << " " << vmprojPHI1Z2_->nTracklets() << endl; 
	out << vmprojPHI2Z1_->getName() << " " << vmprojPHI2Z1_->nTracklets() << endl; 
	out << vmprojPHI2Z2_->getName() << " " << vmprojPHI2Z2_->nTracklets() << endl; 
	out << vmprojPHI3Z1_->getName() << " " << vmprojPHI3Z1_->nTracklets() << endl; 
	out << vmprojPHI3Z2_->getName() << " " << vmprojPHI3Z2_->nTracklets() << endl; 
	if (layer_%2==0) {
	  out << vmprojPHI4Z1_->getName() << " " << vmprojPHI4Z1_->nTracklets() << endl; 
	  out << vmprojPHI4Z2_->getName() << " " << vmprojPHI4Z2_->nTracklets() << endl;
	}

      }
    }




  }
  

private:

  int layer_; 
  int disk_; 
  int dct_;

  vector<FPGATrackletProjections*> inputproj_;
  FPGATrackletProjections* inputplusproj_;
  FPGATrackletProjections* inputminusproj_;

  FPGAAllProjections* allproj_;
  FPGAVMProjections* vmprojPHI1Z1_;
  FPGAVMProjections* vmprojPHI1Z2_;
  FPGAVMProjections* vmprojPHI2Z1_;
  FPGAVMProjections* vmprojPHI2Z2_;
  FPGAVMProjections* vmprojPHI3Z1_;
  FPGAVMProjections* vmprojPHI3Z2_;
  FPGAVMProjections* vmprojPHI4Z1_;
  FPGAVMProjections* vmprojPHI4Z2_;

  FPGAVMProjections* vmprojPHI1R1_;
  FPGAVMProjections* vmprojPHI1R2_;
  FPGAVMProjections* vmprojPHI2R1_;
  FPGAVMProjections* vmprojPHI2R2_;
  FPGAVMProjections* vmprojPHI3R1_;
  FPGAVMProjections* vmprojPHI3R2_;
  FPGAVMProjections* vmprojPHI4R1_;
  FPGAVMProjections* vmprojPHI4R2_;


};

#endif

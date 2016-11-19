//This class implementes the VM router
#ifndef FPGAVMROUTER_H
#define FPGAVMROUTER_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAVMRouter:public FPGAProcessBase{

public:

  FPGAVMRouter(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    layer_=0;
    disk_=0;
    if (name_.substr(0,4)=="VMR_") {
      layer_=name_[6]-'0';
    }
    if (name_.substr(0,6)=="VMRD_F") {
      disk_=name_[7]-'0';
    }
    if (name_.substr(0,6)=="VMRD_B") {
      disk_=-(name_[7]-'0');
    }
    //cout << "name "<<name<<" "<<layer_<<" "<<disk_<<endl;
    assert(disk_!=0||layer_!=0);
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="allstubout"||
	output=="allstuboutn1"||
	output=="allstuboutn2"||
	output=="allstuboutn3"||
	output=="allstuboutn4"||
	output=="allstuboutn5"||
	output=="allstuboutn6"
	){
      FPGAAllStubs* tmp=dynamic_cast<FPGAAllStubs*>(memory);
      assert(tmp!=0);
      allstubs_.push_back(tmp);
      return;
    }
    if (output=="vmstuboutPHI1X1n1"||
	output=="vmstuboutPHI1X1n2"||
	output=="vmstuboutPHI1X1n3"||
	output=="vmstuboutPHI1X1n4"||
	output=="vmstuboutPHI1X1n5"||
	output=="vmstuboutPHI1X1n6"||
	output=="vmstuboutPHI1X1n7"||
	output=="vmstuboutPHI1X1n8"||
	output=="vmstuboutPHI1X1n9"||
	output=="vmstuboutPHI1X1n10"||
	output=="vmstuboutPHI1X1n11"||
	output=="vmstuboutPHI1X1n12") {
      FPGAVMStubs* tmp=dynamic_cast<FPGAVMStubs*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmstubsPHI1Z1_.push_back(tmp);
      }
      if (disk_!=0) {
	vmstubsPHI1R1_.push_back(tmp);
      }
      return;
    }
    if (output=="vmstuboutPHI1X2n1"||
	output=="vmstuboutPHI1X2n2"||
	output=="vmstuboutPHI1X2n3"||
	output=="vmstuboutPHI1X2n4"||
	output=="vmstuboutPHI1X2n5"||
	output=="vmstuboutPHI1X2n6"||
	output=="vmstuboutPHI1X2n7"||
	output=="vmstuboutPHI1X2n8"||
	output=="vmstuboutPHI1X2n9"||
	output=="vmstuboutPHI1X2n10"||
	output=="vmstuboutPHI1X2n11"||
	output=="vmstuboutPHI1X2n12") {
      FPGAVMStubs* tmp=dynamic_cast<FPGAVMStubs*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmstubsPHI1Z2_.push_back(tmp);
      }
      if (disk_!=0) {
	vmstubsPHI1R2_.push_back(tmp);
      }
      return;
    }
    if (output=="vmstuboutPHI2X1n1"||
	output=="vmstuboutPHI2X1n2"||
	output=="vmstuboutPHI2X1n3"||
	output=="vmstuboutPHI2X1n4"||
	output=="vmstuboutPHI2X1n5"||
	output=="vmstuboutPHI2X1n6"||
	output=="vmstuboutPHI2X1n7"||
	output=="vmstuboutPHI2X1n8"||
	output=="vmstuboutPHI2X1n9"||
	output=="vmstuboutPHI2X1n10"||
	output=="vmstuboutPHI2X1n11"||
	output=="vmstuboutPHI2X1n12") {
      FPGAVMStubs* tmp=dynamic_cast<FPGAVMStubs*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmstubsPHI2Z1_.push_back(tmp);
      }
      if (disk_!=0) {
	vmstubsPHI2R1_.push_back(tmp);
      }
      return;
    }
    if (output=="vmstuboutPHI2X2n1"||
	output=="vmstuboutPHI2X2n2"||
	output=="vmstuboutPHI2X2n3"||
	output=="vmstuboutPHI2X2n4"||
	output=="vmstuboutPHI2X2n5"||
	output=="vmstuboutPHI2X2n6"||
	output=="vmstuboutPHI2X2n7"||
	output=="vmstuboutPHI2X2n8"||
	output=="vmstuboutPHI2X2n9"||
	output=="vmstuboutPHI2X2n10"||
	output=="vmstuboutPHI2X2n11"||
	output=="vmstuboutPHI2X2n12") {
      FPGAVMStubs* tmp=dynamic_cast<FPGAVMStubs*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmstubsPHI2Z2_.push_back(tmp);
      }
      if (disk_!=0) {
	vmstubsPHI2R2_.push_back(tmp);
      }
      return;
    }
    if (output=="vmstuboutPHI3X1n1"||
	output=="vmstuboutPHI3X1n2"||
	output=="vmstuboutPHI3X1n3"||
	output=="vmstuboutPHI3X1n4"||
	output=="vmstuboutPHI3X1n5"||
	output=="vmstuboutPHI3X1n6"||
	output=="vmstuboutPHI3X1n7"||
	output=="vmstuboutPHI3X1n8"||
	output=="vmstuboutPHI3X1n9"||
	output=="vmstuboutPHI3X1n10"||
	output=="vmstuboutPHI3X1n11"||
	output=="vmstuboutPHI3X1n12") {
      FPGAVMStubs* tmp=dynamic_cast<FPGAVMStubs*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmstubsPHI3Z1_.push_back(tmp);
      }
      if (disk_!=0) {
	vmstubsPHI3R1_.push_back(tmp);
      }
      return;
    }
    if (output=="vmstuboutPHI3X2n1"||
	output=="vmstuboutPHI3X2n2"||
	output=="vmstuboutPHI3X2n3"||
	output=="vmstuboutPHI3X2n4"||
	output=="vmstuboutPHI3X2n5"||
	output=="vmstuboutPHI3X2n6"||
	output=="vmstuboutPHI3X2n7"||
	output=="vmstuboutPHI3X2n8"||
	output=="vmstuboutPHI3X2n9"||
	output=="vmstuboutPHI3X2n10"||
	output=="vmstuboutPHI3X2n11"||
	output=="vmstuboutPHI3X2n12") {
      FPGAVMStubs* tmp=dynamic_cast<FPGAVMStubs*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmstubsPHI3Z2_.push_back(tmp);
      }
      if (disk_!=0) {
	vmstubsPHI3R2_.push_back(tmp);
      }
      return;
    }
    if (output=="vmstuboutPHI4X1n1"||
	output=="vmstuboutPHI4X1n2"||
	output=="vmstuboutPHI4X1n3"||
	output=="vmstuboutPHI4X1n4"||
	output=="vmstuboutPHI4X1n5"||
	output=="vmstuboutPHI4X1n6"||
	output=="vmstuboutPHI4X1n7"||
	output=="vmstuboutPHI4X1n8"||
	output=="vmstuboutPHI4X1n9"||
	output=="vmstuboutPHI4X1n10"||
	output=="vmstuboutPHI4X1n11"||
	output=="vmstuboutPHI4X1n12") {
      FPGAVMStubs* tmp=dynamic_cast<FPGAVMStubs*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmstubsPHI4Z1_.push_back(tmp);
      }
      if (disk_!=0) {
	vmstubsPHI4R1_.push_back(tmp);
      }
      return;
    }
    if (output=="vmstuboutPHI4X2n1"||
	output=="vmstuboutPHI4X2n2"||
	output=="vmstuboutPHI4X2n3"||
	output=="vmstuboutPHI4X2n4"||
	output=="vmstuboutPHI4X2n5"||
	output=="vmstuboutPHI4X2n6"||
	output=="vmstuboutPHI4X2n7"||
	output=="vmstuboutPHI4X2n8"||
	output=="vmstuboutPHI4X2n9"||
	output=="vmstuboutPHI4X2n10"||
	output=="vmstuboutPHI4X2n11"||
	output=="vmstuboutPHI4X2n12") {
      FPGAVMStubs* tmp=dynamic_cast<FPGAVMStubs*>(memory);
      assert(tmp!=0);
      if (layer_!=0) {
	vmstubsPHI4Z2_.push_back(tmp);
      }
      if (disk_!=0) {
	vmstubsPHI4R2_.push_back(tmp);
      }
      return;
    }
    cout << "Could not find : "<<output<<endl;
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="stubinLink1"||input=="stubinLink2"||input=="stubinLink3"){
      FPGAStubLayer* tmp1=dynamic_cast<FPGAStubLayer*>(memory);
      assert(tmp1!=0);
      if (name_[3]=='_'){
	stubinputs_.push_back(tmp1);
      } else if (name_[3]=='D'){
	stubinputsdisk_.push_back(tmp1);
      } else {
	assert(0);
      }
      return;
    }
    cout << "Could not find input : "<<input<<endl;
    assert(0);
  }

  void execute(){

    //only one of these can be filled!
    assert(stubinputsdisk_.size()*stubinputs_.size()==0);

    unsigned int count=0;

    //cout << "FPGAVMRouter : "<<getName()
    //	 <<" stubinputs_.size() = "<<stubinputs_.size()
    //	 <<" stubinputsdisk_.size() = "<<stubinputsdisk_.size()
    //	 <<endl;

    if (stubinputs_.size()!=0){
      for(unsigned int j=0;j<stubinputs_.size();j++){
	for(unsigned int i=0;i<stubinputs_[j]->nStubs();i++){
	  count++;
	  if (count>MAXVMROUTER) continue;
	  std::pair<FPGAStub*,L1TStub*> stub=stubinputs_[j]->getStub(i);
	  //FIXME Next few lines should be member data in FPGAStub...
	  int iz=4+(stub.first->z().value()>>(stub.first->z().nbits()-3));
	  int iphitmp=stub.first->phi().value();
	  int layer=stub.first->layer().value()+1;
	  if ((layer%2)==1) iphitmp-=(1<<(stub.first->phi().nbits()-3));  
	  assert(iphitmp>=0);
	  int iphi=iphitmp>>(stub.first->phi().nbits()-2);

	  //cout << "iphi iz : "<<iphi<<" "<<iz<<endl;
	  assert(iz>=0);
	  assert(iz<=7);
	  iz=iz%2;
	  assert(iphi>=0);
	  assert(iphi<=3);

	  stub.first->setAllStubIndex(allstubs_[0]->nStubs());
	  stub.second->setAllStubIndex(allstubs_[0]->nStubs());

	  for (unsigned int l=0;l<allstubs_.size();l++){
	    allstubs_[l]->addStub(stub);
	  }

	  bool insert=false;

	  //cout << "FPGAVMRouter "<<getName()<<" "<<iphi<<" "<<iz<<endl;

	  if (iphi==0&&iz==0) {
	    for (unsigned int l=0;l<vmstubsPHI1Z1_.size();l++){
	      vmstubsPHI1Z1_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI1Z1_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  if (iphi==0&&iz==1) {
	    for (unsigned int l=0;l<vmstubsPHI1Z2_.size();l++){
	      vmstubsPHI1Z2_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI1Z2_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  
	  if (iphi==1&&iz==0) {
	    for (unsigned int l=0;l<vmstubsPHI2Z1_.size();l++){
	      vmstubsPHI2Z1_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI2Z1_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  if (iphi==1&&iz==1) {
	    for (unsigned int l=0;l<vmstubsPHI2Z2_.size();l++){
	      vmstubsPHI2Z2_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI2Z2_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  
	  if (iphi==2&&iz==0) {
	    for (unsigned int l=0;l<vmstubsPHI3Z1_.size();l++){
	      vmstubsPHI3Z1_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI3Z1_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  if (iphi==2&&iz==1) {
	    for (unsigned int l=0;l<vmstubsPHI3Z2_.size();l++){
	      vmstubsPHI3Z2_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI3Z2_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  
	  if (iphi==3&&iz==0) {
	    for (unsigned int l=0;l<vmstubsPHI4Z1_.size();l++){
	      vmstubsPHI4Z1_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI4Z1_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  if (iphi==3&&iz==1) {
	    for (unsigned int l=0;l<vmstubsPHI4Z2_.size();l++){
	      vmstubsPHI4Z2_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI4Z2_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  assert(insert);
	}
      }

    }
    if (stubinputsdisk_.size()>0) {
      //cout << "Routing stubs in disk" <<endl;
      for(unsigned int j=0;j<stubinputsdisk_.size();j++){
	for(unsigned int i=0;i<stubinputsdisk_[j]->nStubs();i++){
	  //cout << "Found stub in disk in "<<getName()<<endl;
	  std::pair<FPGAStub*,L1TStub*> stub=stubinputsdisk_[j]->getStub(i);
	  //
          bool isSS = name_[8] == '6' || name_[8] == '8';
          if(isSS && stub.second->r() <60) {
            cout<<"in SS router "<<name_<<" but r = "<<stub.second->r()<<"\n";
            assert(0);
          }
          if(!isSS && stub.second->r() >60) {
            cout<<"in PS router "<<name_<<" but r = "<<stub.second->r()<<"\n";
            assert(0);
          }
	  //FIXME Next few lines should be member data in FPGAStub...
	  int irtmp=stub.first->r().value(); 
	  int ir;
	  if(stub.second->r()>60){//SS module, r encoded differently
	    ir=2+(irtmp>>3);
	    //cout<<"That's SS module! "<<stub.second->r()<<": "<<irtmp<<" with "<<stub.first->r().nbits()<<"\n";
	  }
	  else {
	    ir=irtmp>>(stub.first->r().nbits()-2);
	  }
	  int iphitmp=stub.first->phi().value();
	  int disk=stub.first->disk().value();
	  if ((disk%2)==0) iphitmp-=(1<<(stub.first->phi().nbits()-3));  
	  assert(iphitmp>=0);
	  int iphi=iphitmp>>(stub.first->phi().nbits()-2);

	  //cout << "iphi ir irtmp: "<<iphi<<" "<<ir<<"\t"<<irtmp<<"\t ("<<stub.second->r()<<")"<<endl;
	  assert(ir>=0);
	  assert(ir<=3);
	  ir=ir%2;
	  assert(iphi>=0);
	  assert(iphi<=3);

	  stub.first->setAllStubIndex(allstubs_[0]->nStubs());
	  stub.second->setAllStubIndex(allstubs_[0]->nStubs());
	  
	  for (unsigned int l=0;l<allstubs_.size();l++){
	    allstubs_[l]->addStub(stub);
	  }

	  bool insert=false;
	  
	  if (iphi==0&&ir==0) {
	    for (unsigned int l=0;l<vmstubsPHI1R1_.size();l++){
	      vmstubsPHI1R1_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI1R1_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  if (iphi==0&&ir==1) {
	    for (unsigned int l=0;l<vmstubsPHI1R2_.size();l++){
	      vmstubsPHI1R2_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI1R2_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  
	  if (iphi==1&&ir==0) {
	    for (unsigned int l=0;l<vmstubsPHI2R1_.size();l++){
	      vmstubsPHI2R1_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI2R1_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  if (iphi==1&&ir==1) {
	    for (unsigned int l=0;l<vmstubsPHI2R2_.size();l++){
	      vmstubsPHI2R2_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI2R2_[l]->getName() << endl;	      
	      insert=true;
	    }
	  }
	  
	  if (iphi==2&&ir==0) {
	    for (unsigned int l=0;l<vmstubsPHI3R1_.size();l++){
	      vmstubsPHI3R1_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI3R1_[l]->getName() << endl;
	      insert=true;

	    }
	  }
	  if (iphi==2&&ir==1) {
	    for (unsigned int l=0;l<vmstubsPHI3R2_.size();l++){
	      vmstubsPHI3R2_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI3R2_[l]->getName() << endl;
	      insert=true;
	    }
	  }
	  
	  if (iphi==3&&ir==0) {
	    for (unsigned int l=0;l<vmstubsPHI4R1_.size();l++){
	      vmstubsPHI4R1_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI4R1_[l]->getName() << endl;	   
	      insert=true; 
	    }
	  }
	  if (iphi==3&&ir==1) {
	    for (unsigned int l=0;l<vmstubsPHI4R2_.size();l++){
	      vmstubsPHI4R2_[l]->addStub(stub);
	      //cout << "Adding stub in " << vmstubsPHI4R1_[l]->getName() << endl;	 
	      insert=true;
	    }
	  }
	  assert(insert);
	}
      }

    }

      if (writeVMOccupancy) {
	static ofstream out("vmoccupancy.txt");
	if (stubinputs_.size()!=0){
	  out<<vmstubsPHI1Z1_[0]->getName()<<" "<<vmstubsPHI1Z1_[0]->nStubs()<<endl;
	  out<<vmstubsPHI1Z2_[0]->getName()<<" "<<vmstubsPHI1Z2_[0]->nStubs()<<endl;
	  out<<vmstubsPHI2Z1_[0]->getName()<<" "<<vmstubsPHI2Z1_[0]->nStubs()<<endl;
	  out<<vmstubsPHI2Z2_[0]->getName()<<" "<<vmstubsPHI2Z2_[0]->nStubs()<<endl;
	  out<<vmstubsPHI3Z1_[0]->getName()<<" "<<vmstubsPHI3Z1_[0]->nStubs()<<endl;
	  out<<vmstubsPHI3Z2_[0]->getName()<<" "<<vmstubsPHI3Z2_[0]->nStubs()<<endl;
	  if ((vmstubsPHI3Z2_[0]->getName()[5]-'0')%2==0) {
	    out<<vmstubsPHI4Z1_[0]->getName()<<" "<<vmstubsPHI4Z1_[0]->nStubs()<<endl;
	    out<<vmstubsPHI4Z2_[0]->getName()<<" "<<vmstubsPHI4Z2_[0]->nStubs()<<endl;
	  }
	}
	if (stubinputsdisk_.size()>0) {
	  out<<vmstubsPHI1R1_[0]->getName()<<" "<<vmstubsPHI1R1_[0]->nStubs()<<endl;
	  out<<vmstubsPHI1R2_[0]->getName()<<" "<<vmstubsPHI1R2_[0]->nStubs()<<endl;
	  out<<vmstubsPHI2R1_[0]->getName()<<" "<<vmstubsPHI2R1_[0]->nStubs()<<endl;
	  out<<vmstubsPHI2R2_[0]->getName()<<" "<<vmstubsPHI2R2_[0]->nStubs()<<endl;
	  out<<vmstubsPHI3R1_[0]->getName()<<" "<<vmstubsPHI3R1_[0]->nStubs()<<endl;
	  out<<vmstubsPHI3R2_[0]->getName()<<" "<<vmstubsPHI3R2_[0]->nStubs()<<endl;
	  if ((vmstubsPHI3R2_[0]->getName()[5]-'0')%2==1) {
	    out<<vmstubsPHI4R1_[0]->getName()<<" "<<vmstubsPHI4R1_[0]->nStubs()<<endl;
	    out<<vmstubsPHI4R2_[0]->getName()<<" "<<vmstubsPHI4R2_[0]->nStubs()<<endl;
	  }
	}
      }


    if (writeAllStubs) {
      static ofstream out("allstubs.txt");
      out<<allstubs_[0]->getName()<<" "<<allstubs_[0]->nStubs()<<endl;
    }
 


  }



private:

  int layer_;
  int disk_;

  vector<FPGAStubLayer*> stubinputs_;
  vector<FPGAStubLayer*> stubinputsdisk_;
  vector<FPGAAllStubs*> allstubs_;

  vector<FPGAVMStubs*> vmstubsPHI1Z1_;
  vector<FPGAVMStubs*> vmstubsPHI1Z2_;
  vector<FPGAVMStubs*> vmstubsPHI2Z1_;
  vector<FPGAVMStubs*> vmstubsPHI2Z2_;
  vector<FPGAVMStubs*> vmstubsPHI3Z1_;
  vector<FPGAVMStubs*> vmstubsPHI3Z2_;
  vector<FPGAVMStubs*> vmstubsPHI4Z1_;
  vector<FPGAVMStubs*> vmstubsPHI4Z2_;

  vector<FPGAVMStubs*> vmstubsPHI1R1_;
  vector<FPGAVMStubs*> vmstubsPHI1R2_;
  vector<FPGAVMStubs*> vmstubsPHI2R1_;
  vector<FPGAVMStubs*> vmstubsPHI2R2_;
  vector<FPGAVMStubs*> vmstubsPHI3R1_;
  vector<FPGAVMStubs*> vmstubsPHI3R2_;
  vector<FPGAVMStubs*> vmstubsPHI4R1_;
  vector<FPGAVMStubs*> vmstubsPHI4R2_;


};

#endif


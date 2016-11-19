#ifndef FPGATRACKDERTABLE_H
#define FPGATRACKDERTABLE_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>
#include "FPGATrackDer.hh"

using namespace std;

class FPGATrackDerTable{

public:

  FPGATrackDerTable() {

    Nlay_=6;
    Ndisk_=5;

    LayerMemBits_=5;
    DiskMemBits_=7;
    
    LayerDiskMemBits_=15;

    alphaBits_=alphaBitsTable;

    nextLayerValue_=0;
    nextDiskValue_=0;
    nextLayerDiskValue_=0;
    lastMultiplicity_=(1<<(3*alphaBits_));


    for(int i=0;i<(1<<Nlay_);i++){
      LayerMem_.push_back(-1);
    }

    for(int i=0;i<(1<<(2*Ndisk_));i++){
      DiskMem_.push_back(-1);
    }

    for(int i=0;i<(1<<(LayerMemBits_+DiskMemBits_));i++){
      LayerDiskMem_.push_back(-1);
    }
   
  }

  ~FPGATrackDerTable() {

  }

  FPGATrackDer* getDerivatives(int index){
    return &derivatives_[index];
  }

  FPGATrackDer* getDerivatives(unsigned int layermask, 
			       unsigned int diskmask,
			       unsigned int alphaindex){
    int index=getIndex(layermask,diskmask);
    //if (index<0||index!=17984||alphaindex!=20) {
    if (index<0) {
      return 0;
    }
    //cout << "getDerivatives index alphaindex "<<index<<" "<<alphaindex<<endl;
    return &derivatives_[index+alphaindex];
  }


  int getIndex(unsigned int layermask,unsigned int diskmask) {

    assert(layermask<LayerMem_.size());

    assert(diskmask<DiskMem_.size());

    int layercode=LayerMem_[layermask];
    int diskcode=DiskMem_[diskmask];

    if (diskcode<0||layercode<0) {
      cout << "layermask diskmask : "<<layermask<<" "<<diskmask<<endl;
      return -1;
    }

    assert(layercode>=0);
    assert(layercode<(1<<LayerMemBits_));
    assert(diskcode>=0);
    assert(diskcode<(1<<DiskMemBits_));

    int layerdiskaddress=layercode+(diskcode<<LayerMemBits_);

    assert(layerdiskaddress>=0);
    assert(layerdiskaddress<(1<<(LayerMemBits_+DiskMemBits_)));

    int address=LayerDiskMem_[layerdiskaddress];

    if (address<0) {
      cout << "layermask diskmask : "<<layermask<<" "<<diskmask<<endl;
      return -1;
    }

    assert(address>=0);
    //cout << "address LayerDiskMemBits_ : "<<address<<" "<<LayerDiskMemBits_<<endl;
    assert(address<(1<<LayerDiskMemBits_));

    return address;

  }

  void addEntry(unsigned int layermask, unsigned int diskmask, int multiplicity){

    assert(multiplicity<=(1<<(3*alphaBits_)));

    assert(layermask<(unsigned int)(1<<Nlay_));

    assert(diskmask<(unsigned int)(1<<(2*Ndisk_)));

    if (LayerMem_[layermask]==-1) {
      LayerMem_[layermask]=nextLayerValue_++;
    }
    if (DiskMem_[diskmask]==-1) {
      DiskMem_[diskmask]=nextDiskValue_++;
    }

    int layercode=LayerMem_[layermask];
    int diskcode=DiskMem_[diskmask];

    assert(layercode>=0);
    assert(layercode<(1<<LayerMemBits_));
    assert(diskcode>=0);
    assert(diskcode<(1<<DiskMemBits_));

    int layerdiskaddress=layercode+(diskcode<<LayerMemBits_);

    assert(layerdiskaddress>=0);
    assert(layerdiskaddress<(1<<(LayerMemBits_+DiskMemBits_)));

    int address=LayerDiskMem_[layerdiskaddress];

    if (address!=-1) {
      cout << "Duplicate entry:  layermask="
	   <<layermask<<" diskmaks="<<diskmask<<endl;
    }

    assert(address==-1);  //Should not already have this one!

    LayerDiskMem_[layerdiskaddress]=nextLayerDiskValue_;

    nextLayerDiskValue_+=multiplicity;

    lastMultiplicity_=multiplicity;

    for(int i=0;i<multiplicity;i++) {
      FPGATrackDer tmp;
      tmp.setIndex(layermask,diskmask,i);
      derivatives_.push_back(tmp);
    }

  }

  void readPatternFile(std::string fileName){

    ifstream in(fileName.c_str());
    cout<<"reading fit pattern file "<<fileName<<"\n";
    cout<<"  flags (good/eof/fail/bad): "<<in.good()<<" "<<in.eof()<<" "<<in.fail()<<" "<<in.bad()<<"\n"; 

    while (in.good()) {

      std::string layerstr,diskstr;
      int multiplicity;

      in >>layerstr>>diskstr>>multiplicity;

      if (alphaBits_==2) {
	if (multiplicity==8) multiplicity=4;
	if (multiplicity==64) multiplicity=16;
	if (multiplicity==512) multiplicity=64;
      }

      if (alphaBits_==1) {
	if (multiplicity==8) multiplicity=2;
	if (multiplicity==64) multiplicity=4;
	if (multiplicity==512) multiplicity=8;
      }


      if (!in.good()) continue;
      
      char** tmpptr=0;

      int layers=strtol(layerstr.c_str(), tmpptr, 2); 
      int disks=strtol(diskstr.c_str(), tmpptr, 2); 

      //cout << "adding: "<<layers<<" "<<disks<<" "<<multiplicity<<endl;
   
      addEntry(layers,disks,multiplicity);

    }

  }

  int getEntries() const {
    return nextLayerDiskValue_;
  }

  void fillTable() {

    int nentries=getEntries();

    for (int i=0;i<nentries;i++){
      FPGATrackDer& der=derivatives_[i];
      int layermask=der.getLayerMask();
      int diskmask=der.getDiskMask();
      int alphamask=der.getAlphaMask();

      bool print=getIndex(layermask,diskmask)==17984&&alphamask==20;
      print=false;

      if (print) {
	cout << "i "<<i<<" "<<layermask<<" "<<diskmask<<" "
	     <<alphamask<<" "<<print<<endl;
      }

      int nlayers=0;
      //int layers[6];
      double r[6];

      for (unsigned l=0;l<6;l++){
	if (layermask&(1<<(5-l))) {
	  //layers[nlayers]=l+1;
	  r[nlayers]=rmean[l];
	  //cout << "Hit in layer "<<layers[nlayers]<<" "<<r[nlayers]<<endl;
	  nlayers++;  
	}
      }

      int ndisks=0;
      //int disks[5];
      double z[5];
      double alpha[5];

      //double t=sinh(1.2);
      //double rinv=-0.0057/2.5;

      double t=gett(diskmask);
      double rinv=0.00000001;
      //cout << "layermask diskmask t :"<<layermask<<" "<<diskmask<<" "<<t<<endl;
      
      for (unsigned d=0;d<5;d++){
	if (diskmask&(3<<(2*(4-d)))) {
	  //disks[ndisks]=d+1;
	  z[ndisks]=zmean[d];
	  alpha[ndisks]=0.0;
	  if (diskmask&(1<<(2*(4-d)))) {
	    if (alphaBits_==3) {
	      int ialpha=alphamask&7;
	      alphamask=alphamask>>3;
	      double r=zmean[d]/t;
	      alpha[ndisks]=480*0.009*(ialpha-3.5)/(4.0*r*r);
	    }
	    if (alphaBits_==2) {
	      int ialpha=alphamask&3;
	      alphamask=alphamask>>2;
	      double r=zmean[d]/t;
	      alpha[ndisks]=480*0.009*(ialpha-1.5)/(4.0*r*r);
	    }
	    if (alphaBits_==1) {
	      int ialpha=alphamask&1;
	      alphamask=alphamask>>1;
	      double r=zmean[d]/t;
	      alpha[ndisks]=480*0.009*(ialpha-0.5)/(4.0*r*r);
	    }
	    //if (d==0) alpha[ndisks]=-0.000220;
	    //if (d==1) alpha[ndisks]=0.000244;
	    //if (print) {
	    //  cout << "Hit in disk "<<z[ndisks]<<" "
	    //	   <<ialpha<<" "<<alpha[ndisks]<<endl;
	    //}
	  }
	  ndisks++;  
	}
      }


      double D[4][12];
      double MinvDt[4][12];
      int iMinvDt[4][12];
      

      if (print) {
	for(int ii=0;ii<nlayers;ii++){
	  cout << "Layer r : "<<r[ii]<<endl;
	}
	for(int ii=0;ii<ndisks;ii++){
	  cout << "Disk z alpha : "<<z[ii]<<" "<<alpha[ii]<<endl;
	}
      }

      calculateDerivativesTable(nlayers,r,ndisks,z,alpha,t,rinv,D,MinvDt,iMinvDt);


      for(int j=0;j<nlayers+ndisks;j++){
	if (print) {
	  cout << "Table "<<endl;
	  cout << MinvDt[0][2*j] <<" "
	       << MinvDt[1][2*j] <<" "
	       << MinvDt[2][2*j] <<" "
	       << MinvDt[3][2*j] <<" "
	       <<endl;
	  cout << MinvDt[0][2*j+1] <<" "
	       << MinvDt[1][2*j+1] <<" "
	       << MinvDt[2][2*j+1] <<" "
	       << MinvDt[3][2*j+1] <<" "
	       <<endl;
	}
	
	der.sett(t);

	//integer
	assert(abs(iMinvDt[0][2*j])<(1<<15));
	assert(abs(iMinvDt[0][2*j+1])<(1<<15));
	assert(abs(iMinvDt[1][2*j])<(1<<15));
	assert(abs(iMinvDt[1][2*j+1])<(1<<15));
	assert(abs(iMinvDt[2][2*j])<(1<<15));
	assert(abs(iMinvDt[2][2*j+1])<(1<<15));
	assert(abs(iMinvDt[3][2*j])<(1<<15));
	assert(abs(iMinvDt[3][2*j+1])<(1<<15));

	der.setirinvdphi(j,iMinvDt[0][2*j]); 
	der.setirinvdzordr(j,iMinvDt[0][2*j+1]); 
	der.setiphi0dphi(j,iMinvDt[1][2*j]); 
	der.setiphi0dzordr(j,iMinvDt[1][2*j+1]); 
	der.setitdphi(j,iMinvDt[2][2*j]); 
	der.setitdzordr(j,iMinvDt[2][2*j+1]); 
	der.setiz0dphi(j,iMinvDt[3][2*j]); 
	der.setiz0dzordr(j,iMinvDt[3][2*j+1]); 
	//floating point
	der.setrinvdphi(j,MinvDt[0][2*j]); 
	der.setrinvdzordr(j,MinvDt[0][2*j+1]); 
	der.setphi0dphi(j,MinvDt[1][2*j]); 
	der.setphi0dzordr(j,MinvDt[1][2*j+1]); 
	der.settdphi(j,MinvDt[2][2*j]); 
	der.settdzordr(j,MinvDt[2][2*j+1]); 
	der.setz0dphi(j,MinvDt[3][2*j]); 
	der.setz0dzordr(j,MinvDt[3][2*j+1]); 
      }



    }

	/*    
    if (writeFitDerTable) {
      for (unsigned int seedlayer=1;seedlayer<=5;seedlayer+=2) {
	for(unsigned int j=0;j<8;j++) {
	  std::ostringstream oss;
	  oss << "FitDerTable_L"<<seedlayer<<"_"<<j+1<<".txt";
	  std::string fname=oss.str();
	  ofstream out(fname.c_str());
	  for(unsigned int imask=0;imask<16;imask++) {
	    int nmatches=0;
	    for (unsigned int i=0;i<4;i++) {
	      if (((1<<i)&imask)!=0) nmatches++;
	    }
	    unsigned int ifullmask=0;
	    if (seedlayer==5) ifullmask=(imask<<2)+3;
	    if (seedlayer==3) ifullmask=(imask&12)*4+12+(imask&3);
	    if (seedlayer==1) ifullmask=imask+48;
	    
	    unsigned int layer=j/2;
	    if (seedlayer==1) layer+=2;
	    if (seedlayer==3&&layer>1) layer+=2;
	    assert(layer<6);

	    bool hitlayer=((1<<(5-layer))&ifullmask)!=0;

	    FPGATrackDer* der=0;
	    if (hitlayer&&(nmatches>1)) {
	      der=getDerivatives(ifullmask,0,0); 
	    }
	    

	   
	    if (der!=0) {

	      unsigned int index=0;
	      for (unsigned int i=0;i<layer;i++) {
		if (((1<<(5-i))&ifullmask)!=0) index++;
		//cout << "layer index : "<<layer<<" "<<index<<endl;
	      }
	      
	      //cout << "seedlayer j imask ifullmask layer index "<<seedlayer<<" "
	      //		   <<j<<" "<<imask<<" "<<ifullmask<<" "<<layer<<" "
	      //   <<index<<endl;


	      //for (unsigned int i=0;i<6;i++) {
	      //	cout <<der->getirinvdphi(i)<<" ";
	      //}
	      //cout << endl;

	      assert(index<6);



	      FPGAWord tmp1,tmp2,tmp3,tmp4;
	      if (j%2==0) {
		tmp1.set(der->getirinvdphi(index),15,false); 
		tmp2.set(der->getiphi0dphi(index),15,false); 
		tmp3.set(der->getitdphi(index),15,false); 
		tmp4.set(der->getiz0dphi(index),15,false); 
	      } else {
		tmp1.set(der->getirinvdzordr(index),15,false); 
		tmp2.set(der->getiphi0dzordr(index),15,false); 
		tmp3.set(der->getitdzordr(index),15,false); 
		tmp4.set(der->getiz0dzordr(index),15,false); 
	      }
	      //FPGAWord tmp;
	      //tmp.set(imask,4,true);
	      out //<< tmp.str()<<" "<<index<<" "
		  <<tmp1.str()<<tmp2.str()<<tmp3.str()<<tmp4.str()<<endl;
	   
	    } else {
	      //FPGAWord tmp;
	      //tmp.set(ifullmask,6,true);
	      out //<<tmp.str()<<" "<<hitlayer
		<<"000000000000000000000000000000000000000000000000000000000000"<<endl;
	    }
	  }
	}
      }
    }
	*/


    //New table format with disks
	/*
	if (writeFitDerTable) {

      ofstream outL("FitDerTableNew_LayerMem.txt");
      for (unsigned int i=0;i<LayerMem_.size();i++){
	FPGAWord tmp;
	int tmp1=LayerMem_[i];
	if (tmp1<0) tmp1=(1<<6)-1;
	cout << "i LayerMem_ : "<<i<<" "<<tmp1<<endl;
	tmp.set(tmp1,6,true,__LINE__,__FILE__);
	outL << tmp.str() <<endl;
      }
      outL.close();


      ofstream outD("FitDerTableNew_DiskMem.txt");
      for (unsigned int i=0;i<DiskMem_.size();i++){
	int tmp1=DiskMem_[i];
	if (tmp1<0) tmp1=(1<<10)-1;
	FPGAWord tmp;
	tmp.set(tmp1,10,true,__LINE__,__FILE__);
	outD << tmp.str() <<endl;
      }
      outD.close();


      ofstream outLD("FitDerTableNew_LayerDiskMem.txt");
      for (unsigned int i=0;i<LayerDiskMem_.size();i++){
	int tmp1=LayerDiskMem_[i];
	if (tmp1<0) tmp1=(1<<15)-1;
	FPGAWord tmp;
	tmp.set(tmp1,15,true,__LINE__,__FILE__);
	outLD << tmp.str() <<endl;
      }
      outLD.close();


      unsigned int nderivatives=derivatives_.size();

      cout << "nderivatives = "<<nderivatives<<endl;

      ofstream outrinvdphi[6];
      outrinvdphi[0].open("FitDerTableNew_Rinvdphi_1.txt");
      outrinvdphi[1].open("FitDerTableNew_Rinvdphi_2.txt");
      outrinvdphi[2].open("FitDerTableNew_Rinvdphi_3.txt");
      outrinvdphi[3].open("FitDerTableNew_Rinvdphi_4.txt");
      outrinvdphi[4].open("FitDerTableNew_Rinvdphi_5.txt");
      outrinvdphi[5].open("FitDerTableNew_Rinvdphi_6.txt");

      ofstream outrinvdzordr[6];
      outrinvdzordr[0].open("FitDerTableNew_Rinvdzordr_1.txt");
      outrinvdzordr[1].open("FitDerTableNew_Rinvdzordr_2.txt");
      outrinvdzordr[2].open("FitDerTableNew_Rinvdzordr_3.txt");
      outrinvdzordr[3].open("FitDerTableNew_Rinvdzordr_4.txt");
      outrinvdzordr[4].open("FitDerTableNew_Rinvdzordr_5.txt");
      outrinvdzordr[5].open("FitDerTableNew_Rinvdzordr_6.txt");

      ofstream outphi0dphi[6];
      outphi0dphi[0].open("FitDerTableNew_Phi0dphi_1.txt");
      outphi0dphi[1].open("FitDerTableNew_Phi0dphi_2.txt");
      outphi0dphi[2].open("FitDerTableNew_Phi0dphi_3.txt");
      outphi0dphi[3].open("FitDerTableNew_Phi0dphi_4.txt");
      outphi0dphi[4].open("FitDerTableNew_Phi0dphi_5.txt");
      outphi0dphi[5].open("FitDerTableNew_Phi0dphi_6.txt");

      ofstream outphi0dzordr[6];
      outphi0dzordr[0].open("FitDerTableNew_Phi0dzordr_1.txt");
      outphi0dzordr[1].open("FitDerTableNew_Phi0dzordr_2.txt");
      outphi0dzordr[2].open("FitDerTableNew_Phi0dzordr_3.txt");
      outphi0dzordr[3].open("FitDerTableNew_Phi0dzordr_4.txt");
      outphi0dzordr[4].open("FitDerTableNew_Phi0dzordr_5.txt");
      outphi0dzordr[5].open("FitDerTableNew_Phi0dzordr_6.txt");

      ofstream outtdphi[6];
      outtdphi[0].open("FitDerTableNew_Tdphi_1.txt");
      outtdphi[1].open("FitDerTableNew_Tdphi_2.txt");
      outtdphi[2].open("FitDerTableNew_Tdphi_3.txt");
      outtdphi[3].open("FitDerTableNew_Tdphi_4.txt");
      outtdphi[4].open("FitDerTableNew_Tdphi_5.txt");
      outtdphi[5].open("FitDerTableNew_Tdphi_6.txt");

      ofstream outtdzordr[6];
      outtdzordr[0].open("FitDerTableNew_Tdzordr_1.txt");
      outtdzordr[1].open("FitDerTableNew_Tdzordr_2.txt");
      outtdzordr[2].open("FitDerTableNew_Tdzordr_3.txt");
      outtdzordr[3].open("FitDerTableNew_Tdzordr_4.txt");
      outtdzordr[4].open("FitDerTableNew_Tdzordr_5.txt");
      outtdzordr[5].open("FitDerTableNew_Tdzordr_6.txt");

      ofstream outz0dphi[6];
      outz0dphi[0].open("FitDerTableNew_Z0dphi_1.txt");
      outz0dphi[1].open("FitDerTableNew_Z0dphi_2.txt");
      outz0dphi[2].open("FitDerTableNew_Z0dphi_3.txt");
      outz0dphi[3].open("FitDerTableNew_Z0dphi_4.txt");
      outz0dphi[4].open("FitDerTableNew_Z0dphi_5.txt");
      outz0dphi[5].open("FitDerTableNew_Z0dphi_6.txt");

      ofstream outz0dzordr[6];
      outz0dzordr[0].open("FitDerTableNew_Z0dzordr_1.txt");
      outz0dzordr[1].open("FitDerTableNew_Z0dzordr_2.txt");
      outz0dzordr[2].open("FitDerTableNew_Z0dzordr_3.txt");
      outz0dzordr[3].open("FitDerTableNew_Z0dzordr_4.txt");
      outz0dzordr[4].open("FitDerTableNew_Z0dzordr_5.txt");
      outz0dzordr[5].open("FitDerTableNew_Z0dzordr_6.txt");


      for (unsigned int i=0;i<nderivatives;i++){

	FPGATrackDer der=derivatives_[i];
	
	for (unsigned int j=0;j<6;j++){
	  FPGAWord tmprinvdphi;
	  int itmprinvdphi=der.getirinvdphi(j);
	  if (itmprinvdphi>(1<<13)) itmprinvdphi=(1<<13)-1;
	  tmprinvdphi.set(itmprinvdphi,14,false,__LINE__,__FILE__);
	  outrinvdphi[j] << tmprinvdphi.str() <<endl;
	  FPGAWord tmprinvdzordr;
	  int itmprinvdzordr=der.getirinvdzordr(j);
	  if (itmprinvdzordr>(1<<15)) itmprinvdzordr=(1<<15)-1;
	  tmprinvdzordr.set(itmprinvdzordr,16,false,__LINE__,__FILE__);
	  outrinvdzordr[j] << tmprinvdzordr.str() <<endl;

	  FPGAWord tmpphi0dphi;
	  int itmpphi0dphi=der.getiphi0dphi(j);
	  if (itmpphi0dphi>(1<<13)) itmpphi0dphi=(1<<13)-1;
	  tmpphi0dphi.set(itmpphi0dphi,14,false,__LINE__,__FILE__);
	  outphi0dphi[j] << tmpphi0dphi.str() <<endl;
	  FPGAWord tmpphi0dzordr;
	  int itmpphi0dzordr=der.getiphi0dzordr(j);
	  if (itmpphi0dzordr>(1<<15)) itmpphi0dzordr=(1<<15)-1;
	  tmpphi0dzordr.set(itmpphi0dzordr,16,false,__LINE__,__FILE__);
	  outphi0dzordr[j] << tmpphi0dzordr.str() <<endl;

	  FPGAWord tmptdphi;
	  int itmptdphi=der.getitdphi(j);
	  if (itmptdphi>(1<<13)) itmptdphi=(1<<13)-1;
	  tmptdphi.set(itmptdphi,14,false,__LINE__,__FILE__);
	  outtdphi[j] << tmptdphi.str() <<endl;
	  FPGAWord tmptdzordr;
	  int itmptdzordr=der.getitdzordr(j);
	  if (itmptdzordr>(1<<15)) itmptdzordr=(1<<15)-1;
	  tmptdzordr.set(itmptdzordr,16,false,__LINE__,__FILE__);
	  outtdzordr[j] << tmptdzordr.str() <<endl;

	  FPGAWord tmpz0dphi;
	  int itmpz0dphi=der.getiz0dphi(j);
	  if (itmpz0dphi>(1<<13)) itmpz0dphi=(1<<13)-1;
	  tmpz0dphi.set(itmpz0dphi,14,false,__LINE__,__FILE__);
	  outz0dphi[j] << tmpz0dphi.str() <<endl;
	  FPGAWord tmpz0dzordr;
	  int itmpz0dzordr=der.getiz0dzordr(j);
	  if (itmpz0dzordr>(1<<15)) itmpz0dzordr=(1<<15)-1;
	  tmpz0dzordr.set(itmpz0dzordr,16,false,__LINE__,__FILE__);
	  outz0dzordr[j] << tmpz0dzordr.str() <<endl;

	}

      }

      for (unsigned int j=0;j<6;j++){
	outrinvdphi[j].close();
	outrinvdzordr[j].close();

	outphi0dphi[j].close();
	outphi0dzordr[j].close();

	outtdphi[j].close();
	outtdzordr[j].close();

	outz0dphi[j].close();
	outz0dzordr[j].close();

      }

    }
	*/
		
	// -------------------------------------------------------------
	// Modified new table format with disks

	if (writeFitDerTable) {

		ofstream outL("FitDerTableNew_LayerMem.txt");
		for (unsigned int i=0;i<LayerMem_.size();i++){
			FPGAWord tmp;
			int tmp1=LayerMem_[i];
			if (tmp1<0) tmp1=(1<<5)-1;
			cout << "i LayerMem_ : "<<i<<" "<<tmp1<<endl;
			tmp.set(tmp1,5,true,__LINE__,__FILE__);
			outL << tmp.str() <<endl;
		}
		outL.close();

		ofstream outD("FitDerTableNew_DiskMem.txt");
		for (unsigned int i=0;i<DiskMem_.size();i++){
			int tmp1=DiskMem_[i];
			if (tmp1<0) tmp1=(1<<7)-1;
			FPGAWord tmp;
			tmp.set(tmp1,7,true,__LINE__,__FILE__);
			outD << tmp.str() <<endl;
		}
		outD.close();

		ofstream outLD("FitDerTableNew_LayerDiskMem.txt");
		for (unsigned int i=0;i<LayerDiskMem_.size();i++){
			int tmp1=LayerDiskMem_[i];
			if (tmp1<0) tmp1=(1<<15)-1;
			FPGAWord tmp;
			tmp.set(tmp1,15,true,__LINE__,__FILE__);
			outLD << tmp.str() <<endl;
		}
		outLD.close();

		unsigned int nderivatives = derivatives_.size();

		cout << "nderivatives = " << nderivatives << endl;

		// 6(+3) seeding pairs; each pair has max 4 inputs
		// if (some global parameter == "Forward")
		string seedings[] = {"L1L2", "L3L4", "L5L6", "F1F2", "F3F4", "L1F1"};
		// else
		// seedings = {"L1L2", "L3L4", "L5L6", "B1B2", "B3B4", "L1B1"};
		
		string prefix = "FitDerTableNewer_";
		
		ofstream outrinvdphi[6];
		for (unsigned int i = 0; i < 6; ++i) {
			string seed = seedings[i];
			string fname = prefix + "Rinvdphi_" + seed + ".txt";
			outrinvdphi[i].open(fname.c_str());
		}
		
		ofstream outrinvdzordr[6];
		for (unsigned int i = 0; i < 6; ++i) {
			string seed = seedings[i];
			string fname = prefix + "Rinvdzordr_" + seed + ".txt";
			outrinvdzordr[i].open(fname.c_str());
		}

		ofstream outphi0dphi[6];
		for (unsigned int i = 0; i < 6; ++i) {
			string seed = seedings[i];
			string fname = prefix + "Phi0dphi_"+seed+".txt";
			outphi0dphi[i].open(fname.c_str());
		}
		
		ofstream outphi0dzordr[6];
		for (unsigned int i = 0; i < 6; ++i) {
			string seed = seedings[i];
			string fname = prefix + "Phi0dzordr_" + seed + ".txt";
			outphi0dzordr[i].open(fname.c_str());
		}

		ofstream outtdphi[6];
		for (unsigned int i = 0; i < 6; ++i) {
			string seed = seedings[i];
			string fname = prefix + "Tdphi_" + seed + ".txt";
			outtdphi[i].open(fname.c_str());
		}

		ofstream outtdzordr[6];
		for (unsigned int i = 0; i < 6; ++i) {
			string seed = seedings[i];
			string fname = prefix + "Tdzordr_" + seed + ".txt";
			outtdzordr[i].open(fname.c_str());
		}
		
		ofstream outz0dphi[6];
		for (unsigned int i = 0; i < 6; ++i) {
			string seed = seedings[i];
			string fname = prefix + "Z0dphi_" + seed + ".txt";
			outz0dphi[i].open(fname.c_str());
		}

		ofstream outz0dzordr[6];
		for (unsigned int i = 0; i < 6; ++i) {
			string seed = seedings[i];
			string fname = prefix + "Z0dzordr_" + seed + ".txt";
			outz0dzordr[i].open(fname.c_str());
		}
		
		for (auto & der : derivatives_) {
				  
			unsigned int layerhits = der.getLayerMask();
			unsigned int diskmask = der.getDiskMask();
			unsigned int diskhits = 0;
			if ( diskmask&(3<<8) ) diskhits += 16;
			if ( diskmask&(3<<6) ) diskhits += 8;
			if ( diskmask&(3<<4) ) diskhits += 4;
			if ( diskmask&(3<<2) ) diskhits += 2;
			if ( diskmask&(3<<0) ) diskhits += 1;
			assert(diskhits < 32);
			unsigned int hits = (layerhits<<5) + diskhits;
			assert(hits < 4096);

			// loop over all seedings
			for (unsigned int i = 0; i < 6; ++i) {

				string seed = seedings[i];
				
				int itmprinvdphi[4] = {9999999,9999999,9999999,9999999};
				int itmprinvdzordr[4] = {9999999,9999999,9999999,9999999};
				int itmpphi0dphi[4] = {9999999,9999999,9999999,9999999};
				int itmpphi0dzordr[4] = {9999999,9999999,9999999,9999999};
				int itmptdphi[4] = {9999999,9999999,9999999,9999999};
				int itmptdzordr[4] = {9999999,9999999,9999999,9999999};
				int itmpz0dphi[4] = {9999999,9999999,9999999,9999999};
				int itmpz0dzordr[4] = {9999999,9999999,9999999,9999999};
				
				// check if current hit pattern is possible for the seeding pairs
				bool goodseed = false;
				unsigned int iseed1, iseed2;
				if (seed == "L1L2") {
					goodseed = (layerhits&(1<<5)) and (layerhits&(1<<4));
					iseed1 = 0;
					iseed2 = 1;
				}
				else if (seed == "L3L4") {
					goodseed = (layerhits&(1<<3)) and (layerhits&(1<<2));
					iseed1 = 2;
					iseed2 = 3;
				}
				else if (seed == "L5L6") {
					goodseed = (layerhits&(1<<1)) and (layerhits&(1<<0));
					iseed1 = 4;
					iseed2 = 5;
				}
				else if (seed == "F1F2" or seed == "B1B2") {
					goodseed = (diskmask&(1<<9)) and (diskmask&(1<<7));
					iseed1 = 6;
					iseed2 = 7;
				}
				else if (seed == "F3F4" or seed == "B3B4") {
					goodseed = (diskmask&(1<<5)) and (diskmask&(1<<3));
					iseed1 = 8;
					iseed2 = 9;
				}
				else if (seed == "L1F1" or seed == "L1B1") {
					goodseed = (layerhits&(1<<5)) and (diskmask&(1<<9));
					iseed1 = 0;
					iseed2 = 6;
				}
				else {
					cerr << "Not valid seeding pairs" << endl;
					return;
				}
				
				// loop over each bit in the hit pattern
				int ider = 0;
				int inputI = 99;
				if (goodseed) {
					for (unsigned int ihit = 0; ihit < 11; ++ihit) {
						// skip seeding layers
						if (ihit == iseed1 or ihit == iseed2) {
							ider++;
							continue;
						}
						
						if ( hits&(1<<(10-ihit)) ) {
							//
							if (seed == "L1L2") {
								if (ihit == 2 or ihit == 9) inputI = 0;//L3 or D4
								if (ihit == 3 or ihit == 8) inputI = 1;//L4 or D3
								if (ihit == 4 or ihit == 7) inputI = 2;//L5 or D2
								if (ihit == 5 or ihit == 6) inputI = 3;//L6 or D1
							}
							else if (seed == "L3L4") {
								if (ihit == 0) inputI = 0; //L1
								if (ihit == 1) inputI = 1; //L2
								if (ihit == 4 or ihit == 7) inputI = 2;//L5 or D2
								if (ihit == 5 or ihit == 6) inputI = 3;//L6 or D1
							}
							else if (seed == "L5L6") {
								if (ihit == 0) inputI = 0; //L1
								if (ihit == 1) inputI = 1; //L2
								if (ihit == 2) inputI = 2; //L3
								if (ihit == 3) inputI = 3; //L4
							}
							else if (seed == "F1F2" or seed == "B1B2") {
								if (ihit == 0) inputI = 0; //L1
								if (ihit == 8) inputI = 1; //D3
								if (ihit == 9) inputI = 2; //D4
								if (ihit == 1 or ihit == 10) inputI = 3;//L2 or D5
							}
							else if (seed == "F3F4" or seed == "B3B4") {
								if (ihit == 0) inputI = 0; //L1
								if (ihit == 6) inputI = 1; //D1
								if (ihit == 7) inputI = 2; //D2
								if (ihit == 1 or ihit == 10) inputI = 3;//L2 or D5
							}
							else if (seed == "L1F1" or seed == "L1B1") {
								if (ihit == 7) inputI = 0; //D2
								if (ihit == 8) inputI = 1; //D3
								if (ihit == 9) inputI = 2; //D4
								if (ihit == 10) inputI = 3; //D5
							}

							itmprinvdphi[inputI] = der.getirinvdphi(ider);
							itmprinvdzordr[inputI] = der.getirinvdzordr(ider);
							itmpphi0dphi[inputI] = der.getiphi0dphi(ider);
							itmpphi0dzordr[inputI] = der.getiphi0dzordr(ider);
							itmptdphi[inputI] = der.getitdphi(ider);
							itmptdzordr[inputI] = der.getitdzordr(ider);
							itmpz0dphi[inputI] = der.getiz0dphi(ider);
							itmpz0dzordr[inputI] = der.getiz0dzordr(ider);
							
							ider++;
						} // end if (hits&(1<<(10-ihit)))
						
					} // end of hit loop
				} // end if (goodseed)
				
				FPGAWord tmprinvdphi[4];
				for (unsigned int j = 0; j < 4; ++j) {
					if (itmprinvdphi[j]>(1<<13)) itmprinvdphi[j]=(1<<13)-1;
					tmprinvdphi[j].set(itmprinvdphi[j],14,false,__LINE__,__FILE__);
				}
				outrinvdphi[i] << tmprinvdphi[0].str() << tmprinvdphi[1].str()
							   << tmprinvdphi[2].str() << tmprinvdphi[3].str()
							   << endl;
				
				FPGAWord tmprinvdzordr[4];
				for (unsigned int j = 0; j < 4; ++j) {
					if (itmprinvdzordr[j]>(1<<15)) itmprinvdzordr[j]=(1<<15)-1;
					tmprinvdzordr[j].set(itmprinvdzordr[j],16,false,__LINE__,__FILE__);
				}
				outrinvdzordr[i] << tmprinvdzordr[0].str() << tmprinvdzordr[1].str()
								 << tmprinvdzordr[2].str() << tmprinvdzordr[3].str()
								 << endl;
				
				FPGAWord tmpphi0dphi[4];
				for (unsigned int j = 0; j < 4; ++j) {
					if (itmpphi0dphi[j]>(1<<13)) itmpphi0dphi[j]=(1<<13)-1;
					tmpphi0dphi[j].set(itmpphi0dphi[j],14,false,__LINE__,__FILE__);
				}
				outphi0dphi[i] << tmpphi0dphi[0].str() << tmpphi0dphi[1].str()
							   << tmpphi0dphi[2].str() << tmpphi0dphi[3].str()
							   << endl;
				
				FPGAWord tmpphi0dzordr[4];
				for (unsigned int j = 0; j < 4; ++j) {
					if (itmpphi0dzordr[j]>(1<<15)) itmpphi0dzordr[j]=(1<<15)-1;
					tmpphi0dzordr[j].set(itmpphi0dzordr[j],16,false,__LINE__,__FILE__);
				}
				outphi0dzordr[i] << tmpphi0dzordr[0].str() << tmpphi0dzordr[1].str()
								 << tmpphi0dzordr[2].str() << tmpphi0dzordr[3].str()
								 << endl;
				
				FPGAWord tmptdphi[4];
				for (unsigned int j = 0; j < 4; ++j) {
					if (itmptdphi[j]>(1<<13)) itmptdphi[j]=(1<<13)-1;
					tmptdphi[j].set(itmptdphi[j],14,false,__LINE__,__FILE__);
				}
				outtdphi[i] << tmptdphi[0].str() << tmptdphi[1].str()
							<< tmptdphi[2].str() << tmptdphi[3].str()
							<< endl;
				
				FPGAWord tmptdzordr[4];
				for (unsigned int j = 0; j < 4; ++j) {
					if (itmptdzordr[j]>(1<<15)) itmptdzordr[j]=(1<<15)-1;
					tmptdzordr[j].set(itmptdzordr[j],16,false,__LINE__,__FILE__);
				}
				outtdzordr[i] << tmptdzordr[0].str() << tmptdzordr[1].str()
							  << tmptdzordr[2].str() << tmptdzordr[3].str()
							  << endl;
			
				FPGAWord tmpz0dphi[4];
				for (unsigned int j = 0; j < 4; ++j) {
					if (itmpz0dphi[j]>(1<<13)) itmpz0dphi[j]=(1<<13)-1;
					tmpz0dphi[j].set(itmpz0dphi[j],14,false,__LINE__,__FILE__);
				}
				outz0dphi[i] << tmpz0dphi[0].str() << tmpz0dphi[1].str()
							 << tmpz0dphi[2].str() << tmpz0dphi[3].str()
							 << endl;
				
				FPGAWord tmpz0dzordr[4];
				for (unsigned int j = 0; j < 4; ++j) {
					if (itmpz0dzordr[j]>(1<<15)) itmpz0dzordr[j]=(1<<15)-1;
					tmpz0dzordr[j].set(itmpz0dzordr[j],16,false,__LINE__,__FILE__);
				}
				outz0dzordr[i] << tmpz0dzordr[0].str() << tmpz0dzordr[1].str()
							   << tmpz0dzordr[2].str() << tmpz0dzordr[3].str()
							   << endl;
								
			} // end of seeding loop
			
			
		} // end of derivatives_ loop

		// close files
		for (unsigned int i = 0; i < 6; ++i) {
			outrinvdphi[i].close();
			outrinvdzordr[i].close();
			outphi0dphi[i].close();
			outphi0dzordr[i].close();
			outtdphi[i].close();
			outtdzordr[i].close();
			outz0dphi[i].close();
			outz0dzordr[i].close();
		}
		
	} // end of if(writeFitDerTable)
	// --------------------------------------------------------

  }


  void invert(double M[4][8],unsigned int n){

    assert(n<=4);

    unsigned int i,j,k;
    double ratio,a;

    for(i = 0; i < n; i++){
      for(j = n; j < 2*n; j++){
	if(i==(j-n))
	  M[i][j] = 1.0;
	else
	  M[i][j] = 0.0;
      }
    }

    for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){
	if(i!=j){
	  ratio = M[j][i]/M[i][i];
	  for(k = 0; k < 2*n; k++){
	    M[j][k] -= ratio * M[i][k];
	  }
	}
      }
    }

    for(i = 0; i < n; i++){
      a = M[i][i];
      for(j = 0; j < 2*n; j++){
	M[i][j] /= a;
      }
    }
  }



  void calculateDerivativesTable(unsigned int nlayers,
				 double r[6],
				 unsigned int ndisks,
				 double z[5],
				 double alpha[5],
				 double t,
				 double rinv,
				 double D[4][12],
				 double MinvDt[4][12],
				 int iMinvDt[4][12]){


    double sigmax=0.01/sqrt(12.0);
    double sigmaz=0.15/sqrt(12.0);
    double sigmaz2=5.0/sqrt(12.0);

    unsigned int n=nlayers+ndisks;
    
    assert(n<=6);

    double rnew[6];

    int j=0;

    double M[4][8];


    //here we handle a barrel hit
    for(unsigned int i=0;i<nlayers;i++) {

      double ri=r[i];

      rnew[i]=ri;

      //double rinv=0.000001; //should simplify?

      //first we have the phi position
      D[0][j]=-0.5*ri*ri/sqrt(1-0.25*ri*ri*rinv*rinv)/sigmax;
      D[1][j]=ri/sigmax;
      D[2][j]=0.0;
      D[3][j]=0.0;
      j++;
      //second the z position
      D[0][j]=0.0;
      D[1][j]=0.0;
      if (ri<60.0) {
	D[2][j]=(2/rinv)*asin(0.5*ri*rinv)/sigmaz;
	D[3][j]=1.0/sigmaz;
      } else {
	D[2][j]=(2/rinv)*asin(0.5*ri*rinv)/sigmaz2;
	D[3][j]=1.0/sigmaz2;
      }

      //cout << "0 D "<<ri<<" "<<rinv<<" "
      //   <<D[0][j]<<" "<<D[1][j]<<" "<<D[2][j]<<" "<<D[3][j]<<endl;

      j++;

    }


    for(unsigned int i=0;i<ndisks;i++) {

      double zi=z[i];


      //double rinv=0.000001;
      double z0=0.0;
 
      double rmultiplier=alpha[i]*zi/t;
      //rmultiplier=0.0;
      double phimultiplier=zi/t;
      
      //cout << "i zi t : "<<i<<" "<<zi<<" "<<t<<endl;

      double drdrinv=-2.0*sin(0.5*rinv*(zi-z0)/t)/(rinv*rinv)
      +(zi-z0)*cos(0.5*rinv*(zi-z0)/t)/(rinv*t);
      double drdphi0=0;
      double drdt=-(zi-z0)*cos(0.5*rinv*(zi-z0)/t)/(t*t);
      double drdz0=-cos(0.5*rinv*(zi-z0)/t)/t;


      double dphidrinv=-0.5*(zi-z0)/t;
      double dphidphi0=1.0;
      double dphidt=0.5*rinv*(zi-z0)/(t*t);
      double dphidz0=0.5*rinv/t;

      double r=(zi-z0)/t;

      //cout << "r alpha : "<<r<<" "<<alpha[i]<<endl;

      rnew[i+nlayers]=r;

      //cout << "rnew["<<i+nlayers<<"] = "<<rnew[i+nlayers]
      //	   <<" "<<zi<<" "<<z0<<" "<<t<<endl;


      //cout << "FITLINNEW r = "<<r<<endl; 

      //second the rphi position
      double sigmaxtmp=sigmax;
      if (fabs(alpha[i])>1e-10) {
	//cout << "Table for outer disks : "<<r<<" "<<zi<<endl;
	sigmaxtmp*=errfac;
      }

      D[0][j]=(phimultiplier*dphidrinv+rmultiplier*drdrinv)/sigmaxtmp;
      D[1][j]=(phimultiplier*dphidphi0+rmultiplier*drdphi0)/sigmaxtmp;
      D[2][j]=(phimultiplier*dphidt+rmultiplier*drdt)/sigmaxtmp;
      D[3][j]=(phimultiplier*dphidz0+rmultiplier*drdz0)/sigmaxtmp;

      //cout << "1 D "<<D[0][j]<<" "<<D[1][j]<<" "<<D[2][j]<<" "<<D[3][j]<<endl;

      j++;

      //cout << "alpha i r = "<<alpha[i]<<" "<<i<<" "<<r<<endl;
      if (fabs(alpha[i])<1e-10) {
	D[0][j]=drdrinv/sigmaz;
	D[1][j]=drdphi0/sigmaz;
	D[2][j]=drdt/sigmaz;
	D[3][j]=drdz0/sigmaz;
      }
      else {
	D[0][j]=drdrinv/sigmaz2;
	D[1][j]=drdphi0/sigmaz2;
	D[2][j]=drdt/sigmaz2;
	D[3][j]=drdz0/sigmaz2;
      }

      //cout << "2 D "<<D[0][j]<<" "<<D[1][j]<<" "<<D[2][j]<<" "<<D[3][j]<<endl;

      j++;
      

    }

    //cout << "j n : "<<j<<" "<<n<<endl;
    
    for(unsigned int i1=0;i1<4;i1++){
      for(unsigned int i2=0;i2<4;i2++){
	M[i1][i2]=0.0;
	for(unsigned int j=0;j<2*n;j++){
	  M[i1][i2]+=D[i1][j]*D[i2][j];	  
	}
      }
    }

    invert(M,4);

    for(unsigned int j=0;j<12;j++) {
      for(unsigned int i1=0;i1<4;i1++) {
	MinvDt[i1][j]=0.0;
	iMinvDt[i1][j]=0;
      }
    }  

    for(unsigned int j=0;j<2*n;j++) {
      for(unsigned int i1=0;i1<4;i1++) {
	for(unsigned int i2=0;i2<4;i2++) {
	  MinvDt[i1][j]+=M[i1][i2+4]*D[i2][j];
	  //cout << "0 M D = "<<M[i1][i2+4]<<" "<<D[i2][j]<<endl;
	
	}
      }
    }


    for (unsigned int i=0;i<n;i++) {


      //First the barrel
      if (i<nlayers) {

	MinvDt[0][2*i]*=rnew[i]/sigmax;
	MinvDt[1][2*i]*=rnew[i]/sigmax;
	MinvDt[2][2*i]*=rnew[i]/sigmax;
	MinvDt[3][2*i]*=rnew[i]/sigmax;

	//cout << "1 MinvDt[0][2*i] = "<<MinvDt[0][2*i]<<endl;

      
	iMinvDt[0][2*i]=(1<<fitrinvbitshift)*MinvDt[0][2*i]*kphi1/krinvpars;
	iMinvDt[1][2*i]=(1<<fitphi0bitshift)*MinvDt[1][2*i]*kphi1/kphi0pars;
	iMinvDt[2][2*i]=(1<<fittbitshift)*MinvDt[2][2*i]*kphi1/ktpars;
	iMinvDt[3][2*i]=(1<<fitz0bitshift)*MinvDt[3][2*i]*kphi1/kzpars;

	if (rnew[i]<60.0) {
	  MinvDt[0][2*i+1]/=sigmaz;
	  MinvDt[1][2*i+1]/=sigmaz;
	  MinvDt[2][2*i+1]/=sigmaz;
	  MinvDt[3][2*i+1]/=sigmaz;

	  iMinvDt[0][2*i+1]=(1<<fitrinvbitshift)*MinvDt[0][2*i+1]*kzproj/krinvpars;
	  iMinvDt[1][2*i+1]=(1<<fitphi0bitshift)*MinvDt[1][2*i+1]*kzproj/kphi0pars;
	  iMinvDt[2][2*i+1]=(1<<fittbitshift)*MinvDt[2][2*i+1]*kzproj/ktpars;
	  iMinvDt[3][2*i+1]=(1<<fitz0bitshift)*MinvDt[3][2*i+1]*kzproj/kzpars;
	} else {
	  MinvDt[0][2*i+1]/=sigmaz2;
	  MinvDt[1][2*i+1]/=sigmaz2;
	  MinvDt[2][2*i+1]/=sigmaz2;
	  MinvDt[3][2*i+1]/=sigmaz2;

	  iMinvDt[0][2*i+1]=(1<<fitrinvbitshift)*MinvDt[0][2*i+1]*kzproj/krinvpars;
	  iMinvDt[1][2*i+1]=(1<<fitphi0bitshift)*MinvDt[1][2*i+1]*kzproj/kphi0pars;
	  iMinvDt[2][2*i+1]=(1<<fittbitshift)*MinvDt[2][2*i+1]*kzproj/ktpars;
	  iMinvDt[3][2*i+1]=(1<<fitz0bitshift)*MinvDt[3][2*i+1]*kzproj/kzpars;
	}
      }

      //Secondly the disks
      else {


	if (fabs(alpha[i-nlayers])<1e-10) {
	  MinvDt[0][2*i]*=(rnew[i]/sigmax);
	  MinvDt[1][2*i]*=(rnew[i]/sigmax);
	  MinvDt[2][2*i]*=(rnew[i]/sigmax);
	  MinvDt[3][2*i]*=(rnew[i]/sigmax);
	} else {
	  MinvDt[0][2*i]*=(rnew[i]/(errfac*sigmax));
	  MinvDt[1][2*i]*=(rnew[i]/(errfac*sigmax));
	  MinvDt[2][2*i]*=(rnew[i]/(errfac*sigmax));
	  MinvDt[3][2*i]*=(rnew[i]/(errfac*sigmax));
	}      

	//cout << "2 MinvDt[0][2*i] = "<<MinvDt[0][2*i]<<endl;

	assert(MinvDt[0][2*i]==MinvDt[0][2*i]);

	iMinvDt[0][2*i]=(1<<fitrinvbitshift)*MinvDt[0][2*i]*kphiprojdisk/krinvparsdisk;
	iMinvDt[1][2*i]=(1<<fitphi0bitshift)*MinvDt[1][2*i]*kphiprojdisk/kphi0parsdisk;
	iMinvDt[2][2*i]=(1<<fittbitshift)*MinvDt[2][2*i]*kphiprojdisk/ktparsdisk;
	iMinvDt[3][2*i]=(1<<fitz0bitshift)*MinvDt[3][2*i]*kphiprojdisk/kzdisk;

	if (alpha[i-nlayers]==0.0) {
	  MinvDt[0][2*i+1]/=sigmaz;
	  MinvDt[1][2*i+1]/=sigmaz;
	  MinvDt[2][2*i+1]/=sigmaz;
	  MinvDt[3][2*i+1]/=sigmaz;
	} else {
	  MinvDt[0][2*i+1]/=sigmaz2;
	  MinvDt[1][2*i+1]/=sigmaz2;
	  MinvDt[2][2*i+1]/=sigmaz2;
	  MinvDt[3][2*i+1]/=sigmaz2;
	}

	iMinvDt[0][2*i+1]=(1<<fitrinvbitshift)*MinvDt[0][2*i+1]*krprojshiftdisk/krinvparsdisk;
	iMinvDt[1][2*i+1]=(1<<fitphi0bitshift)*MinvDt[1][2*i+1]*krprojshiftdisk/kphi0parsdisk;
	iMinvDt[2][2*i+1]=(1<<fittbitshift)*MinvDt[2][2*i+1]*krprojshiftdisk/ktparsdisk;
	iMinvDt[3][2*i+1]=(1<<fitz0bitshift)*MinvDt[3][2*i+1]*krprojshiftdisk/kzdisk;
      
      }

    }
    

  }
  

  double gett(int diskmask) { //should use layers also..

    if (diskmask==0) return 0.0;

    double tmax=1000.0;
    double tmin=0.0;

    for(int d=1;d<=5;d++) {

      if (diskmask&(1<<(2*(5-d)+1))) { //PS hit
	double dmax=zmean[d-1]/22.0;
	if (dmax>sinh(2.4)) dmax=sinh(2.4);
	double dmin=zmean[d-1]/60.0;
	if (dmax<tmax) tmax=dmax;
	if (dmin>tmin) tmin=dmin;
      } 

      if (diskmask&(1<<(2*(5-d)))) { //2S hit
	double dmax=zmean[d-1]/60.0;
	double dmin=zmean[d-1]/105.0;
	if (dmax<tmax) tmax=dmax;
	if (dmin>tmin) tmin=dmin;	
      } 

    }

    return 0.5*(tmax+tmin);

  }


private:


  vector<int> LayerMem_;
  vector<int> DiskMem_;

  vector<int> LayerDiskMem_;

  unsigned int LayerMemBits_;
  unsigned int DiskMemBits_;
  unsigned int LayerDiskMemBits_;
  unsigned int alphaBits_;

  unsigned int Nlay_;
  unsigned int Ndisk_;

  vector<FPGATrackDer> derivatives_;
  
  int nextLayerValue_;
  int nextDiskValue_;
  int nextLayerDiskValue_;
  int lastMultiplicity_;

};



#endif




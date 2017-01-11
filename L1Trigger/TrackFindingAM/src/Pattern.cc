#include "../interface/Pattern.h"
#include "../interface/Detector.h"

Pattern::Pattern(int nb){
  nb_layer = nb;
  strips=NULL;
  nb_strips=NULL;
  nbFakeSSKnown=false;
  nbFakeSS=0;
  order_in_chip = -1;
  for(int i=0;i<nb_layer;i++){
    layer_strips.push_back(NULL);
  }
}

Pattern::Pattern(){
  nb_layer = 3;
  strips=NULL;
  nb_strips=NULL;
  nbFakeSSKnown=false;
  nbFakeSS=0;
  order_in_chip = -1;
  for(int i=0;i<nb_layer;i++){
    layer_strips.push_back(NULL);
  }
}

Pattern::Pattern(const Pattern& p){
  nb_layer = p.nb_layer;
  nb_strips = NULL;
  strips = NULL;
  nbFakeSSKnown=false;
  nbFakeSS=0;
  order_in_chip = p.order_in_chip;

  if(p.strips!=NULL){
    nb_strips=new char[nb_layer];
    strips = new SuperStrip**[nb_layer];
    for(int i=0;i<nb_layer;i++){
      strips[i] = new SuperStrip*[p.nb_strips[i]];
    }
  }
  for(int i=0;i<nb_layer;i++){
    layer_strips.push_back(p.layer_strips[i]->clone());
    if(p.strips!=NULL){
      nb_strips[i] = p.nb_strips[i];
      for(int j=0;j<nb_strips[i];j++){
	strips[i][j]=p.strips[i][j];
      }
    }
  }
}

Pattern::~Pattern(){
  if(nb_strips!=NULL){
    delete [] nb_strips;
  }
  for(int i=0;i<nb_layer;i++){
    if(layer_strips[i]!=NULL)
      delete layer_strips[i];
    if(strips!=NULL)
      delete [] strips[i];
  }
  delete [] strips;
  strips = NULL;
  nb_strips = NULL;
}

void Pattern::setLayerStrip(int layer, PatternLayer* strip){
  if(layer<nb_layer && layer>=0){
    layer_strips[layer]=strip->clone();
  }
}

PatternLayer* Pattern::getLayerStrip(int layer){
  if(layer<nb_layer && layer>=0)
    return layer_strips[layer];
  return NULL;
}

PatternLayer* Pattern::operator[](int layer){
  if(layer<nb_layer && layer>=0)
    return layer_strips[layer];
  return NULL;
}

bool Pattern::isActive(int active_threshold){
  int score=0;
  if(strips!=NULL){
    for(int i=0;i<nb_layer;i++){
      for(int j=0;j<(int)nb_strips[i];j++){
	if(strips[i][j]->isHit()){
	  score++;
	  break;
	}
      }
    }
  }
  if(score>=active_threshold)
    return true;
  else
    return false;
}

bool Pattern::isActiveUsingMissingHit(int nb_allowed_missing_hit, int active_threshold){
  int score=0;
  if(strips!=NULL){
    for(int i=0;i<nb_layer;i++){
      for(int j=0;j<(int)nb_strips[i];j++){
	if(strips[i][j]->isHit()){
	  score++;
	  break;
	}
      }
    }
  }
  int limit=nb_layer-getNbFakeSuperstrips()-nb_allowed_missing_hit;
  if(limit<active_threshold)
    limit=active_threshold;
  if(score>=limit){
    return true;
  }
  else
    return false;
}

string Pattern::getKey() const{
  string key="";
  for(int i=0;i<nb_layer;i++){
    key.append(layer_strips[i]->getCode());
  }
  return key;
}

int Pattern::getNbLayers() const {
  return nb_layer;
}

void Pattern::unlink(){ 
  if(nb_strips!=NULL)
    delete [] nb_strips;
  
  if(strips!=NULL){
    for(int i=0;i<nb_layer;i++){
      delete [] strips[i];
    }
  }

  delete [] strips;

  strips = NULL;
  nb_strips = NULL;
}

void Pattern::link(Detector& d){
  if(strips!=NULL){// already linked
    unlink();
  }
  
  nb_strips=new char[nb_layer];
  strips = new SuperStrip**[nb_layer];
  for(int i=0;i<nb_layer;i++){
    vector<SuperStrip*> tmp_strips = layer_strips[i]->getSuperStrip(i, d);
    nb_strips[i] = tmp_strips.size();
    strips[i] = new SuperStrip*[nb_strips[i]];
    for(unsigned int j=0;j<tmp_strips.size();j++){
      if(tmp_strips[j]==NULL)
	cout<<"linking to a NULL SuperStrip!!"<<endl;
      strips[i][j] = tmp_strips[j];
    }
  }
}

#ifdef IPNL_USE_CUDA
void Pattern::linkCuda(patternBank* p, deviceDetector* d, int pattern_index, const vector< vector<int> >& sec, const vector<map<int, vector<int> > >& modules, vector<int> layers,
                       unsigned int* cache){
  memset(cache,PATTERN_UNUSED,PATTERN_LAYERS*PATTERN_SSTRIPS*sizeof(unsigned int));
  for(int i=0;i<nb_layer;i++){
    layer_strips[i]->getSuperStripCuda(i, sec[i], modules[i], layers[i], cache+i*PATTERN_SSTRIPS);
  }
  cudaSetLink(p,pattern_index*PATTERN_SIZE,cache);
}
#endif

bool compareHits(Hit* h1, Hit* h2){
  return h1->getID()<h2->getID();
}

vector<Hit*> Pattern::getHits(){
  vector<Hit*> hits;
  for(int i=0;i<nb_layer;i++){
    vector<Hit*> local_hits;
    for(int j=0;j<(int)nb_strips[i];j++){
      vector<Hit*> l_hits = strips[i][j]->getHits();
      local_hits.insert(local_hits.end(), l_hits.begin(), l_hits.end());
    }
    sort(local_hits.begin(),local_hits.end(),compareHits);
    int nb_hits=(int)local_hits.size();
    if(Detector::hw_limit_stub_per_layer>0 && nb_hits>Detector::hw_limit_stub_per_layer)
      nb_hits=Detector::hw_limit_stub_per_layer;
    hits.insert(hits.end(),local_hits.begin(),local_hits.begin()+nb_hits);
  }
  return hits;
}

vector<Hit*> Pattern::getHits(int layerPosition){
  vector<Hit*> hits;
  if(layerPosition>=0 && layerPosition<nb_layer){
    for(int j=0;j<(int)nb_strips[layerPosition];j++){
      vector<Hit*> l_hits = strips[layerPosition][j]->getHits();
      hits.insert(hits.end(), l_hits.begin(), l_hits.end());
    }
  }
  sort(hits.begin(),hits.end(),compareHits);
  if(Detector::hw_limit_stub_per_layer>0 && (int)hits.size()>Detector::hw_limit_stub_per_layer){ 
    hits.erase (hits.begin()+Detector::hw_limit_stub_per_layer,hits.end());
  }
  return hits;
}

bool Pattern::contains(Pattern* hdp){
  if(nb_layer!=hdp->getNbLayers())
    return false;
  
  for(int i=0;i<nb_layer;i++){
    int base_index = layer_strips[i]->getStripCode()<<layer_strips[i]->getDCBitsNumber();
    vector<short> positions=layer_strips[i]->getPositionsFromDC();
    bool found = false;
    short reference = hdp->getLayerStrip(i)->getStripCode();
    if(positions.size()==0){//no DC bits used
      if(reference==base_index){
	found=true;
      }
    }
    else{
      for(unsigned int j=0;j<positions.size();j++){
	int index = base_index | positions[j];
	if(reference==index){
	  found=true;
	  break;
	}
      }
    }
    if(!found){
      return false;
    }
  }
  return true;
}

ostream& operator<<(ostream& out, const Pattern& s){
  out<<s.getOrderInChip()<<" : ";
  for(int i=0;i<s.getNbLayers();i++){
    out<<s.layer_strips[i]->toString()<<" - ";
  }
  out<<endl;
  return out;
}

int Pattern::getNbFakeSuperstrips(){
  if(nbFakeSSKnown)
    return nbFakeSS;
  else{
    int score = 0;
    for(int k=0;k<getNbLayers();k++){
      PatternLayer* mp = getLayerStrip(k);
      if(mp->isFake())
	score++;
    }
    nbFakeSSKnown=true;
    nbFakeSS=score;
    return score;
  }
}

void Pattern::setOrderInChip(int i){
  order_in_chip = i;
}

int Pattern::getOrderInChip() const{
  return order_in_chip;
}

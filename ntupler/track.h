#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

struct Particle{
 private:
  double pt;
  double eta;
  double phi;
 public:
  Particle(double pt, double eta, double phi){
   this->pt=pt;
   this->eta=eta;
   this->phi=phi;
  }
};

class track{
 private:
  double pt;
  double eta;
  double phi;
  int nmatched;
  double matching_dR;
  particle **matched_particles;
 public:
  track(double pt, double eta, double phi){
   matching_dR=0.1;
   nmatched=0;
   this->pt=pt;
   this->eta=eta;
   this->phi=phi;   
  }
  void add_matched_particle(Particle *particle){
   if(sqrt(pow(particle->eta-eta,2)+pow(particle->phi-phi,2)<matching_dR){
    if(matched_particles==NULL){
     matched_particles = (Particle **)calloc(2,sizeof(Particle *));
     matched_particles[0]=particle;
     matched_particles[1]=NULL;
    }
    else
    {
     matched_particles = (Particle **)realloc(matched_particles,(nmatched+2)*sizeof(Particle *));
     matched_particles[nmatched]=particle;
     matched_particles[nmatched+1]=NULL;
    }
    nmatched++;
   }
  }
};


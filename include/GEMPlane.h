#ifndef GEMPlane_H
#define GEMPlane_H
/*This is for a plane (either X/Y or U/V
 *
 *
 *
 */
#include <iostream>
#include <vector>
#include <map>
#include <TROOT.h>
//#include "TH1D.h"
//#include "TRandom3.h"

class GEMPlane {
 public:
  GEMPlane();
  GEMPlane(short layer, short mod, int nstrips, int nsamples = 6, double thr = 100, double offset = 0, double roangle = 0);
  virtual ~GEMPlane();
  void Clear();
  
  //void SetStripThreshold(double thr){Striphr = ;};

  Short_t Layer(){return fLayer;};
  Short_t Module(){return fModule;};
  Double_t dX(){return fdX;};
  Double_t Xoffset(){return fXoffset;};
  Double_t ROangle(){return fROangle;};
  Int_t GetNStrips(){return fNStrips;};
  Short_t GetADC(int strip, int samp){return fStripADC[strip*fNSamples+samp];};
  Int_t GetADCSum(int strip){return fStripADCsum[strip];};
  void SetADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripADCsum[strip]+= adc-fStripADC[strip*fNSamples+samp];// CG: Not clear why to subtract the 0 vector (as this stage the fStripADC is set to zero. It only set to adc value in the next line. However, it gives the correct sum
      fStripADC[strip*fNSamples+samp] = adc;
    }
  };
  void AddADC(int strip, int samp, int adc){
    if(strip<fNStrips){
      fStripADC[strip*fNSamples+samp]+=adc; //CG: why do we need to do here - not clear!
      fStripADCsum[strip]+=adc;
    }//else{
    //printf("strip = %d / %d, sample %d /%d \n", strip, fNStrips, samp, fNSamples);
    //}
  };
  
 private:
  // ADC sampled value of strip array of each axis
  Short_t fLayer;
  Short_t fModule;
  Int_t fNStrips;
  Int_t fNSamples;
  Double_t fStripThr;//threshold for ADC sum
  Int_t* fStripADCsum;
  Short_t* fStripADC;
  double fdX;
  double fXoffset;
  double fROangle;
  //ClassDef(GEMPlane, 1)
};
#endif

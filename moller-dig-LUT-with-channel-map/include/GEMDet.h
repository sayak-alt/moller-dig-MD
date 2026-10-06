#ifndef GEMDet_H
#define GEMDet_H

#include <iostream>
#include <vector>
#include <map>
#include "GEMPlane.h"

//________________________________
class GEMDet {
 public:
  GEMDet();
  GEMDet(UShort_t uinqueid, UInt_t nplanes, int* layer, int* nstrips, double* offset, double* roangle, int nsamp, double zsup_thr);
  virtual ~GEMDet();
  void Clear();

  struct gemhit{
    int source;
    int module;
    double edep;
    //double tmin;
    //double tmax;
    double t;
    double xin;
    double yin;
    double zin;
    double xout;
    double yout;
    double zout;
  };

  std::vector<gemhit> fGEMhits;
  
  //private:
  UShort_t fUniqueID;
  UInt_t fNPlanes;
  Double_t fGateWidth;
  std::vector<Double_t> fZLayer;
  //std::map<int, GEMPlane> GEMPlanes;
  std::vector<GEMPlane> GEMPlanes;
};

#endif

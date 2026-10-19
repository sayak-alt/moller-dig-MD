#include "GEMDet.h"
#include "TMath.h"
#define DEBUG 0

using namespace std;

GEMDet::GEMDet()
{
}

GEMDet::GEMDet(UShort_t uniqueid, UInt_t nplanes, int* layer, int* nstrips, double* offset, double* roangle, int nsamp, double zsup_thr):
  fUniqueID(uniqueid), fNPlanes(nplanes)
{
  //for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i] = DigGEMPlane(nstrips[i], nsamp, zsup_thr);
  for(uint i = 0; i<fNPlanes; i++){
    cout <<"Plane Number:  " << ", Layer # " << ", strip = " << ", Offset =  " << ", Readout Angle = " << endl; 
    cout << i << " " << layer[i] << " " << nstrips[i] << " " << offset[i] << " " << roangle[i] << endl; 
    GEMPlanes.push_back(GEMPlane(layer[i], i/2, nstrips[i], nsamp, zsup_thr, offset[i], roangle[i]));
  }
}

GEMDet::~GEMDet()
{
  
}

void GEMDet::Clear()
{
  for(uint i = 0; i<fNPlanes; i++)GEMPlanes[i].Clear();
  fGEMhits.clear();
}

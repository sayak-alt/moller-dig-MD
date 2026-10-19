#include "GEMPlane.h"
#include "TMath.h"

using namespace std;

//
// Class GEMPlane
//
GEMPlane::GEMPlane() :
  fLayer(0), fNStrips(605), fNSamples(6), fStripThr(100), fXoffset(0), fROangle(0)
{
}

GEMPlane::GEMPlane(short layer, short mod, int nstrips, int nsamples, double thr, double offset, double roangle) :
  fLayer(layer), fModule(mod), fNStrips(nstrips), fNSamples(nsamples), fStripThr(thr), fXoffset(offset), fROangle(roangle)
{
  //fModule = mod;
  fdX = fNStrips*8.e-4;
  
  fStripADCsum = new Int_t[fNStrips];
  fStripADC = new Short_t[fNStrips*fNSamples];
  Clear();
}

GEMPlane::~GEMPlane()
{
  Clear();
}


void GEMPlane::Clear()
{
  memset(fStripADCsum, 0, fNStrips*sizeof(Int_t));
  memset(fStripADC, 0, fNStrips*fNSamples*sizeof(Short_t));
}

//ClassImp(GEMPlane);


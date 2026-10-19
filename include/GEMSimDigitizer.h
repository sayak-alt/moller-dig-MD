#ifndef GEMSIMDIGITIZER_H
#define GEMSIMDIGITIZER_H

#include "TRandom3.h"
#include "TVector3.h"
//#include "GEMPlane.h"
//#include "GEMDet.h"
//#include "TArrayS.h"
//#include "TArrayI.h"
//#include "TArrayD.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>
#include <vector>
#include <chrono>

#define DBG_HISTOS 0

class GEMDet;
class GEMPlane;
//class gmn_tree;
class remoll_tree;

class GEMSimDigitizer {
 public:
  //Constructor and destructor
  GEMSimDigitizer();
  GEMSimDigitizer(int nchambers, double* trigoffset, double* gain, double zsup_thr, int napv = 0, double* commonmode_array = 0);
  virtual ~GEMSimDigitizer();
  void Print();
  
  Int_t Digitize (GEMDet* gemdet, TRandom3* R, bool bkgdonly = false);//, gmn_tree* T);
  //void CheckOut(GEMDet* gemdet, TRandom3* R, gmn_tree* T);
  void CheckOut(GEMDet* gemdet, const int uniqueid, TRandom3* R, remoll_tree* T, bool sigonly = false);
  //void FillBBGEMTree(const GEMPlane pl, gmn_tree* T, int j);
  void write_histos();
  void print_time_execution();
  
  struct IonPar_t {
    Double_t X;       // position of the point on the projection
    Double_t Y;
    Double_t Charge;  // Charge deposited by this ion
    Double_t SNorm;   // 3 x radius of ion diffusion area at readout
    Double_t R2;      // = SNorm^2 : radius of numerical integration area
    Double_t ggnorm;  // = Charge/R2/pi : charge per unit area
  };

  //private:
  void AvaModel(const int ic, //module number
		GEMDet* gemdet,
		TRandom3* R,
		const TVector3& xi,
		const TVector3& xo,
		const Double_t t0);
  
  void AvaModel_2(const int ic, //module number
		  GEMDet* gemdet,
		  TRandom3* R,
		  const TVector3& xi,
		  const TVector3& xo,
		  const Double_t t0);
  
  void Integration_semiana(double roangle, 
			    double xl, double xr, double yb, double yt, 
			    int nx, double xbw);
  
  void Integration_fastappx(TRandom3* R, double roangle, 
			    double nions_strip,
			    double xc_hit, double yc_hit,
			    double dx2_hit, double dy2_hit,
			    double xl, double xr, double yb, double yt, 
			    int nx, double xbw, int ny, double ybw);
  
  void IonModel (TRandom3* R,
		 const TVector3& xi,
		 const TVector3& xo,
		 const Double_t elost);
  
  std::vector<Double_t> fTriggerOffset; // trigger offset (ns), incl latency & readout offset
  //UInt_t fNChambers;  // # chambers
  //UInt_t* fNROPlanes;   // # planes in each chamber
  UInt_t   fRNIon;    // number of ions
  std::vector<IonPar_t> fRIon;
  Double_t fRSMax;
  Double_t fRTotalCharge;
  Double_t fRTime0;
  //Double_t fTimeZero;
  
  std::vector<Double_t> fSumA;
  //std::vector<Short_t>  fDADC;

  Double_t fAvaGain;
  Int_t fNSamples;
  
  //zero suppression and common mode
  Bool_t fDoZeroSup;
  Double_t fZeroSup;
  Bool_t fDoCommonMode;
  std::vector<Double_t> fCommonModeArray;
  
  Short_t ADCConvert(Double_t val, Double_t off, Double_t gain, Int_t bits);
  Double_t PulseShape(Double_t t, 
		      Double_t C,  // normalization factor
		      Double_t Tp); // shaping time 
  //ClassDef (GEMSimDigitizer, 0) 
#if DBG_HISTOS>0
  TH2D* h2D_nplanesV_ava_dx; 
  TH2D* h2D_nplanesV_ava_dxs; 
  TH2D* h2D_nplanesV_ava_dy; 
  TH2D* h2D_nplanesV_ava_dys; 
  TH2D* h2D_nplanesV_ava_nstrips; 
  TH2D* h2D_nplanesV_ava_nx; 
  TH2D* h2D_nplanesVnActiveStrips;
  TH2D* h2D_nplanesVnAllHitStrips;
  TH2D* h2D_nplanesVnADCSum;
#endif
  /*
  TH2D* h2D_edepVdr;
  
  TH2D* h1_AvaSizeYvsX_SemiAna;
  TH2D* h1_AvaSizeYvsX_FastAppx;
  TH2D* h1_AvaSizeVsZion_SemiAna;
  TH2D* h1_AvaSizeVsTTime_SemiAna;
  
  TH1D* h1_SumweightsSemiAna;
  TH2D* h2D_SumweightsFastAppx;

  TH1D* h1_GammaEffSemiAna;
  TH2D* h2D_GammaEffFastAppx;
  
  TH1D* h1_NormSemiAna;
  TH1D* h1_NormFastAppx;

  TH1D* h1_SigmaEff;
  TH1D* h1_NionsPix;
 
  TH1D** h1_nbins_X;
  TH1D** h1_nbins_Y;

  TH1D** h1_binw_X;
  TH1D** h1_binw_Y;
  
  TH1D* h1_modhit_s;
  TH1D* h1_xhit_s;
  TH1D* h1_yhit_s;
  TH1D* h1_zhit_s;
  TH1D* h1_xdiff_s;
  TH1D* h1_ydiff_s;
  TH1D* h1_thit_s;
  TH1D* h1_edep_s;
  TH1D* h1_nions_s;
  
  TH1D* h1_modhit_b;
  TH1D* h1_xhit_b;
  TH1D* h1_yhit_b;
  TH1D* h1_zhit_b;
  TH1D* h1_xdiff_b;
  TH1D* h1_ydiff_b;
  TH1D* h1_thit_b;
  TH1D* h1_edep_b;
  TH1D* h1_nions_b;
  
  TH1D** h1_nstrips_X;
  TH1D** h1_nstrips_Y;

  TH2D** h1_ds_X;
  TH2D** h1_ds_Y;

  TH1D** h1_nbins_X;
  TH1D** h1_nbins_Y;

  TH1D** h1_fSumA_X;
  TH1D** h1_fSumA_Y;
  
  TH2D* h1_QvsX_ion;
  TH2D* h1_QvsY_ion;
  TH2D* h1_QnormvsX_ion;
  TH2D* h1_QnormvsY_ion;
  TH2D* h1_QareavsX_ion;
  TH2D* h1_QareavsY_ion;
  TH2D* h1_QintvsX_ion;
  TH2D* h1_QintvsY_ion;

  TH2D* h1_QvsX_ava;
  TH2D* h1_QvsY_ava;
  TH2D* h1_QintYvsX_ava;
  TH2D* h1_QintYvsY_ava;
  
  TH1D* h1_yGEM_preion;
  TH1D* h1_yGEM_preava;
  TH1D* h1_yGEM_inava;
  TH1D* h1_yGEM_inava_2;
  TH1D* h1_yGEM_inava_3;
  TH1D* h1_yGEM_inava_4;
  TH2D* h1_xGEMvsADC_inava_4;
  TH2D* h1_yGEMvsADC_inava_4;
  TH1D* h1_yGEM_incheckout;
  
  TH2D* h1_ava_yint;
  TH1D* h1_ava_int;
  */
  
  //std::chrono::steady_clock fClock_dbg;
  double fTotalTime_ion;
  double fTotalTime_ava;
  double fTotalTime_int;
  
  std::chrono::time_point<std::chrono::steady_clock> fStart;
  std::chrono::time_point<std::chrono::steady_clock> fEnd;
  std::chrono::duration<double> fDiff;
};

#endif



#include "GEMSimDigitizer.h"
#include "GEMDet.h"
//#include "gmn_tree.h"
#include "remoll_tree.h"

//gas parameters
#define fGasWion 26.         // eV
#define fGasDiffusion 1.e5       // mm2/s
#define fGasDriftVelocity 5.5e7   // mm/s
#define fAvalancheFiducialBand 10. // number of sigma defining the band around the avalanche in readout plane
#define fAvalancheChargeStatistics 1 // 0 Furry, 1 Gaussian
#define fGainMean 8.e3
#define fGain0 20.
#define fMaxNIon 1.e4               //maximum amount of ion pairs allowed in the digitization

#define fSNormNsigma 18.          //fSNormNsigma is an arbitrary multiplicative fact  
//#define fAvaGain 20.
#define fLateralUncertainty 0.

//electronics parameters
#define fAPVTimeJitter 25.    // time jitter associated with the APV internal clock
  
#define fEleSamplingPoints 6
#define fEleSamplingPeriod 25. // ns
#define fADCoffset 0.         // ADC offset
#define fADCgain 1.          // ADC gain
#define fADCbits 12         // ADC resolutions in bits
#define fGateWidth 400.    // to be changed , ns - pulse shape width at ~1/10 max

//parameter for GEM pedestal noise
#define fPulseNoiseSigma 20.   // additional sigma term of the pedestal noise
#define fPulseNoisePeriod 200. // period of the pedestal noise, assuming sinusoidal function
#define fPulseNoiseAmpConst 0  // constant term of the pedestal noise amplitude
#define fPulseNoiseAmpSigma 0  // sigma term of the pedestal noise amplitude

//paramters for cross-talk
#define fNCStripApart 32  // # of strips the induced signal is away from the mean signal
#define fCrossFactor 0.1  //reduction factor for the induced signal
#define fCrossSigma 0.03  //uncertainty of the reduction factor

// Pulse shaping parameters
#define fPulseShapeTau 56.   // [ns] GEM model 0 = 50. in SiD model
//#define fPulseShapeTau1 0.   // [ns] GEM model only; if negative assume SiD model

#define fEntranceRef -1.5  // z position of the copper layer right before the first GEM gas layer,             // relative to the center of the GEM chamber
#define fExitRef 1.5  // z position of the copper layer right before the first GEM gas layer,             // relative to the center of the GEM chamber
                         // which introduce additional uncertainty in the lateral direction
#define fRoutZ 9.185   // z-distance hit entrance to readout plane [mm]

//numerical integration parameters
#define fYIntegralStepsPerPitch 4
#define fXIntegralStepsPerPitch 4

#define fNROPlanes 2
#define fStripPitch 8.e-4

#include <iomanip>

using namespace std;

//stupid but: let's do the defautl constructor:
GEMSimDigitizer::GEMSimDigitizer()// :
  //fGasWion(), fGasDiffusion(), fGasDriftVelocity(), fAvalancheFiducialBand(), fAvalancheChargeStatistics(), fGainMean(), fGain0(), fMaxNIon(), fSNormNsigma(), fAvaGain(), fAPVTimeJitter() 
{
  fAvaGain = 20.0;
  fNSamples = 6;
}

GEMSimDigitizer::GEMSimDigitizer(int nchambers, double* trigoffset, double* gain, double zsup_thr, int napv, double* commonmode_array) : fZeroSup(zsup_thr) 
{
  //TODO: pass as parameters.
  fAvaGain = gain[0];
  fNSamples = 6;
  
  for(int i = 0; i<nchambers; i++){
    fTriggerOffset.push_back(trigoffset[i]);
    cout << i << "/" << nchambers << ": " << fTriggerOffset[i] << endl;
  }
  if(fZeroSup>0)fDoZeroSup = true;
  if(napv){
    fDoCommonMode = true;
    for(int i = 0; i<napv; i++){
      fCommonModeArray.push_back(commonmode_array[i]);
      cout << i << "/" << napv << ": " << fCommonModeArray[i] << endl;
    }
  }
  fRIon.resize((int)fMaxNIon);

#if DBG_HISTOS > 0
  h2D_nplanesV_ava_dx = new TH2D("h2D_nplanesV_ava_dx", "FT;AVA_dx;Module", 100, 0, 20, nchambers*2, 0, nchambers*2);
  h2D_nplanesV_ava_dxs = new TH2D("h2D_nplanesV_ava_dxs", "FT;AVA_dxs;Module", 100, 0, 20, nchambers*2, 0, nchambers*2);
  h2D_nplanesV_ava_dy = new TH2D("h2D_nplanesV_ava_dy", "FT;AVA_dy;Module", 100, 0, 20, nchambers*2, 0, nchambers*2);
  h2D_nplanesV_ava_dys = new TH2D("h2D_nplanesV_ava_dys", "FT;AVA_dys;Module", 100, 0, 20, nchambers*2, 0, nchambers*2);
  h2D_nplanesV_ava_nstrips = new TH2D("h2D_nplanesV_ava_nstrips", "FT;AVA_nstrips;Module", 100, 0, 100, nchambers*2, 0, nchambers*2);
  h2D_nplanesV_ava_nx = new TH2D("h2D_nplanesV_ava_nx", "FT;AVA_nx;Module", 100, 0, 100, nchambers*2, 0, nchambers*2);
  h2D_nplanesVnActiveStrips = new TH2D("h2D_nplanesVnActiveStrips", "FT;Strips above threshold;Module", 50, 0, 50, nchambers*2, 0, nchambers*2);
  h2D_nplanesVnAllHitStrips = new TH2D("h2D_nplanesVnAllHitStrips", "FT;All Hit Strips;Module", 100, 0, 100, nchambers*2, 0, nchambers*2);
  h2D_nplanesVnADCSum = new TH2D("h2D_nplanesVnADCSum", "FT;ADC sum;Module", 250, 0, 10000, nchambers*2, 0, nchambers*2);
#endif
  /*
  h2D_edepVdr = new TH2D("h2D_edepVdr", ";sqrt(dx2_hit+dy2_hit);edep;", 100, 0., 50., 100, 0., 1.e5);

  h1_AvaSizeYvsX_SemiAna = new TH2D("h1_AvaSizeYvsX_SemiAna", "", 100, 0.0, 100, 100, 0.0, 100);
  h1_AvaSizeYvsX_FastAppx = new TH2D("h1_AvaSizeYvsX_FastAppx", "", 100, 0.0, 100, 100, 0.0, 100);
  h1_AvaSizeVsZion_SemiAna = new TH2D("h1_AvaSizeVsZion_SemiAna", "", 100, -2.0, 2.0, 100, 0.1, 0.2);
  h1_AvaSizeVsTTime_SemiAna = new TH2D("h1_AvaSizeVsTTime_SemiAna", "", 150, 0., 300, 100, 0.1, 0.2);
  
  h1_SumweightsSemiAna = new TH1D("h1_SumweightsSemiAna", "", 100, 0., 100);
  h2D_SumweightsFastAppx = new TH2D("h2D_SumweightsFastAppx", "", 100, 0., 50., 750, 0., 750);
  
  h1_GammaEffSemiAna = new TH1D("h1_GammaEffSemiAna", "", 100, 0.1, 0.2);
  h2D_GammaEffFastAppx = new TH2D("h2D_GammaEffFastAppx", "", 100, 0., 50., 100, 0.0, 10.0);
  
  h1_NormSemiAna = new TH1D("h1_NormSemiAna", "", 1000, 0, 1.e5);
  h1_NormFastAppx = new TH1D("h1_NormFastAppx", "", 1000, 0, 1.e5);

  h1_SigmaEff = new TH1D("h1_SigmaEff", "", 150, 0., 0.30);
  h1_NionsPix = new TH1D("h1_NionsPix", "", 100, 0., 1000.0);
  h1_nbins_X = new TH1D*[2];
  h1_nbins_Y = new TH1D*[2];
  h1_binw_X = new TH1D*[2];
  h1_binw_Y = new TH1D*[2];
  for(int i = 0; i<2; i++){
    h1_nbins_X[i] = new TH1D(Form("h1_nbins_X_%d", i), "", 100, 0, 500);
    h1_nbins_Y[i] = new TH1D(Form("h1_nbins_Y_%d", i), "", 100, 0, 500);
    h1_binw_X[i] = new TH1D(Form("h1_binw_X_%d", i), "", 100, 0, 0.2);
    h1_binw_Y[i] = new TH1D(Form("h1_binw_Y_%d", i), "", 100, 0, 0.2);
  }
  
  h1_ava_yint = new TH2D("h1_ava_yint", "", 100, -5.0, 5.0, 100, -1., 1.);
  h1_ava_int = new TH1D("h1_ava_int", "", 100, -1., 1.);
  
  h1_modhit_s = new TH1D("h1_modhit_s", "", 36, 0, 36);
  h1_xhit_s = new TH1D("h1_xhit_s", "", 205, -1.025, 1.025);
  h1_yhit_s = new TH1D("h1_yhit_s", "", 62, -0.31, 0.31);
  h1_zhit_s = new TH1D("h1_zhit_s", "", 100, -0.01, 0.01);
  h1_xdiff_s = new TH1D("h1_xdiff_s", "", 100, -0.1, 0.1);
  h1_ydiff_s = new TH1D("h1_ydiff_s", "", 100, -0.1, 0.1);
  h1_thit_s = new TH1D("h1_thit_s", "", 1000, -500, 500);
  h1_edep_s = new TH1D("h1_edep_s", "", 1000, 0, 1.e6);
  h1_nions_s = new TH1D("h1_nions_s", "", 100, 0, 1000);
  
  h1_modhit_b = new TH1D("h1_modhit_b", "", 36, 0, 36);
  h1_xhit_b = new TH1D("h1_xhit_b", "", 205, -1.025, 1.025);
  h1_yhit_b = new TH1D("h1_yhit_b", "", 62, -0.31, 0.31);
  h1_zhit_b = new TH1D("h1_zhit_b", "", 100, -0.01, 0.01);
  h1_xdiff_b = new TH1D("h1_xdiff_b", "", 100, -0.1, 0.1);
  h1_ydiff_b = new TH1D("h1_ydiff_b", "", 100, -0.1, 0.1);
  h1_thit_b = new TH1D("h1_thit_b", "", 1000, -500, 500);
  h1_edep_b = new TH1D("h1_edep_b", "", 1000, 0, 1.e6);
  h1_nions_b = new TH1D("h1_nions_b", "", 100, 0, 1000);
  
  h1_nstrips_X = new TH1D*[2];
  h1_nstrips_Y = new TH1D*[2];
  
  h1_ds_X = new TH2D*[2];
  h1_ds_Y = new TH2D*[2];
  
  h1_nbins_X = new TH1D*[2];
  h1_nbins_Y = new TH1D*[2];
  
  h1_fSumA_X = new TH1D*[2];
  h1_fSumA_Y = new TH1D*[2];

  for(int i = 0; i<2; i++){
    h1_nstrips_X[i] = new TH1D(Form("h1_nstrips_X_%d", i), "", 100, 0, 100);
    h1_nstrips_Y[i] = new TH1D(Form("h1_nstrips_Y_%d", i), "", 100, 0, 100);
    
    h1_ds_X[i] = new TH2D(Form("h1_ds_X_%d", i), "", 100, -5.e1, 5.e1, 100, -5.e1, 5.e1);
    h1_ds_Y[i] = new TH2D(Form("h1_ds_Y_%d", i), "", 100, -5.e1, 5.e1, 100, -5.e1, 5.e1);
    
    h1_nbins_X[i] = new TH1D(Form("h1_nbins_X_%d", i), "", 100, 0, 10000);
    h1_nbins_Y[i] = new TH1D(Form("h1_nbins_Y_%d", i), "", 100, 0, 10000);
    
    h1_fSumA_X[i] = new TH1D(Form("h1_fSumA_X_%d", i), "", 100, 0, 1.0e6);
    h1_fSumA_Y[i] = new TH1D(Form("h1_fSumA_Y_%d", i), "", 100, 0, 1.0e6);
  }
  
  h1_QvsX_ion = new TH2D("h1_QvsX_ion", "", 250, -0.25, 0.25, 200, 0, 2.e4);
  h1_QvsY_ion = new TH2D("h1_QvsY_ion", "", 200, -0.2, 0.2, 200, 0, 2.e4);
  h1_QnormvsX_ion = new TH2D("h1_QnormvsX_ion", "", 250, -0.25, 0.25, 200, 0, 2.e4);
  h1_QnormvsY_ion = new TH2D("h1_QnormvsY_ion", "", 200, -0.2, 0.2, 200, 0, 2.e4);
  h1_QareavsX_ion = new TH2D("h1_QareavsX_ion", "", 750, -0.75, 0.75, 100, 0, 100.);
  h1_QareavsY_ion = new TH2D("h1_QareavsY_ion", "", 200, -0.2, 0.2, 100, 0, 100.);
  h1_QintvsX_ion = new TH2D("h1_QintvsX_ion", "", 250, -0.25, 0.25, 1000, 0, 1.e6);
  h1_QintvsY_ion = new TH2D("h1_QintvsY_ion", "", 200, -0.2, 0.2, 1000, 0, 1.e6);

  h1_QvsX_ava = new TH2D("h1_QvsX_ava", "", 750, -0.75, 0.75, 1000, 0, 1.e5);
  h1_QvsY_ava = new TH2D("h1_QvsY_ava", "", 200, -0.2, 0.2, 1000, 0, 1.e5);
  h1_QintYvsX_ava = new TH2D("h1_QintYvsX_ava", "", 750, -0.75, 0.75, 1000, 0, 1.e5);
  h1_QintYvsY_ava = new TH2D("h1_QintYvsY_ava", "", 200, -0.2, 0.2, 1000, 0, 1.e5);
  
  h1_yGEM_preion = new TH1D("h1_yGEM_preion", "", 200, -0.2, 0.2);
  h1_yGEM_preava = new TH1D("h1_yGEM_preava", "", 200, -0.2, 0.2);
  h1_yGEM_inava = new TH1D("h1_yGEM_inava", "", 200, -0.2, 0.2);
  h1_yGEM_inava_2 = new TH1D("h1_yGEM_inava_2", "", 200, -0.2, 0.2);
  h1_yGEM_inava_3 = new TH1D("h1_yGEM_inava_3", "", 200, -0.2, 0.2);
  h1_yGEM_inava_4 = new TH1D("h1_yGEM_inava_4", "", 200, -0.2, 0.2);
  h1_xGEMvsADC_inava_4 = new TH2D("h1_xGEMvsADC_inava_4", "", 750, -0.75, 0.75, 1000, -19, 20000-19);
  h1_yGEMvsADC_inava_4 = new TH2D("h1_yGEMvsADC_inava_4", "", 200, -0.2, 0.2, 1000, -19, 20000-19);
  h1_yGEM_incheckout = new TH1D("h1_yGEM_incheckout", "", 200, -0.2, 0.2);
  */
  
  fTotalTime_ion = 0;
  fTotalTime_ava = 0;
  fTotalTime_int = 0;
}


GEMSimDigitizer::~GEMSimDigitizer()
{
}


//.......................................................
// ionization Model
//
void
GEMSimDigitizer::IonModel(TRandom3* R,
			  const TVector3& xi,
			  const TVector3& xo,
			  const Double_t elost) // eV
{
#define DBG_ION 0

  TVector3 vseg = xo-xi; // mm
  
  // ---- extract primary ions from Poisson
  fRNIon = R->Poisson(elost/fGasWion);

  if (fRNIon <=0)
    return;

#if DBG_ION > 0
  cout << "E lost = " << elost << ", " << fRNIon << " ions" << endl;
#endif
  if (fRNIon > fMaxNIon) {
#if DBG_ION > 0
    cout << __FUNCTION__ << ": WARNING: too many primary ions " << fRNIon << " limit to "
	 << fMaxNIon << endl;
#endif
    fRNIon = fMaxNIon;
  }

  fRSMax = 0.;
  fRTotalCharge = 0;
  fRTime0 = 999999.; // minimum time of drift

  for (UInt_t i=0;i<fRNIon;i++) { // first loop used to evaluate quantities
    IonPar_t ip;

    Double_t lion = R->Uniform(0.,1.); // position of the hit along the track segment (fraction)

    //In principle, the lateral uncertainty should have been put in the Ava model, but not here
    //But since we are not simulating the details of the avalanche, I think it is ok (Weizhi)
    ip.X = vseg.X()*lion+xi.X() + R->Gaus(0., fLateralUncertainty);
    ip.Y = vseg.Y()*lion+xi.Y() + R->Gaus(0., fLateralUncertainty);

    // Note the definition of fRoutZ is the distance from xi.Z() to xrout.Z():
    //        xi               xo   xrout
    // |<-LD->|<-----vseg----->|    |
    // |<-------fRoutZ---------|--->|
    // |      |<-lion*vseg->   |    |
    // |      |             <--LL-->|
    
    Double_t LD = TMath::Abs(xi.Z() - fEntranceRef);//usually should be 0,
                                            //unless particle is produced inside the gas layer

    Double_t LL = TMath::Abs(fRoutZ - LD - vseg.Z()*lion);
    Double_t ttime = LL/fGasDriftVelocity; // traveling time from the drift gap to the readout
    
    //cout << " rout Z  (mm?) " << fRoutZ << ", LD (mm?) " << LD << " vseg Z (mm?) " << vseg.Z()  << endl;
    //cout << " travelling length (mm?) " << LL << ", travelling time:  " <<  ttime << endl;
    
    fRTime0 = TMath::Min(ttime, fRTime0); // minimum traveling time [s]

    ip.SNorm = TMath::Sqrt(2.*fGasDiffusion*ttime); // spatial sigma on readout [mm] - width spread due to diffusion
    // cout<<"ip.SNorm: "<<ip.SNorm<<endl;
    //  cout<<"vseg: "<<vseg.X()<<" deltaX: "<<vseg.X()*lion<<endl;
    if( fAvalancheChargeStatistics == 1 ) {
      Double_t gnorm = fGainMean/TMath::Sqrt(fGain0); // overall gain TBC
      ip.Charge = R->Gaus(fGainMean, gnorm); // Gaussian distribution of the charge
    }
    else {
      ip.Charge = R->Exp(fGainMean); // Furry distribution
    }

    if( ip.Charge > 0 )
      fRTotalCharge += ip.Charge;
    else
      ip.Charge = 0;

    //h1_AvaSizeVsZion_SemiAna->Fill(xi.Z()+vseg.Z()*lion, ip.SNorm);
    //h1_AvaSizeVsTTime_SemiAna->Fill(ttime*1.e9, ip.SNorm);
    
    fRSMax = TMath::Max(ip.SNorm, fRSMax);

    // Derived quantities needed by the numerical integration in AvaModel
    ip.SNorm *= fSNormNsigma;
    ip.R2 = ip.SNorm * ip.SNorm;
    ip.ggnorm = ip.Charge * TMath::InvPi() / ip.R2; // normalized charge

#if DBG_ION > 1
    printf("z coords %f %f %f %f lion %f LL %lf\n",
	   xi.Z(), xo.Z(), vseg.Z(), lion, LL);
    printf("ttime = %e\n", ttime);
#endif
#if DBG_ION > 0
    cout << " x, y = " << ip.X << ", " << ip.Y << " snorm = "
	 << ip.SNorm/fSNormNsigma << " charge " << ip.Charge << endl;
    cout << "fRTime0 = " << fRTime0 << endl;
    cout << "fRion size " << fRIon.size() << " " << i << endl;
#endif
    /*
    h1_QvsX_ion->Fill(ip.X, ip.Charge);
    h1_QvsY_ion->Fill(ip.Y, ip.Charge);
    h1_QvsX_ion->Fill(ip.X, ip.ggnorm);
    h1_QvsY_ion->Fill(ip.Y, ip.ggnorm);
    */
    fRIon[i] = ip;
  }
  return;
}

Short_t 
GEMSimDigitizer::ADCConvert(Double_t val, Double_t off, Double_t gain, Int_t bits)
{
  // Convert analog value 'val' to integer ADC reading with 'bits' resolution
  assert( bits >= 0 && bits <= 12 );

  if( val < 0. )
    val = 0.;
  Double_t vvv = (val - off)/gain;
  //std::cout<<val<<" : "<<vvv<<std::endl;
  //printf("offset = %1.3f, gain = %1.3f, input value = %1.3f , output value = %1.3f  \n", off, gain, val, vvv);
  Double_t saturation = static_cast<Double_t>( (1<<bits)-1 );
  if( vvv > saturation )
    vvv = saturation;

  Short_t dval =
    static_cast<Short_t>( TMath::Floor( (vvv>saturation) ? saturation : vvv ));

  //  cerr << val << " dval = " << dval << endl;
  if( dval < 0 ) dval = 0;
  return dval;

}

// Pulse Shape SiD model
// APV25 time function from  M. Friedl et al NIMA 572 (2007) pg 385-387 (APV25 for silicon!)
//

Double_t 
GEMSimDigitizer::PulseShape(Double_t t, 
			    Double_t C,  // normalization factor
			    Double_t Tp) // shaping time 
{

  Double_t v;
  Double_t x;
  x = t/Tp;
  v = C/Tp * x * TMath::Exp(-x);
  
  return ( v>0. ) ? v : 0.;

}


#define DBG_AVA 0
//.......................................................
// avalanche model
//
// integration methods
//
void GEMSimDigitizer::Integration_semiana(double roangle, 
					   double xl, double xr, 
					   double yb, double yt, 
					   int nx, double xbw)
{
  double amplitude_sum = 0;
#if DBG_AVA > 1
  cout << "Integration_semi_ana: roangle (deg)" << roangle*TMath::RadToDeg() << endl;
  cout << "x limits: l " << xl << " r " << xr << " y limits: b " << yb << " " << yt << endl;
#endif
    
  for (UInt_t i = 0; i < fRNIon; i++){
	//This is just a rotation from the Cartesian to the readout direction
    Double_t frxs = fRIon[i].X*cos(roangle) - fRIon[i].Y*sin(roangle);
    Double_t frys = fRIon[i].X*sin(roangle) + fRIon[i].Y*cos(roangle);
    
    Int_t ix = (frxs-xl) / xbw;// the bin number
    Int_t dx = fRIon[i].SNorm / xbw  + 1; // # of bins the ion cloud contains.
#if DBG_AVA > 2
    cout << "frxs " << frxs << " frys " << frys << endl;
    cout << "ix dx " << ix << " " << dx << endl;
#endif
    
    //
    // NL change:
    //
    // ggnorm is the avalance charge for the i^th ion, and R2 is the square of the radius of the diffusion 
    // circle, mutiplied by the kSNormNsigma factor: (ip.SNorm * ip.SNorm)*kSNormNsigma*kSNormNsigma. All 
    // strips falling within this circle are considered in charge summing. 
    //
    // The charge contribution to a given strip by the i^th ion is evaluated by a Lorentzian (or Gaussian)
    // distribution; the sigma for this distribution is eff_sigma, which is the actual avalance sigma. 
    //
    Double_t ggnorm = fRIon[i].ggnorm;
    Double_t r2 = fRIon[i].R2;
    Double_t eff_sigma_square = r2/(fSNormNsigma*fSNormNsigma);
    Double_t eff_sigma = TMath::Sqrt(eff_sigma_square);
    Double_t current_ion_amplitude = fAvaGain*ggnorm*(1./(TMath::Pi()*eff_sigma))*(eff_sigma*eff_sigma);
   // (1./(TMath::Pi()*eff_sigma))comes from a normalization factor that has its origin in a 2D (or effectively 2D) distribution and/or the analytic form of the ion-induced current shape.
   //(eff_sigma*eff_sigma) comes from integrating over the transverse extent of the ion cloud (or equivalently, from the fact that the total charge scales with the area  of the cloud).
    amplitude_sum+= current_ion_amplitude;
    
    double inte4 = 0.;
    double ceff, yinte;
    
    // xc and yc are center of current bin
    //double sumA = 0;
    // Loop over bins
    Int_t min_binNb_x = max(ix-dx,0);
    Int_t max_binNb_x = min(ix+dx+1,nx);
    Int_t jx = min_binNb_x;
    Double_t xc = xl + (jx+0.5) * xbw;
    for (; jx < max_binNb_x; ++jx, xc += xbw){
      Double_t xd2 = frxs-xc; xd2 *= xd2;
      if( xd2 > r2 ){
	if( (xc - frxs)>0 )
	  break;
	else
	  continue;
      }
      
      ceff = sqrt(xd2 + eff_sigma_square); //It is sum of eff_sigma_square and the square of charge center and bin center
      yinte = 1./ceff *(atan((yt-frys)/ceff) - atan((yb-frys)/ceff));//Analytic integral over Y - holding x-fixed
      //These lines compute the analytic integral of a 2D Gaussian ion cloud over the rectangular area of a strip, using the arctangent formula
      inte4+= yinte*xbw;
      fSumA[jx] += (yinte*xbw)*current_ion_amplitude;//This gives the total charge fraction landing on strip jx.
    }
    //h1_SumweightsSemiAna->Fill(inte4);
    //h1_GammaEffSemiAna->Fill(eff_sigma);
  }
// #if DBG_AVA >1
// #endif
  //h1_NormSemiAna->Fill(amplitude_sum);
}

void GEMSimDigitizer::Integration_fastappx(TRandom3* R, double roangle, 
					   double nions_strip,
					   double xc_hit, double yc_hit,
					   double dx2_hit, double dy2_hit,
					   double xl, double xr, 
					   double yb, double yt, 
					   int nx, double xbw, 
					   int ny, double ybw)
{
#if DBG_AVA > 0
  cout << " integration_fastappx " << fRNIon << " " << fRSMax << " " << fRTime0 << endl;
#endif
  //so we do the numerical stuff, but only on the average hit... and then what?
  Double_t frxs = xc_hit*cos(roangle) - yc_hit*sin(roangle);
  Double_t frys = xc_hit*sin(roangle) + yc_hit*cos(roangle);
  //for (UInt_t i = 0; i < fRNIon; i++){
  //
  // bin containing center and # bins each side to process
  Int_t ix = (frxs-xl) / xbw;
  Int_t iy = (frys-yb) / ybw;
  //Int_t dx = fRIon[i].SNorm / xbw  + 1;
  //Int_t dy = fRIon[i].SNorm / ybw  + 1;
  Double_t LL = TMath::Abs(fRoutZ - R->Uniform(fEntranceRef,fExitRef));
  Double_t ttime = LL/fGasDriftVelocity;
  Double_t SNorm = TMath::Sqrt(2.*fGasDiffusion*ttime)*fSNormNsigma;
  Int_t dx =  SNorm / xbw  + 1;
  Int_t dy =  SNorm / ybw  + 1;
  
#if DBG_AVA > 1
  cout << "ix dx iy dy " << ix << " " << dx << " " << iy << " " << dy << endl;
#endif
  
  //double NionsStrip = fRNIon/sqrt(dx2_hit);
  //double Ld_ion = fRNIon/sqrt(dx2_hit+dy2_hit);
  //
  // NL change:
  //
  // ggnorm is the avalance charge for the i^th ion, and R2 is the square of the radius of the diffusion 
  // circle, mutiplied by the kSNormNsigma factor: (ip.SNorm * ip.SNorm)*kSNormNsigma*kSNormNsigma. All 
  // strips falling within this circle are considered in charge summing. 
  //
  // The charge contribution to a given strip by the i^th ion is evaluated by a Lorentzian (or Gaussian)
  // distribution; the sigma for this distribution is eff_sigma, which is the actual avalance sigma. 
  //
  Double_t Charge = R->Gaus(fGainMean, fGainMean/TMath::Sqrt(fGain0));
  Double_t r2 = SNorm*SNorm+dx2_hit+dy2_hit;// ?
  Double_t ggnorm = Charge * TMath::InvPi() / r2;
  Double_t eff_sigma_square = r2/(fSNormNsigma*fSNormNsigma);
  Double_t eff_sigma = TMath::Sqrt(eff_sigma_square);
  Double_t Amplitude = fAvaGain*fRNIon*ggnorm*(1./(TMath::Pi()*eff_sigma))*(eff_sigma*eff_sigma);

  //double NionsStrip = Ld_ion*eff_sigma*R->Gaus(4., 1.);//
  
  //h1_NormFastAppx->Fill(Amplitude);
  //h1_SigmaEff->Fill(eff_sigma);
  //h1_NionsPix->Fill(NionsPix);
  
  // xc and yc are center of current bin
  // Loop over bins
  Int_t min_binNb_x = max(ix-dx,0);
  Int_t min_binNb_y = max(iy-dy,0);
  Int_t max_binNb_x = min(ix+dx+1,nx);
  Int_t max_binNb_y = min(iy+dy+1,ny);
  Int_t jx = min_binNb_x;
  Double_t xc = xl + (jx+0.5) * xbw;
  double sumweights = 0;
  //double sumweights_reg = 0;
  double weight;
  for (; jx < max_binNb_x; ++jx, xc += xbw){
    Double_t xd2 = frxs-xc; xd2 *= xd2;
    
    Int_t jx_base = jx * ny;
    Int_t jy = min_binNb_y;
    Double_t yc = yb + (jy+0.5) * ybw;
    
    //double b_smear = R->Poisson(nions_strip)/nions_strip;//should fluctuate around 1.

    for (; jy < max_binNb_y; ++jy, yc += ybw){
      Double_t yd2 = frys-yc; yd2 *= yd2;
      
      if( xd2 + yd2 <= r2 ) {
	weight = 1./ // b_smear / 
	  (eff_sigma_square + xd2/(1+dx2_hit) +yd2/(1+dy2_hit) );
	//(eff_sigma_square + xd2 +yd2 );
	fSumA[jx_base+jy] += weight*Amplitude;
	sumweights+= weight;
	//sumweights_reg+= 1. / (eff_sigma_square + xd2 +yd2 );
	//current_ion_amplitude / ((xd2+yd2)+eff_sigma_square);
      }
    }
  }
  /**/
  //cout << sumweights << " " << NionsPix << endl;
  //h2D_SumweightsFastAppx->Fill(sqrt(dx2_hit+dy2_hit), sumweights*xbw*ybw);
  //h2D_GammaEffFastAppx->Fill(sqrt(dx2_hit+dy2_hit), sqrt(eff_sigma_square*(1+dx2_hit)*(1+dy2_hit)));
  //reloop to normalize the individual weights.
  /*
  jx = min_binNb_x;
  for (; jx < max_binNb_x; ++jx){
    Int_t jx_base = jx * ny;
    Int_t jy = min_binNb_y;
    for (; jy < max_binNb_y; ++jy){
      //fSumA[jx_base+jy]*=fAvaGain*Charge;
      fSumA[jx_base+jy]*= sumweights_reg/sumweights;
    }
  }
  */
  //}
  /*
  //Loop on strips instead of bins...
  double strippitch_mm = fStripPitch*1.e3;
  int nstrips_x = (xr-xl)/strippitch_mm;
  int nstrips_y = (yt-yb)/strippitch_mm;
  fSumA.resize(nstrips_x*nstrips_y);
  Double_t xc = xl + fStripPitch/2.0;
  double sumweights = 0;
  cout << nstrips_x << " " << nstrips_y << endl;
  double weight;
  for(int ix = 0; ix<nstrips_x; xc+= strippitch_mm, ix++){
    Double_t yc = yb + fStripPitch/2.0;
    Double_t xd2 = frxs-xc; xd2 *= xd2;
    
    Int_t jx_base = ix * nstrips_y;
    
    for(int jy = 0; jy<nstrips_y; yc+= strippitch_mm, jy++){
      Double_t yd2 = frys-yc; yd2 *= yd2;
      
      if( xd2 + yd2 <= r2 ) {
	weight = 1. / // R->Poisson(NionsPix) / 
	  (eff_sigma_square + xd2/(1+dx2_hit) +yd2/(1+dy2_hit) );
	  //(eff_sigma_square + xd2 +yd2 );
	fSumA[jx_base+jy] += weight*Amplitude;
	sumweights+= weight;
	//current_ion_amplitude / ((xd2+yd2)+eff_sigma_square);
      }
      
    }
  }
  h1_Sumweights->Fill(dx2_hit+dy2_hit, sumweights*fStripPitch*fStripPitch);
  */
  
}



//TGEMSBSGEMHit **
void GEMSimDigitizer::AvaModel(const int ic,
			       GEMDet* gemdet, 
			       TRandom3* R,
			       const TVector3& xi,
			       const TVector3& xo,
			       const Double_t t0)
{
#if DBG_AVA > 0
  cout << "Chamber " << ic << "----------------------------------" << endl;
  cout << "In  " << xi.X() << " " << xi.Y() << " " << xi.Z() << endl;
  cout << "Out " << xo.X() << " " << xo.Y() << " " << xo.Z() << endl;
#endif

  // xi, xo are in chamber frame, in mm

  Double_t nsigma = fAvalancheFiducialBand; // coverage factor
  
  //trying something:
  int integral_steps_x = fXIntegralStepsPerPitch;
  int integral_steps_y = fYIntegralStepsPerPitch;
  
#if DBG_AVA > 0
  cout << "fRSMax, nsigma " << fRSMax << " " << nsigma << endl;
#endif
  
  Double_t x0,y0,x1,y1; // lower and upper corners of avalanche diffusion area
  
  if (xi.X()<xo.X()) {
    x0 = xi.X()-nsigma*fRSMax;
    x1 = xo.X()+nsigma*fRSMax;
  } else {
    x1 = xi.X()+nsigma*fRSMax;
    x0 = xo.X()-nsigma*fRSMax;
  }

  if (xi.Y()< xo.Y()) {
    y0 = xi.Y()-nsigma*fRSMax;
    y1 = xo.Y()+nsigma*fRSMax;
  } else {
    y1 = xi.Y()+nsigma*fRSMax;
    y0 = xo.Y()-nsigma*fRSMax;
  }
  
  //h1_AvaSizeYvsX_SemiAna->Fill(x1-x0, y1-y0);
  
  // Check if any part of the avalanche region is in the active area of the sector.
  // Here, "active area" means the chamber's *bounding box*, which is
  // larger than the wedge's active area (section of a ring)
  
  //const TGEMSBSGEMChamber& chamber = spect.GetChamber(ic);
  Double_t glx = (-gemdet->GEMPlanes[ic*2].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;//(chamber.GetPlane(0).GetStripLowerEdge(0)+chamber.GetPlane(0).GetSPitch()/2.0) * 1000.0;
  Double_t gly = (-gemdet->GEMPlanes[ic*2+1].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;
  Double_t gux = (gemdet->GEMPlanes[ic*2].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;//(chamber.GetPlane(0).GetStripUpperEdge(chamber.GetPlane(0).GetNStrips()-1) -chamber.GetPlane(0).GetSPitch()/2.0) * 1000.0;
  Double_t guy = (gemdet->GEMPlanes[ic*2+1].dX())/2.*1.e3;//+fStripPitch)/2.0*1000.;//(chamber.GetPlane(1).GetStripUpperEdge(chamber.GetPlane(1).GetNStrips()-1) -chamber.GetPlane(1).GetSPitch()/2.0) * 1000.0;
 /* 
  if (x1<glx || x0>gux ||
      y1<gly || y0>guy) { // out of the sector's bounding box
    cerr << __FILE__ << " " << __FUNCTION__ << ": out of sector, "
	 << "chamber " << ic << " sector " << ic/30 << " plane " << ic%30 << endl
	 << "Following relations should hold:" << endl
	 << "(x1 " << x1 << ">glx " << glx << ") (x0 " << x0 << "<gux " << gux << ")" << endl
	 << "(y1 " << y1 << ">gly " << gly << ") (y0 " << y0 << "<guy " << guy << ")" << endl;
    //return 0;
  }
  */
  //bool bb_clipped = (x0<glx||y0<gly||x1>gux||y1>guy);
  if(x0<glx) x0=glx;
  if(y0<gly) y0=gly;
  if(x1>gux) x1=gux;
  if(y1>guy) y1=guy;

  // Loop over chamber planes
  double roangle_mod, dx_mod, xoffset_mod;
  int GEMstrips;
  
  double xt_factor;
  int isLeft;
  //TGEMSBSGEMHit **virs;
  //virs = new TGEMSBSGEMHit *[fNROPlanes[ic]];
  for (UInt_t ipl = 0; ipl < fNROPlanes; ++ipl){
#if DBG_AVA > 1
    cout << "coordinate " << ipl << " =========================" << endl;
#endif
     
    xt_factor = R->Gaus(fCrossFactor, fCrossSigma);
    isLeft = R->Uniform(1.) < 0.5 ? -1 : 1;
    
    //cout << xt_factor << endl;
    
    // Compute strips affected by the avalanche
    //const GEMPlane& pl = gemdet->GEMPlanes[ic*2+ipl];
    // strip angle = - strip angle for plane -> strip rotation purposes.
    roangle_mod = -gemdet->GEMPlanes[ic*2+ipl].ROangle();
    dx_mod = gemdet->GEMPlanes[ic*2+ipl].dX();
    xoffset_mod = gemdet->GEMPlanes[ic*2+ipl].Xoffset();
    GEMstrips = gemdet->GEMPlanes[ic*2+ipl].GetNStrips();
    
    // Positions in strip frame
    Double_t xsm = ( (x0+x1)*cos(roangle_mod) - (y0+y1)*sin(roangle_mod) ) * 0.5;//x of the charge center in the rotated strip coordinate system.
    Double_t ysm = ( (x0+x1)*sin(roangle_mod) + (y0+y1)*cos(roangle_mod) ) * 0.5;
    Double_t xs0 = xsm-(x1-x0)*0.5;//Lower x-boundary in the rotated strip coordinate system
    Double_t ys0 = ysm-(y1-y0)*0.5;
    Double_t xs1 = xsm+(x1-x0)*0.5;
    Double_t ys1 = ysm+(y1-y0)*0.5;
    // Double_t xs0 = x0*cos(roangle_mod) - y0*sin(roangle_mod);
    // Double_t ys0 = x0*sin(roangle_mod) + y0*cos(roangle_mod);
    // Double_t xs1 = x1*cos(roangle_mod) - y1*sin(roangle_mod); 
    // Double_t ys1 = x1*sin(roangle_mod) + y1*cos(roangle_mod);
#if DBG_AVA > 1
    cout << "glx gly gux guy " << glx << " " << gly << " " << gux << " " << guy << endl;
    cout << ic << " " << ipl << " " << roangle_mod 
	 << " xs0 ys0 xs1 ys1 " << xs0 << " " << ys0 << " " << xs1 << " " << ys1 << endl;
    cout << " ref: xs0 ys0 xs1 ys1 " << x0*cos(roangle_mod) - y0*sin(roangle_mod) << " " << x0*sin(roangle_mod) + y0*cos(roangle_mod)
	 << " " << x1*cos(roangle_mod) - y1*sin(roangle_mod) << " " << x1*sin(roangle_mod) + y1*cos(roangle_mod) << endl;
#endif
    //if(ipl==1 && ic<12)h1_yGEM_inava->Fill(xs0*1.e-3);
    //if(ipl==0 && ic<4){
    //h1_ds_X[int(ic>0)]->Fill(xs0-xs1, ys0-ys1);
    //}
    //if(ipl==1 && ic<4){
    //h1_ds_Y[int(ic>0)]->Fill(xs0-xs1, ys0-ys1);
    //}
#if DBG_HISTOS > 0
    h2D_nplanesV_ava_dx->Fill(x1-x0, min(ic,3)*2+ipl);
    h2D_nplanesV_ava_dxs->Fill(abs(xs1-xs0), min(ic,3)*2+ipl);
    h2D_nplanesV_ava_dy->Fill(y1-y0, min(ic,3)*2+ipl);
    h2D_nplanesV_ava_dys->Fill(abs(ys1-ys0), min(ic,3)*2+ipl);
#endif
    Int_t iL = max(0, Int_t((xs0*1.e-3+dx_mod/2.)/fStripPitch) );
    iL = min(iL, GEMstrips);
    //pl.GetStrip (xs0 * 1e-3, ys0 * 1e-3);
    Int_t iU = min(Int_t((xs1*1.e-3+dx_mod/2.)/fStripPitch), GEMstrips);
    iU = max(0, iU);
    //pl.GetStrip (xs1 * 1e-3, ys1 * 1e-3);
     
    // Check for (part of) the avalanche area being outside of the strip region
    //if( (iL <= 0 && iU <= 0) || (iL>=pl.GetNStrips() && iU>=pl.GetNStrips()) ) {
      // All of the avalanche outside -> nothing to do
      // TODO: what if this happens for only one strip coordinate (ipl)?
#if DBG_AVA > 1
      // cerr << __FILE__ << " " << __FUNCTION__ << ": out of active area, "
      // 	   << "chamber " << ic << " sector " << ic%30 << " plane " << ic/30 << endl
      // 	   << "iL_raw " << pl.GetStripUnchecked(xs0*1e-3) << " "
      // 	   << "iU_raw " << pl.GetStripUnchecked(xs1*1e-3) << endl
      // 	   << endl << endl;
    cout << "Low strip " << iL << " Up strip " << iU << " N strips? " << abs(iU-iL) << endl;
#endif
    if(iL==iU){//nothing to do
      return;
    }

    if(iU<iL)swap(iU, iL);


    //
    // Bounds of rectangular avalanche region, in strip frame
    //

    // Limits in x are low edge of first strip to high edge of last
    Double_t xl = (iL*fStripPitch-dx_mod/2.)*1.e3;//pl.GetStripLowerEdge (iL) * 1000.0;//mm
    Double_t xr = ((iU+1)*fStripPitch-dx_mod/2.)*1.e3;//pl.GetStripUpperEdge (iU) * 1000.0;//mm

#if DBG_AVA > 1
    cout << "iL gsle " << iL << " " << xl << endl;
    cout << "iU gsue " << iU << " " << xr << endl;
#endif
    
    //if(ipl==1 && ic<12)h1_yGEM_inava_2->Fill((xl+xr)*5.e-4);

    // Limits in y are y limits of track plus some reasonable margin
    // We do this in units of strip pitch for convenience (even though
    // this is the direction orthogonal to the pitch direction)

    // Use y-integration step size of 1/10 of strip pitch (in mm)
    Double_t yq = fStripPitch * 1000.0 / integral_steps_y;//fYIntegralStepsPerPitch;
    Double_t yb = ys0, yt = ys1;
    if (yb > yt)
      swap( yb, yt );
    yb = yq * TMath::Floor (yb / yq);
    yt = yq * TMath::Ceil  (yt / yq);
    
    // We should also allow x to have variable bin size based on the db
    // the new avalanche model (Cauchy-Lorentz) has a very sharp full width
    // half maximum, so if the bin size is too large, it can introduce
    // fairly large error on the charge deposition. Setting fXIntegralStepsPerPitch
    // to 1 will go back to the original version -- Weizhi Xiong

    Int_t nstrips = iU - iL + 1;
    Int_t nx = (iU - iL + 1) * integral_steps_x;//fXIntegralStepsPerPitch;
    Int_t ny = TMath::Nint( (yt - yb)/yq );
#if DBG_AVA > 1
    cout << "xr xl yt yb nx ny "
	 << xr << " " << xl << " " << yt << " " << yb
	 << " " << nx << " " << ny << endl;
#endif
    assert( nx > 0 && ny > 0 );
    
#if DBG_HISTOS > 0
    h2D_nplanesV_ava_nstrips->Fill(nstrips, min(ic,3)*2+ipl);
    h2D_nplanesV_ava_nx->Fill(nx, min(ic,3)*2+ipl);
#endif
    //if(ipl==0 && ic<4){
    //h1_nstrips_X[int(ic>0)]->Fill(nstrips);
    //h1_ds_X[int(ic>0)]->Fill(yt-yb);
    //}
    //if(ipl==1 && ic<4){
    //h1_nstrips_Y[int(ic>0)]->Fill(nstrips);
    //h1_ds_Y[int(ic>0)]->Fill(yt-yb);
    //}

    // define function, gaussian and sum of gaussian

    Double_t xbw = (xr - xl) / nx;
    Double_t ybw = (yt - yb) / ny;
#if DBG_AVA > 1
    cout << "xbw ybw " << xbw << " " << ybw << endl;
#endif
    
    Int_t sumASize = nx;
#if DBG_AVA > 1
    cout<<nx<<" : "<<ny<< ", nx*ny " << sumASize <<" Nstrips: "<<nstrips<<endl;
#endif
    fSumA.resize(sumASize);
    memset (&fSumA[0], 0, fSumA.size() * sizeof (Double_t));
    //Double_t fSumA[sumASize];
    //memset (fSumA, 0, sumASize * sizeof (Double_t));
    //for(Int_t i=0; i<sumASize; i++){
    //fSumA[i] = 0;
    //}
    
    //    fSumA.resize(nx*ny);
    //memset (&fSumA[0], 0, fSumA.size() * sizeof (Double_t));
#if DBG_AVA > 1
    cout << fRNIon << " " << fRIon.size() << endl;
#endif
    /*
    h1_nbins_X[ipl]->Fill(nx);
    h1_nbins_Y[ipl]->Fill(ny);
    
    h1_binw_X[ipl]->Fill(xbw);
    h1_binw_Y[ipl]->Fill(ybw);
    */
    fStart = std::chrono::steady_clock::now();
    
    Integration_semiana(roangle_mod, xl, xr, yb, yt, nx, xbw);
    
    /*
    for (UInt_t i = 0; i < fRNIon; i++){
      Double_t frxs = fRIon[i].X*cos(roangle_mod) - fRIon[i].Y*sin(roangle_mod);
      Double_t frys = fRIon[i].X*sin(roangle_mod) + fRIon[i].Y*cos(roangle_mod);

      //pl.PlaneToStrip (frxs, frys);
      //frxs *= 1e3; frys *= 1e3;
      //  cout<<"IonStrip: "<<pl.GetStrip(frxs*1e-3,frys*1e-3)<<endl;
      // bin containing center and # bins each side to process
      Int_t ix = (frxs-xl) / xbw;
      //Int_t iy = (frys-yb) / ybw;
      Int_t dx = fRIon[i].SNorm / xbw  + 1;
      //Int_t dy = fRIon[i].SNorm / ybw  + 1;
#if DBG_AVA > 1
      cout << "ix dx iy dy " << ix << " " << dx << " " << iy << " " << dy << endl;
#endif
      
      //if(ipl==1 && ic<12)h1_yGEM_inava_3->Fill(frxs*1.e-3);
       
    
      //
      // NL change:
      //
      // ggnorm is the avalance charge for the i^th ion, and R2 is the square of the radius of the diffusion 
      // circle, mutiplied by the kSNormNsigma factor: (ip.SNorm * ip.SNorm)*kSNormNsigma*kSNormNsigma. All 
      // strips falling within this circle are considered in charge summing. 
      //
      // The charge contribution to a given strip by the i^th ion is evaluated by a Lorentzian (or Gaussian)
      // distribution; the sigma for this distribution is eff_sigma, which is the actual avalance sigma. 
      //
      Double_t ggnorm = fRIon[i].ggnorm;
      Double_t r2 = fRIon[i].R2;
      Double_t eff_sigma_square = r2/(fSNormNsigma*fSNormNsigma);
      Double_t eff_sigma = TMath::Sqrt(eff_sigma_square);
      Double_t current_ion_amplitude = fAvaGain*ggnorm*(1./(TMath::Pi()*eff_sigma))*(eff_sigma*eff_sigma);
      if(isbkgd)current_ion_amplitude*=5;
      //if(ipl==0 && ic<12)h1_QnormvsX_ion->Fill(frxs*1.e-3, current_ion_amplitude);
      //if(ipl==1 && ic<12)h1_QnormvsY_ion->Fill(frxs*1.e-3, current_ion_amplitude);
      
      //h1_SigmaEff->Fill(eff_sigma);
      
      //double inte0 = 0.;
      //double intey0 = 0.;
      double inte4 = 0.;
      double ceff, yinte;
      
      // xc and yc are center of current bin
      //double sumA = 0;
      // Loop over bins
      Int_t min_binNb_x = max(ix-dx,0);
      //Int_t min_binNb_y = max(iy-dy,0);
      Int_t max_binNb_x = min(ix+dx+1,nx);
      //Int_t max_binNb_y = min(iy+dy+1,ny);
      Int_t jx = min_binNb_x;
      Double_t xc = xl + (jx+0.5) * xbw;
      for (; jx < max_binNb_x; ++jx, xc += xbw){
	Double_t xd2 = frxs-xc; xd2 *= xd2;
	//if( xd2 > r2 ){
	//  if( (xc - frxs)>0 )
	//    break;
	//  else
	//    continue;
	//}
	//Int_t jx_base = jx * ny;
	//Int_t jy = min_binNb_y;
	//Double_t yc = yb + (jy+0.5) * ybw;
	
	//cout << eff_sigma_square << " " << eff_sigma << endl;
	
	if( xd2 <= r2 ){
	  ceff = sqrt(xd2 + eff_sigma_square);
	  yinte = 1./ceff *(atan((yt-frys)/ceff) - atan((yb-frys)/ceff));
	  //inte4 += (yinte*xbw)*current_ion_amplitude;
	  fSumA[jx] += (yinte*xbw)*current_ion_amplitude;
	}
	//intey0 = 0.;

	for (; jy < max_binNb_y; ++jy, yc += ybw){
	  Double_t yd2 = frys-yc; yd2 *= yd2;
	  //if( yd2 > r2 ){
	  //  if( (yc - frys)>0 )
	  //    break;
	  //  else
	  //    continue;
	  // }
	  
	  // if(ipl==0 && ic==0)h1_QareavsX_ion->Fill(frxs*1.e-3, (xd2+yd2)+eff_sigma_square);
	  // if(ipl==1 && ic==0)h1_QareavsY_ion->Fill(frxs*1.e-3, (xd2+yd2)+eff_sigma_square);
	  
	  if( xd2 + yd2 <= r2 ) {
	    //cout << frxs-xc << " " << current_ion_amplitude << " " << xd2+yd2 << " " << eff_sigma_square << endl;
	    fSumA[jx_base+jy] += current_ion_amplitude / ((xd2+yd2)+eff_sigma_square);
	    //intey0+= current_ion_amplitude*ybw / ((xd2+yd2)+eff_sigma_square);
	    //inte0+= current_ion_amplitude*ybw / ((xd2+yd2)+eff_sigma_square);
	    //sumA+= current_ion_amplitude / ((xd2+yd2)+eff_sigma_square);
	  }
	}//cout<<endl;
	
	//h1_ava_yint->Fill(frxs-xc, (intey0-yinte*current_ion_amplitude)/(intey0+yinte*current_ion_amplitude));

      }//cout<<"##########################################################################"<<endl<<endl;getchar();

      //cout << " i_chamber " << ic << " iplane " << ipl << endl << " ava x: lower bound " << xl << " center " << frxs << " upper bound " << xr << endl << " ava y: lower bound " << yb << " center " << frys << " upper bound " << yt << endl << " nstrips " << nstrips << " min bin x " << min_binNb_x << " max bin x " << max_binNb_x << " nx " << nx << " xbw " << xbw << endl << " integral " << inte4 << endl; 

      //h1_ava_int->Fill((inte0*xbw-inte4)/(inte0*xbw+inte4));
    }
    */
    fEnd = std::chrono::steady_clock::now();
    fDiff = fEnd-fStart;
    
    fTotalTime_int+= fDiff.count();
    
#if DBG_AVA > 1
    cout << "t0 = " << t0 << " plane " << ipl 
	 << endl;
#endif

    //Int_t ai=0;
    //Double_t area = xbw * ybw;

    //when we integrate in order to get the signal pulse, we want all charge
    //deposition on the area of a single strip -- Weizhi
    
    
    for (Int_t j = 0; j < nstrips; j++){
      //  cout<<"strip: "<<iL+j<<":    ";
      //Int_t posflag = 0;
      Double_t us = 0.;
      //for (UInt_t k=0; k<fXIntegralStepsPerPitch; k++){
      for (UInt_t k=0; k<integral_steps_x; k++){
	int kx = (j * integral_steps_x + k);// * ny;
	us += fSumA[kx++];
      }
      
#if DBG_AVA > 2
      cout << "strip " << iL+j << " us " << us << endl;
#endif
      

      for (Int_t b = 0; b < fEleSamplingPoints; b++){
	Double_t pulse = PulseShape (fEleSamplingPeriod * b - t0,
				     us,
				     fPulseShapeTau);
	//fPulseShapeTau0, fPulseShapeTau1 );
	
	Short_t dadc = ADCConvert( pulse,
				   0,// fADCoffset,
				   fADCgain,
				   fADCbits );

#if DBG_AVA > 2
	if(pulse>0)
	  cout << "strip number " << iL+j << ", sampling number " << b << ", t0 = " << t0 << endl
	       << "pulse = " << pulse << ", (val - off)/gain = " 
	       << (pulse-fADCoffset)/fADCgain << ", dadc = " << dadc << endl;
#endif

	gemdet->GEMPlanes[ic*2+ipl].AddADC(iL+j, b, dadc);
	
	//cross talk here
	if(xt_factor>0){
	  if(iL+j+isLeft*fNCStripApart>=0 && iL+j+isLeft*fNCStripApart<GEMstrips){
	    gemdet->GEMPlanes[ic*2+ipl].AddADC(iL+j+isLeft*fNCStripApart, b, TMath::Nint(dadc*xt_factor));
	  }
	}
      }
      
    }//end loop on strips
    
  }//end loop on planes
}


//TGEMSBSGEMHit **
void GEMSimDigitizer::AvaModel_2(const int ic,
				 GEMDet* gemdet, 
				 TRandom3* R,
				 const TVector3& xi,
				 const TVector3& xo,
				 const Double_t t0)
{
  // xi, xo are in chamber frame, in mm
  
#if DBG_AVA > 0
  cout << " avamodel_2 " << fRNIon << " " << fRSMax << " " << fRTime0 << endl;
#endif
  Double_t nsigma = fAvalancheFiducialBand; // coverage factor
  
  //trying something:
  int integral_steps_x = 1;
  int integral_steps_y = 1;
  
  double xc_hit = (xi.X()+xo.X())/2.;
  double yc_hit = (xi.Y()+xo.Y())/2.;
  double dx2_hit = (xo.X()-xi.X())/3.0;//(xo.Z()-xi.Z());
  double dy2_hit = (xo.Y()-xi.Y())/3.0;//(xo.Z()-xi.Z());
/*Why 3.0?
  Because:
       • 	the induced charge footprint is much narrower than the full track projection
       • 	but not zero
       • 	dividing by 3 empirically matches Garfield++ and test‑beam data for GEMs
  It’s a phenomenological scaling:
  a long track → elongated charge cloud, but only by a fraction of the track length.
*/
  double NionsStrip[2] = {fRNIon*fabs(xo.X()-xi.X()), fRNIon*fabs(xo.Y()-xi.Y())};
  dx2_hit*= dx2_hit;
  dy2_hit*= dy2_hit;
  
  double rsmax = fRSMax*sqrt(1.0+dx2_hit+dy2_hit);//?
  
  Double_t x0,y0,x1,y1; // lower and upper corners of avalanche diffusion area
  
  if (xi.X()<xo.X()) {
    x0 = xi.X()-nsigma*rsmax;//*fRSMax;
    x1 = xo.X()+nsigma*rsmax;//*fRSMax;
  } else {
    x1 = xi.X()+nsigma*rsmax;//*fRSMax;
    x0 = xo.X()-nsigma*rsmax;//*fRSMax;
  }

  if (xi.Y()< xo.Y()) {
    y0 = xi.Y()-nsigma*rsmax;//*fRSMax;
    y1 = xo.Y()+nsigma*rsmax;//*fRSMax;
  } else {
    y1 = xi.Y()+nsigma*rsmax;//*fRSMax;
    y0 = xo.Y()-nsigma*rsmax;//*fRSMax;
  }
  
  //h1_AvaSizeYvsX_FastAppx->Fill(x1-x0, y1-y0);

  // Check if any part of the avalanche region is in the active area of the sector.
  // Here, "active area" means the chamber's *bounding box*, which is
  // larger than the wedge's active area (section of a ring)
  
  Double_t glx = (-gemdet->GEMPlanes[ic*2].dX())/2.*1.e3;
  Double_t gly = (-gemdet->GEMPlanes[ic*2+1].dX())/2.*1.e3;
  Double_t gux = (gemdet->GEMPlanes[ic*2].dX())/2.*1.e3;
  Double_t guy = (gemdet->GEMPlanes[ic*2+1].dX())/2.*1.e3;
  
  if(x0<glx) x0=glx;
  if(y0<gly) y0=gly;
  if(x1>gux) x1=gux;
  if(y1>guy) y1=guy;

  // Loop over chamber planes
  double roangle_mod, dx_mod, xoffset_mod;
  int GEMstrips;
  
  double xt_factor;
  int isLeft;
  for (UInt_t ipl = 0; ipl < fNROPlanes; ++ipl){
    
    xt_factor = R->Gaus(fCrossFactor, fCrossSigma);
    isLeft = R->Uniform(1.) < 0.5 ? -1 : 1;
    
    roangle_mod = -gemdet->GEMPlanes[ic*2+ipl].ROangle();
    dx_mod = gemdet->GEMPlanes[ic*2+ipl].dX();
    xoffset_mod = gemdet->GEMPlanes[ic*2+ipl].Xoffset();
    GEMstrips = gemdet->GEMPlanes[ic*2+ipl].GetNStrips();
    
    // Positions in strip frame
    Double_t xsm = ( (x0+x1)*cos(roangle_mod) - (y0+y1)*sin(roangle_mod) ) * 0.5;
    Double_t ysm = ( (x0+x1)*sin(roangle_mod) + (y0+y1)*cos(roangle_mod) ) * 0.5;
    Double_t xs0 = xsm-(x1-x0)*0.5;
    Double_t ys0 = ysm-(y1-y0)*0.5;
    Double_t xs1 = xsm+(x1-x0)*0.5;
    Double_t ys1 = ysm+(y1-y0)*0.5;
    // Double_t xs0 = x0*cos(roangle_mod) - y0*sin(roangle_mod);
    // Double_t ys0 = x0*sin(roangle_mod) + y0*cos(roangle_mod);
    // Double_t xs1 = x1*cos(roangle_mod) - y1*sin(roangle_mod); 
    // Double_t ys1 = x1*sin(roangle_mod) + y1*cos(roangle_mod);

    Int_t iL = max(0, Int_t((xs0*1.e-3+dx_mod/2.)/fStripPitch) );
    iL = min(iL, GEMstrips);
    Int_t iU = min(Int_t((xs1*1.e-3+dx_mod/2.)/fStripPitch), GEMstrips);
    iU = max(0, iU);
    
    if(iL==iU){//nothing to do
      return;
    }

    if(iU<iL)swap(iU, iL);


    //
    // Bounds of rectangular avalanche region, in strip frame
    //

    // Limits in x are low edge of first strip to high edge of last
    Double_t xl = (iL*fStripPitch-dx_mod/2.)*1.e3;
    Double_t xr = ((iU+1)*fStripPitch-dx_mod/2.)*1.e3;

#if DBG_AVA > 0
    cout << "iL gsle " << iL << " " << xl << endl;
    cout << "iU gsue " << iU << " " << xr << endl;
#endif
    

    // Limits in y are y limits of track plus some reasonable margin
    // We do this in units of strip pitch for convenience (even though
    // this is the direction orthogonal to the pitch direction)

    // Use y-integration step size of 1/10 of strip pitch (in mm)
    Double_t yq = fStripPitch * 1000.0 / integral_steps_y;//fYIntegralStepsPerPitch;
    Double_t yb = ys0, yt = ys1;
    if (yb > yt)
      swap( yb, yt );
    yb = yq * TMath::Floor (yb / yq);
    yt = yq * TMath::Ceil  (yt / yq);
    
    // We should also allow x to have variable bin size based on the db
    // the new avalanche model (Cauchy-Lorentz) has a very sharp full width
    // half maximum, so if the bin size is too large, it can introduce
    // fairly large error on the charge deposition. Setting fXIntegralStepsPerPitch
    // to 1 will go back to the original version -- Weizhi Xiong

    Int_t nstrips = iU - iL + 1;
    Int_t nx = (iU - iL + 1) * integral_steps_x;//fXIntegralStepsPerPitch;
    Int_t ny = TMath::Nint( (yt - yb)/yq );

    // define function, gaussian and sum of gaussian

    Double_t xbw = (xr - xl) / nx;
    Double_t ybw = (yt - yb) / ny;
    
    Int_t sumASize = nx * ny;
    
    fSumA.resize(sumASize);
    memset (&fSumA[0], 0, fSumA.size() * sizeof (Double_t));
    
    fStart = std::chrono::steady_clock::now();
    
    //TEST
    //double ph = R->Uniform(-TMath::Pi(), TMath::Pi());
    //double r_ = R->Gaus(0.55, 0.1);
    //dx2_hit = r_*cos(ph);dx2_hit*= dx2_hit;
    //dy2_hit = r_*sin(ph);dy2_hit*= dy2_hit;
    
    //if(sqrt(dx2_hit+dy2_hit)<1.0){
    Integration_fastappx(R, roangle_mod, NionsStrip[ipl], xc_hit, yc_hit, dx2_hit, dy2_hit, xl, xr, yb, yt, nx, xbw, ny, ybw);
    //}
    
    fEnd = std::chrono::steady_clock::now();
    fDiff = fEnd-fStart;
    
    fTotalTime_int+= fDiff.count();
    
    
    Double_t area = xbw * ybw;

    //when we integrate in order to get the signal pulse, we want all charge
    //deposition on the area of a single strip -- Weizhi
    
    for (Int_t j = 0; j < nstrips; j++){
      Double_t us = 0.;
      for (UInt_t k=0; k<integral_steps_x; k++){
	int kx = (j * integral_steps_x + k);// * ny;
	
	kx*= ny;
	Double_t integralY_tmp = 0;
	for( Int_t jy = ny; jy != 0; --jy ){
	  integralY_tmp += fSumA[kx++];
	}
	us += integralY_tmp * area;
      }
      
      for (Int_t b = 0; b < fEleSamplingPoints; b++){
	Double_t pulse = PulseShape (fEleSamplingPeriod * b - t0,
				     us,
				     fPulseShapeTau);
	//fPulseShapeTau0, fPulseShapeTau1 );
	
	Short_t dadc = ADCConvert( pulse,
				   0,// fADCoffset,
				   fADCgain,
				   fADCbits );

	gemdet->GEMPlanes[ic*2+ipl].AddADC(iL+j, b, dadc);
	
	if(xt_factor>0){
	  if(iL+j+isLeft*fNCStripApart>=0 && iL+j+isLeft*fNCStripApart<GEMstrips){
	    
	    gemdet->GEMPlanes[ic*2+ipl].AddADC(iL+j+isLeft*fNCStripApart, b, TMath::Nint(dadc*xt_factor));
	    
	  }
	}
      }

    }//end loop on strips
    
  }//end loop on planes
  
  //return virs;
}


Int_t
GEMSimDigitizer::Digitize (GEMDet* gemdet,
			   TRandom3* R, bool bkgdonly)//, 
//gmn_tree* T)
{
  // Digitize event. Add results to any existing digitized data.

  //UInt_t nh = gdata.GetNHit();
  bool is_background = false;
  Float_t event_time=0,time_zero=0;
  Double_t trigger_jitter = R->Uniform(-fAPVTimeJitter/2, fAPVTimeJitter/2);
  
  // For signal data, determine the sector of the primary track
  
  for(size_t ih = 0; ih<gemdet->fGEMhits.size(); ih++){
    is_background = (gemdet->fGEMhits[ih].source!=0);
    //cout << (bkgdonly) << " " << (!is_background) << " " << (bkgdonly && !is_background) << endl;
    //if(bkgdonly && !is_background)continue;
    UInt_t igem = gemdet->fGEMhits[ih].module;
    //UInt_t igem = iplane/2;
    
    //if(igem>=16)cout << igem << endl;
    //cout<<igem<<":"<<imodule<<":"<<iplane<<endl;
    //Short_t itype = (gdata.GetParticleType(ih)==1) ? 1 : 2; // primary = 1, secondaries = 2
    // if(gdata.GetParticleType(ih)!=1){cout<<"x"<<endl;getchar();}

    TVector3 vv1(gemdet->fGEMhits[ih].xin*1.e3,
		 gemdet->fGEMhits[ih].yin*1.e3, 
		 gemdet->fGEMhits[ih].zin*1.e3);
    TVector3 vv2(gemdet->fGEMhits[ih].xout*1.e3,
		 gemdet->fGEMhits[ih].yout*1.e3, 
		 gemdet->fGEMhits[ih].zout*1.e3);
    /*
    if(is_background){
      h1_modhit_b->Fill(igem);
      h1_xhit_b->Fill( (gemdet->fGEMhits[ih].xin+gemdet->fGEMhits[ih].xout)/2. );
      h1_yhit_b->Fill( (gemdet->fGEMhits[ih].yin+gemdet->fGEMhits[ih].yout)/2. );
      h1_zhit_b->Fill( (gemdet->fGEMhits[ih].zin+gemdet->fGEMhits[ih].zout)/2. );
      h1_xdiff_b->Fill( gemdet->fGEMhits[ih].xout-gemdet->fGEMhits[ih].xin );
      h1_ydiff_b->Fill( gemdet->fGEMhits[ih].yout-gemdet->fGEMhits[ih].yin );
      h1_edep_b->Fill( gemdet->fGEMhits[ih].edep );
    }else{
      h1_modhit_s->Fill(igem);
      h1_xhit_s->Fill( (gemdet->fGEMhits[ih].xin+gemdet->fGEMhits[ih].xout)/2. );
      h1_yhit_s->Fill( (gemdet->fGEMhits[ih].yin+gemdet->fGEMhits[ih].yout)/2. );
      h1_zhit_s->Fill( (gemdet->fGEMhits[ih].zin+gemdet->fGEMhits[ih].zout)/2. );
      h1_xdiff_s->Fill( gemdet->fGEMhits[ih].xout-gemdet->fGEMhits[ih].xin );
      h1_ydiff_s->Fill( gemdet->fGEMhits[ih].yout-gemdet->fGEMhits[ih].yin );
      h1_edep_s->Fill( gemdet->fGEMhits[ih].edep );
    }
    */
    
    if(abs(vv1.X()-vv2.X())>50 || abs(vv1.Y()-vv2.Y())>50){//in mm
      //cout<<abs(vv1.X()-vv2.X())<<endl;
      //getchar();
      continue;
    }
    //if(igem<12)h1_yGEM_preion->Fill(vv1.Y()*1.e-3);

    //test
    double dx2_hit = (vv2.X()-vv1.X())/3.0;dx2_hit*= dx2_hit;
    double dy2_hit = (vv2.Y()-vv1.Y())/3.0;dx2_hit*= dy2_hit;
    //h2D_edepVdr->Fill(sqrt(dx2_hit+dy2_hit), gemdet->fGEMhits[ih].edep);
    
    fStart = std::chrono::steady_clock::now();
    // if(!is_background){
    //fRNIon = R->Poisson(gemdet->fGEMhits[ih].edep/fGasWion);
    //fRSMax = TMath::Sqrt(2.*fGasDiffusion*TMath::Abs(fRoutZ-fEntranceRef)/fGasDriftVelocity);//maximum avalanche spread =  maximum time = maximum distance
    //fRTime0 = TMath::Abs(fRoutZ-fExitRef)/fGasDriftVelocity;//minimum time = minimum distance
    //cout << " Digitize " << gemdet->fGEMhits[ih].edep << " " << fRNIon << " " << fRSMax << " " << fRTime0 << endl;
    // }else{
    IonModel (R, vv1, vv2, gemdet->fGEMhits[ih].edep);
    //cout << " Digitize " << gemdet->fGEMhits[ih].edep << " " << fRNIon << " " << fRSMax << " " << fRTime0 << endl;
    // }
    fEnd = std::chrono::steady_clock::now();
    fDiff = fEnd-fStart;
    
    
    fTotalTime_ion+= fDiff.count();
    // Get Signal Start Time 'time_zero'
    //if( is_background ) {
    // For background data, uniformly randomize event time between
    // -fGateWidth to +75 ns (assuming 3 useful 25 ns samples).
    // Not using HitTime from simulation file but randomize HitTime to cycle use background files
    //event_time = m(-fGateWidth, 6*fEleSamplingPeriod);
    //event_time = fTimeZero;//fTrnd.Uniform(-fGateWidth/2.-fEleSamplingPeriod, fGateWidth-fEleSamplingPeriod);
    //event_time = fTrnd.Uniform((-fGateWidth+2*fEleSamplingPeriod), 8*fEleSamplingPeriod);
    //} else {
    // Signal events occur at t = 0, 
    event_time = //fTimeZero+
      gemdet->fGEMhits[ih].t;
    //}
    //  cout<<event_time<<"  "<<ih<<endl;
    // Adding drift time and trigger_jitter
    time_zero = event_time - fTriggerOffset[igem] + fRTime0*1e9 - trigger_jitter;
    
    //if( is_background ){
    //h1_thit_b->Fill(time_zero);
    //h1_nions_b->Fill(fRNIon);
    //}else{
    //h1_thit_s->Fill(time_zero);
    //h1_nions_s->Fill(fRNIon);
    //}
    
    //cout << time_zero << " " << fTimeZero << " " << gemdet->fGEMhits[ih].t 
    //<< " " << trigger_jitter << " " << fRTime0*1e9 << endl;
    
#if DBG_AVA > 0
    if(time_zero>200.0)
      cout << "time_zero " << time_zero 
	   << "; evt time " << event_time 
	   << "; hit time " << gemdet->fGEMhits[ih].t
	   << "; drift time " << fRTime0*1e9
	   << endl;
#endif
    if (fRNIon > 0) {
      //cout << "AvaModel..." << endl;
      //if(igem<12)h1_yGEM_preava->Fill(vv1.Y()*1.e-3);
      //if(!is_background)
      
      std::chrono::time_point<std::chrono::steady_clock> fStart2 = 
	std::chrono::steady_clock::now();
      if(is_background){
       	AvaModel_2 (igem, gemdet, R, vv1, vv2, time_zero);
      }else{
       	AvaModel (igem, gemdet, R, vv1, vv2, time_zero);
      }
      fEnd = std::chrono::steady_clock::now();
      fDiff = fEnd-fStart2;
      
      fTotalTime_ava+= fDiff.count();
      //cout << "Done!" << endl;
      //cout << " hou " << gemdet->GEMPlanes[4].GetADCSum(400) << endl;
      //CheckOut(gemdet, R, T);
    }
    
  }//end loop on hits
  //fFilledStrips = false;
  
  //Cumulate(gemdet, R);
    
  return 0;
}

//-------------------------------------------------------
// Helper functions for integration in AvaModel
// inline static
// Double_t IntegralY( Double_t* a, Int_t ix, Int_t nx, Int_t ny )
// {
//   register double sum = 0.;
//   register int kx = ix*ny;
//   for( Int_t jy = ny; jy != 0; --jy )
//     sum += a[kx++];

//   return sum;
// }


void GEMSimDigitizer::CheckOut(GEMDet* gemdet,
			       const int uniqueid, 
			       TRandom3* R, 
			       remoll_tree* T, 
			       bool sigonly)
			       //gmn_tree* T)
{
  //int test = gemdet->GEMPlanes[4].GetADCSum(400);
  //cout << " hou hou " << test << endl;
  // for(int j = 570; j<580; j++){
  //   cout << j << "   ";
  //   for(int b = 0; b<6; b++){
  //     cout << " " << gemdet->GEMPlanes[30].GetADC(j, b);
  //   }cout << endl;
  // }
  short strip;
  double commonmode = 0;
  if(fDoCommonMode){commonmode = fCommonModeArray[0];}
  //cout << "commonmode " << commonmode << endl;
  int apv_ctr;
  for(size_t i = 0; i<gemdet->GEMPlanes.size(); i++){
#if DBG_AVA >0 
    cout << "GEM plane/RO " << i << " ( chamber " << i/2 << ", proj " << i%2 << "); ";
#endif
    double ADC_sum = 0.0;
    int nstripshit_total = 0;
    int nstripshit_abovethr = 0;
    for(int j = 0; j<gemdet->GEMPlanes[i].GetNStrips(); j++){
      //if(gemdet->GEMPlanes[4].GetADCSum(400)!=test){
      //cout << gemdet->GEMPlanes[4].GetADCSum(400) << "!=" << test << ": " << i << " " << j << endl;
      //test = gemdet->GEMPlanes[4].GetADCSum(400);
      //}
      //cout << fDoCommonMode << endl;
      if(fDoCommonMode && !sigonly){
	//cout << " hou " << endl;
	if(j%128==0 && apv_ctr<fCommonModeArray.size()){
	  commonmode = fCommonModeArray[apv_ctr++];
	  //cout << commonmode << endl;
	}
      }
      //if(i==4 && j==400){
      //cout << gemdet->GEMPlanes[i].GetADCSum(j) << endl;
      //}
      if(gemdet->GEMPlanes[i].GetADCSum(j)>0){
#if DBG_AVA >0
#endif
	nstripshit_total++;
	ADC_sum+=gemdet->GEMPlanes[i].GetADCSum(j);
	//if(i%2==1 && i<24)h1_yGEM_incheckout->Fill(j*fStripPitch-gemdet->GEMPlanes[i].dX()/2.);
	if(!sigonly){
	  for(int k = 0; k<fNSamples; k++){
	    //cout << i << " " << j << " " << k << " " << gemdet->GEMPlanes[i].GetADC(j, k) << " " << gemdet->GEMPlanes[i].GetADCSum(j) << " = = > ";
	    //int ped = TMath::Nint(R->Gaus(commonmode, fPulseNoiseSigma));
	    //if(gemdet->GEMPlanes[i].GetADC(j, k)<0 || gemdet->GEMPlanes[i].GetADC(j, k)>4096)cout << i << " " << j << " " << k << " " << gemdet->GEMPlanes[i].GetADC(j, k) << " " << ped << " " << commonmode << " " << fPulseNoiseSigma << endl;
	    gemdet->GEMPlanes[i].AddADC(j, k, TMath::Nint(R->Gaus(commonmode, fPulseNoiseSigma)));
	    //cout << gemdet->GEMPlanes[i].GetADC(j, k) << " " << gemdet->GEMPlanes[i].GetADCSum(j)<< endl;
	    //handle saturation
	    if(gemdet->GEMPlanes[i].GetADC(j, k)>pow(2, fADCbits) ){
	      //cout << gemdet->GEMPlanes[i].GetADC(j, k) << " => ";
	      gemdet->GEMPlanes[i].SetADC(j, k, TMath::Nint(pow(2, fADCbits)) );
	      //cout << gemdet->GEMPlanes[i].GetADC(j, k) << endl;
	    }
	  }
	}
	//#ifdef 
	if( (fDoZeroSup && gemdet->GEMPlanes[i].GetADCSum(j)-commonmode*6>fZeroSup) || !fDoZeroSup) {
	  //if(i<4)cout << i << " " << gemdet->GEMPlanes[i].GetNStrips() << " " << commonmode << endl;
	  if(uniqueid==0){
	    /*
	    if(sigonly){
	      cout << T->mollergem_dig_sig.nstrips << endl;
	      for(int k = 0; k<fNSamples; k++){
		T->mollergem_dig_sig.nstrips++;
		T->mollergem_dig_sig.module->push_back(i);
		T->mollergem_dig_sig.strip->push_back(j);
	      
	      //T->mollergem_dig_sig.adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      //T->mollergem_dig_sig.adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      //T->mollergem_dig_sig.adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      //T->mollergem_dig_sig.adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      //T->mollergem_dig_sig.adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      //T->mollergem_dig_sig.adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
		T->mollergem_dig_sig.samp->push_back(k);
		T->mollergem_dig_sig.adc->push_back(gemdet->GEMPlanes[i].GetADC(j, k));
	      }

	    }else{
	      */
	      for(int k = 0; k<fNSamples; k++){
		T->mollergem_dig.nstrips++;
		T->mollergem_dig.module->push_back(i);
		T->mollergem_dig.strip->push_back(j);
	      /*
	      T->mollergem_dig.adc_0->push_back(gemdet->GEMPlanes[i].GetADC(j, 0));
	      T->mollergem_dig.adc_1->push_back(gemdet->GEMPlanes[i].GetADC(j, 1));
	      T->mollergem_dig.adc_2->push_back(gemdet->GEMPlanes[i].GetADC(j, 2));
	      T->mollergem_dig.adc_3->push_back(gemdet->GEMPlanes[i].GetADC(j, 3));
	      T->mollergem_dig.adc_4->push_back(gemdet->GEMPlanes[i].GetADC(j, 4));
	      T->mollergem_dig.adc_5->push_back(gemdet->GEMPlanes[i].GetADC(j, 5));
	      */
		T->mollergem_dig.samp->push_back(k);
		T->mollergem_dig.adc->push_back(gemdet->GEMPlanes[i].GetADC(j, k));
	      }
	      //}
	  }
	  
	  
	}//end if(...)
	
      }
    }
#if DBG_HISTOS > 0
    if(nstripshit_total>0){
      if(i>9){
	h2D_nplanesVnActiveStrips->Fill(nstripshit_abovethr, i-4);
	h2D_nplanesVnAllHitStrips->Fill(nstripshit_total, i-4);
	h2D_nplanesVnADCSum->Fill(ADC_sum, i-4);
      }else if(i>7){
	h2D_nplanesVnActiveStrips->Fill(nstripshit_abovethr, i-2);
	h2D_nplanesVnAllHitStrips->Fill(nstripshit_total, i-2);
	h2D_nplanesVnADCSum->Fill(ADC_sum, i-2);
      }else{
	h2D_nplanesVnActiveStrips->Fill(nstripshit_abovethr, i);
	h2D_nplanesVnAllHitStrips->Fill(nstripshit_total, i);
	h2D_nplanesVnADCSum->Fill(ADC_sum, i);
      }
    }
#endif
#if DBG_AVA >0 
    if(nstripshit_total>0){cout << " N hit strips (above thr): " << nstripshit_abovethr << " (total) " << nstripshit_total << " ADCsum " << ADC_sum << endl;}else{cout << endl;}
#endif 
  }
}

//___________________________________________________________________________________
void GEMSimDigitizer::Print()
{
  cout << "GEM digitization:" << endl;
  cout << "  Gas parameters:" << endl;
  cout << "    Gas ion width: " << fGasWion << endl;
  cout << "    Gas diffusion: " << fGasDiffusion << endl;
  cout << "    Gas drift velocity: " << fGasDriftVelocity << endl;
  cout << "    Avalanche fiducial band: " << fAvalancheFiducialBand << endl;
  cout << "    Avalanche charge statistics: " << fAvalancheChargeStatistics << endl;
  cout << "    Gain mean: " << fGainMean << endl;
  cout << "    Gain 0: " << fGain0 << endl;

  cout << "  Electronics parameters:" << endl;
  //cout << "    Trigger offsets: "; //<< fTriggerOffset 
  //for(int i = 0; i<fManager->GetNChamber(); i++)cout << fTriggerOffset[i] << " ";
  //cout << endl;
  cout << "    APV time jitter: " << fAPVTimeJitter << endl;
  cout << "    Sampling Period: " << fEleSamplingPeriod << endl;
  cout << "    Sampling Points: " << fEleSamplingPoints   << endl;
  cout << "    Pulse Noise width: " << fPulseNoiseSigma << endl;
  cout << "    ADC offset: " << fADCoffset << endl;
  cout << "    ADC gain: " << fADCgain << endl;
  cout << "    ADC bits: " << fADCbits << endl;
  cout << "    Gate width: " << fGateWidth << endl;

  cout << "  Pulse shaping parameters:" << endl;
  cout << "    Pulse shape tau: " << fPulseShapeTau << endl;
  //cout << "    Pulse shape tau1: " << fPulseShapeTau1 << endl;
}

void GEMSimDigitizer::write_histos()
{
#if DBG_HISTOS > 0
  h2D_nplanesV_ava_dx->Write();
  h2D_nplanesV_ava_dxs->Write();
  h2D_nplanesV_ava_dy->Write();
  h2D_nplanesV_ava_dys->Write();
  h2D_nplanesV_ava_nstrips->Write();
  h2D_nplanesV_ava_nx->Write();
  h2D_nplanesVnActiveStrips->Write();
  h2D_nplanesVnAllHitStrips->Write();
  h2D_nplanesVnADCSum->Write();

  /*
  h2D_edepVdr->Write();

  h1_AvaSizeYvsX_SemiAna->Write();
  h1_AvaSizeYvsX_FastAppx->Write();
  h1_AvaSizeVsZion_SemiAna->Write();
  h1_AvaSizeVsTTime_SemiAna->Write();
  
  h1_SumweightsSemiAna->Write();
  h2D_SumweightsFastAppx->Write();
  
  h1_GammaEffSemiAna->Write();
  h2D_GammaEffFastAppx->Write();
  
  h1_NormSemiAna->Write();
  h1_NormFastAppx->Write();
  
  h1_SigmaEff->Write();
  h1_NionsPix->Write();
  
  for(int i = 0; i<2; i++){
    h1_nbins_X[i]->Write();
    h1_nbins_Y[i]->Write();
    h1_binw_X[i]->Write();
    h1_binw_Y[i]->Write();
  }
  
  h1_ava_yint->Write();
  h1_ava_int->Write();
  
  h1_modhit_s->Write();
  h1_xhit_s->Write();
  h1_yhit_s->Write();
  h1_zhit_s->Write();
  h1_xdiff_s->Write();
  h1_ydiff_s->Write();
  h1_thit_s->Write();
  h1_edep_s->Write();
  h1_nions_s->Write();
  
  h1_modhit_b->Write();
  h1_xhit_b->Write();
  h1_yhit_b->Write();
  h1_zhit_b->Write();
  h1_xdiff_b->Write();
  h1_ydiff_b->Write();
  h1_thit_b->Write();
  h1_edep_b->Write();
  h1_nions_b->Write();
  
  for(int i = 0; i<2; i++){
    h1_nstrips_X[i]->Write();
    h1_nstrips_Y[i]->Write();
    
    h1_ds_X[i]->Write();
    h1_ds_Y[i]->Write();
    
    h1_nbins_X[i]->Write();
    h1_nbins_Y[i]->Write();
    
    h1_fSumA_X[i]->Write();
    h1_fSumA_Y[i]->Write();
  }
  
  h1_QvsX_ion->Write();
  h1_QvsY_ion->Write();
  h1_QnormvsX_ion->Write();
  h1_QnormvsY_ion->Write();
  h1_QareavsX_ion->Write();
  h1_QareavsY_ion->Write();
  h1_QintvsX_ion->Write();
  h1_QintvsY_ion->Write();

  h1_QvsX_ava->Write();
  h1_QvsY_ava->Write();
  h1_QintYvsX_ava->Write();
  h1_QintYvsY_ava->Write();
  
  h1_yGEM_preion->Write();
  h1_yGEM_preava->Write();
  h1_yGEM_inava->Write();
  h1_yGEM_inava_2->Write();
  h1_yGEM_inava_3->Write();
  h1_yGEM_inava_4->Write();
  h1_xGEMvsADC_inava_4->Write();
  h1_yGEMvsADC_inava_4->Write();
  h1_yGEM_incheckout->Write();
  */
#endif
}

void GEMSimDigitizer::print_time_execution()
{
  cout << " Ionization total duration " << std::setprecision(9) << fTotalTime_ion << " s;" << endl;
  cout << " Avalanche total duration " << std::setprecision(9) << fTotalTime_ava << " s;" << endl;
  cout << " Integration total duration " << std::setprecision(9) << fTotalTime_int << " s;" << endl;
}

/*
Double_t
TGEMSBSSimDigitization::CommonMode(UInt_t i_mpd)
{
  if(fCommonModeArray.size() && fDoCommonMode){
    i_mpd = (i_mpd<fCommonModeArray.size() ? i_mpd: 0);
    return fCommonModeArray[i_mpd];
  }else{
    return 0;
  }
}

Double_t
TGEMSBSSimDigitization::ZeroSupThreshold(UInt_t i_mpd)
{
  //i_mpd
  if(fDoZeroSup){
    if(fDoCommonMode){
      return fZeroSup+CommonMode(i_mpd)*fEleSamplingPoints;
    }else{
      return fZeroSup;
    }
  }else{
    return -1000;
  }
}
*/

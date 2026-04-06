#include "mollerDigPMTSignal.h"

#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <limits>

mollerDigPMTSignal::mollerDigPMTSignal()
  : fDetID(-1)
{
  Clear(true);
}

mollerDigPMTSignal::mollerDigPMTSignal(int detid)
  : fDetID(detid)
{
  Clear(true);
}

void mollerDigPMTSignal::Clear(bool clearSamples)
{
  fEdepMeV = 0.0;
  fHitTimesNS.clear();
  fEarliestTimeNS = 0.0;
  fHasEarliestTime = false;

  fWX = 0.0;
  fWY = 0.0;
  fWZ = 0.0;
  fWSum = 0.0;
  fPosX = 0.0;
  fPosY = 0.0;
  fPosZ = 0.0;
  fHasPosition = false;

  fLeffCM = 0.0;
  fMeanPE = 0.0;
  fNPEPoiss = 0;
  fADCInt = 0;
  fADCIntPedSub = 0;
  if (clearSamples) {
    fADCSamples.clear();
  }
}

void mollerDigPMTSignal::AddSumEdep(double edep_mev)
{
  if (edep_mev > 0.0) {
    fEdepMeV += edep_mev;
  }
}

void mollerDigPMTSignal::AddHit(double t_ns, double x, double y, double z, double edep_mev)
{
  fHitTimesNS.push_back(t_ns);

  if (!fHasEarliestTime || t_ns < fEarliestTimeNS) {
    fEarliestTimeNS = t_ns;
    fHasEarliestTime = true;
  }

  double w = std::max(edep_mev, 0.0);
  if (w == 0.0) {
    w = 1.0e-12;
  }

  fWX += w * x;
  fWY += w * y;
  fWZ += w * z;
  fWSum += w;

  if (fWSum > 0.0) {
    fPosX = fWX / fWSum;
    fPosY = fWY / fWSum;
    fPosZ = fWZ / fWSum;
    fHasPosition = true;
  }
}

double mollerDigPMTSignal::FrankTammPhotonsPerMeter(double n, double beta,
                                                    double lambda1_m, double lambda2_m)
{
  const double alpha = 1.0 / 137.036;
  if (beta <= 0.0) return 0.0;
  const double term = 1.0 - 1.0 / (beta * beta * n * n);
  if (term <= 0.0) return 0.0;
  return 2.0 * M_PI * alpha * term * (1.0 / lambda1_m - 1.0 / lambda2_m);
}

double mollerDigPMTSignal::QLSBCoulomb(const MollerDigFADCConfig& c)
{
  const double Vlsb = c.v_range_volt / std::pow(2.0, c.nbits);
  const double dt_s = c.dt_ns * 1.0e-9;
  return (Vlsb * dt_s) / c.R_ohm;
}

double mollerDigPMTSignal::ChargeInSampleFromOnePE(const MollerDigFADCConfig& c,
                                                   double t0_ns,
                                                   double t1_ns,
                                                   double t2_ns)
{
  if (t2_ns <= t0_ns) return 0.0;
  const double a = std::max(t1_ns, t0_ns);
  const double b = t2_ns;
  return c.Qpe_C * (std::exp(-(a - t0_ns) / c.tau_ns) -
                    std::exp(-(b - t0_ns) / c.tau_ns));
}

void mollerDigPMTSignal::Digitize(TRandom3& rng,
                                  const MollerDigOpticsConfig& optics,
                                  const MollerDigFADCConfig& fadc)
{
  fLeffCM = std::max(0.0, fEdepMeV / optics.dEdx_mev_per_cm);

  const double yield_per_m = FrankTammPhotonsPerMeter(optics.n_quartz,
                                                      optics.beta,
                                                      optics.lambda1_m,
                                                      optics.lambda2_m);
  const double yield_per_cm = 0.01 * yield_per_m;
  const double mean_phot = yield_per_cm * fLeffCM;
  fMeanPE = optics.eps_photon_to_pe * mean_phot;
  fNPEPoiss = rng.Poisson(fMeanPE);

  const int nsamp = std::max(0, int(std::lround(fadc.gate_ns / fadc.dt_ns)));
  fADCSamples.assign(nsamp, 0);
  fADCInt = 0;
  fADCIntPedSub = 0;

  if (fNPEPoiss <= 0 || nsamp <= 0) {
    return;
  }

  const double t_ref = fHasEarliestTime ? fEarliestTimeNS : 0.0;
  const double qlsb = QLSBCoulomb(fadc);
  const int adc_max = (1 << fadc.nbits) - 1;
  const double tmin = -0.5 * fadc.gate_ns + fadc.t_offset_ns;

  std::vector<double> q_samp(nsamp, 0.0);

  for (int ipe = 0; ipe < fNPEPoiss; ++ipe) {
    const double tpe = t_ref + rng.Gaus(0.0, fadc.sigma_time_ns);
    for (int is = 0; is < nsamp; ++is) {
      const double t1 = tmin + is * fadc.dt_ns;
      const double t2 = t1 + fadc.dt_ns;
      q_samp[is] += ChargeInSampleFromOnePE(fadc, tpe, t1, t2);
    }
  }

  for (int is = 0; is < nsamp; ++is) {
    double adc_f = q_samp[is] / qlsb;
    double ped = static_cast<double>(fadc.pedestal_mean);
    if (fadc.pedestal_sigma > 0.0) {
      ped += rng.Gaus(0.0, fadc.pedestal_sigma);
    }
    adc_f += ped;

    int adc_i = static_cast<int>(std::llround(adc_f));
    adc_i = std::max(0, std::min(adc_i, adc_max));

    fADCSamples[is] = static_cast<unsigned short>(adc_i);
    fADCInt += adc_i;
    fADCIntPedSub += (adc_i - fadc.pedestal_mean);
  }
}

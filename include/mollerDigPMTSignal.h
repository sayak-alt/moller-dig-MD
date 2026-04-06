#ifndef MollerDigPMTSignal_H
#define MollerDigPMTSignal_H

#include <vector>
#include <cstdint>

class TRandom3;

struct MollerDigFADCConfig {
  double R_ohm = 50.0;
  double dt_ns = 4.0;
  int    nbits = 12;
  double v_range_volt = 1.0;

  int    pedestal_mean = 300;
  double pedestal_sigma = 5.0;

  double tau_ns = 10.0;
  double sigma_time_ns = 1.0;

  double gate_ns = 400.0;
  double t_offset_ns = 30.0;

  double Qpe_C = 1.6e-13;
};

struct MollerDigOpticsConfig {
  double n_quartz = 1.47;
  double beta = 1.0;
  double lambda1_m = 200e-9;
  double lambda2_m = 600e-9;
  double eps_photon_to_pe = 0.03;
  double dEdx_mev_per_cm = 4.0;
};

class mollerDigPMTSignal {
public:
  mollerDigPMTSignal();
  explicit mollerDigPMTSignal(int detid);
  virtual ~mollerDigPMTSignal() = default;

  void SetDetID(int detid) { fDetID = detid; }
  int  GetDetID() const { return fDetID; }

  void Clear(bool clearSamples = true);

  void AddSumEdep(double edep_mev);
  void AddHit(double t_ns, double x, double y, double z, double edep_mev);

  void Digitize(TRandom3& rng,
                const MollerDigOpticsConfig& optics,
                const MollerDigFADCConfig& fadc);

  // accessors for event-level bookkeeping
  double GetEdepMeV() const { return fEdepMeV; }
  double GetLeffCM() const { return fLeffCM; }
  double GetMeanPE() const { return fMeanPE; }
  int    GetNPEPoiss() const { return fNPEPoiss; }
  int    GetADCInt() const { return fADCInt; }
  int    GetADCIntPedSub() const { return fADCIntPedSub; }
  double GetEarliestTime() const { return fEarliestTimeNS; }
  bool   HasEarliestTime() const { return fHasEarliestTime; }
  double GetPosX() const { return fPosX; }
  double GetPosY() const { return fPosY; }
  double GetPosZ() const { return fPosZ; }
  bool   HasPosition() const { return fHasPosition; }

  const std::vector<double>& GetHitTimes() const { return fHitTimesNS; }
  const std::vector<unsigned short>& GetADCSamples() const { return fADCSamples; }

  static double FrankTammPhotonsPerMeter(double n, double beta,
                                         double lambda1_m, double lambda2_m);
  static double QLSBCoulomb(const MollerDigFADCConfig& c);
  static double ChargeInSampleFromOnePE(const MollerDigFADCConfig& c,
                                        double t0_ns,
                                        double t1_ns,
                                        double t2_ns);

private:
  int fDetID;

  // accumulated raw event info
  double fEdepMeV;
  std::vector<double> fHitTimesNS;
  double fEarliestTimeNS;
  bool   fHasEarliestTime;

  double fWX;
  double fWY;
  double fWZ;
  double fWSum;
  double fPosX;
  double fPosY;
  double fPosZ;
  bool   fHasPosition;

  // digitized outputs
  double fLeffCM;
  double fMeanPE;
  int    fNPEPoiss;
  int    fADCInt;
  int    fADCIntPedSub;
  std::vector<unsigned short> fADCSamples;
};

#endif

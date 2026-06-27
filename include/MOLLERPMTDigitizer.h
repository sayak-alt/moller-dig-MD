#ifndef MOLLER_PMT_DIGITIZER_H
#define MOLLER_PMT_DIGITIZER_H

#include <TTree.h>
#include <TRandom3.h>

#include <vector>
#include <unordered_map>

class MOLLERPMTDigitizer {
public:
  struct Config {
    double R_ohm = 50.0;
    double dt_ns = 4.0;
    int nbits = 12;
    double v_range_volt = 1.0;

    int pedestal_mean = 300;
    double pedestal_sigma = 5.0;

    double tau_ns = 10.0;
    double sigma_time_ns = 1.0;

    double gate_ns = 400.0;
    double t_offset_ns = 30.0;

    double Qpe_C = 1.6e-13;

    double n_quartz = 1.47;
    double beta = 1.0;
    double lambda_min_nm = 200.0;
    double lambda_max_nm = 600.0;
    double eps_photon_to_pe = 0.03;
    double dEdx_mev_per_cm = 4.0;

    bool store_waveforms = false;
  };

  MOLLERPMTDigitizer();

  void SetConfig(const Config& cfg);
  void SetRandomSeed(unsigned int seed);
  void SetEventNumber(int iev);

  void BookBranches(TTree* tree);
  void Clear();

  void DigitizeEvent(const std::vector<int>& hit_det,
                     const std::vector<double>& hit_t,
                     const std::vector<double>& hit_x,
                     const std::vector<double>& hit_y,
                     const std::vector<double>& hit_z,
                     const std::vector<double>& hit_edep,
                     const std::vector<int>& sum_det,
                     const std::vector<double>& sum_edep,
                     bool has_sum);

private:
  Config fCfg;
  TRandom3 fRand;
  int fEvnum = 0;

  std::vector<int> detid;
  std::vector<double> edep_mev;
  std::vector<double> leff_cm;
  std::vector<double> meanpe;
  std::vector<int> npe_poiss;

  std::vector<int> adc_int;
  std::vector<int> adc_int_pedsub;

  std::vector<int> hit_detid;
  std::vector<double> hit_time;

  std::vector<int> t0_detid;
  std::vector<double> t0_time;

  // Optional waveform output
  std::vector<int> wf_detid;
  std::vector<int> wf_samp;
  std::vector<unsigned short> wf_adc;

  double FrankTammPhotonsPerCm() const;
  double QLSB() const;
  double ChargeInSampleFromOnePE(double t0_ns, double t1_ns, double t2_ns) const;
};

#endif

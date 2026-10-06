#ifndef MOLLER_SHOWERMAX_DIGITIZER_H
#define MOLLER_SHOWERMAX_DIGITIZER_H

#include <TTree.h>
#include <TRandom3.h>

#include "MOLLERShowerMaxLookupResponse.h"
#include <string>
#include <vector>
#include <unordered_map>

class MOLLERShowerMaxDigitizer {
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

    bool store_waveforms = false;
    std::string detector_map_file; // Optional CSV electronics map
    std::string lut_directory = "showermax_tables";
    bool legacy_pion_names = true;
  };

  MOLLERShowerMaxDigitizer();

  void SetConfig(const Config& cfg);
  void SetRandomSeed(unsigned int seed);
  void SetEventNumber(int iev);
  // Input rate in Hz; divisor follows the analysis file normalization.
  void SetEventRate(double rate_hz, double divisor = 1.0);

  void BookBranches(TTree* tree);
  void Clear();

  void DigitizeEvent(const std::vector<MOLLERShowerMaxLookupResponse::Hit>& hits);

private:
  MOLLERShowerMaxLookupResponse fShowerMaxLookup;
  void DigitizePrepared(const std::unordered_map<int,double>& energy,
                        const std::unordered_map<int,double>& times,
                        const std::unordered_map<int,double>* lut_means = nullptr);
  struct ChannelInfo {
    int segment, ring, rocid, slot, channel;
    std::string segment_group, orientation, subdivision;
    bool enabled;
  };
  void LoadDetectorMap(const std::string& filename);
  void AppendMapping(int did);
  std::unordered_map<int,ChannelInfo> fDetectorMap;
  std::vector<int> rocid, slot, channel, segment, ring, mapping_valid, channel_enabled;
  std::vector<std::string> segment_group, orientation, subdivision;
  std::vector<int> wf_rocid, wf_slot, wf_channel, wf_mapping_valid;
  Config fCfg;
  TRandom3 fRand;
  int fEvnum = 0;
  double event_rate_hz = 0.0, rate_GHz = 0.0;
  std::vector<double> rate_x_meanpe, rate_x_npe;

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

  double QLSB() const;
  double ChargeInSampleFromOnePE(double t0_ns, double t1_ns, double t2_ns) const;
};

#endif

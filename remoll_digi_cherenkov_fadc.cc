//      - hit level timing vectors per detid (all hit times, and earliest hit time)
//      - poisson smeared npe actually used (npe_poiss)
//      - pedestal subtracted integrated ADC (adc_int_pedsub)
//      - effective length per tile per event (leff_cm) used in PE model
//      - meanpe (pre-Poisson expectation) per detid per event
//
// Notes:
//  - Frank–Tamm photons/cm computed for 200–450 nm, n=1.47, beta~1.
//  - Convert Edep(MeV) -> L_eff(cm) by L_eff = Edep / (dE/dx). dEdx_mev_per_cm = ~ 4.0 (google)!
//  - One-number photon->PE eps_photon_to_pe = 0.02 (need to tune).
//  - FADC is charge-integrating: ADC sample counts = Q_sample / Q_LSB + pedestal.
//    Q_LSB = (V_LSB * dt / R), with V_LSB = Vrange/2^nbits.
//
// Compile (bash):
//   g++ -O2 -std=c++17 remoll_digi_fulloptics.cpp $(root-config --cflags --libs) -o remoll_digi_fulloptics
//
// Run:
//   ./remoll_digi_fulloptics input.root output.root T
// -----------------------------------------------------------------------------

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// -----------------------------
// Frank–Tamm yield integrated over wavelength band [λ1, λ2]
// dN/dx = 2π α (1 - 1/(β^2 n^2)) (1/λ1 - 1/λ2)
// λ in meters, x in meters
// -----------------------------
static inline double frank_tamm_photons_per_meter(double n, double beta,
                                                  double lambda1_m, double lambda2_m)
{
  const double alpha = 1.0 / 137.036;
  if (beta <= 0) return 0.0;
  double term = 1.0 - 1.0 / (beta * beta * n * n);
  if (term <= 0.0) return 0.0;
  return 2.0 * M_PI * alpha * term * (1.0 / lambda1_m - 1.0 / lambda2_m);
}

// -----------------------------
// Charge-integrating FADC model
// -----------------------------
struct FADCConfig {
  double R_ohm          = 50.0;
  double dt_ns          = 4.0;
  int    nbits          = 12;
  double v_range_volt   = 1.0;

  int    pedestal_mean  = 300;     // in ADC counts
  double pedestal_sigma = 5.0;   // in ADC counts

  double tau_ns         = 10.0;  // simple exponential decay constant
  double sigma_time_ns  = 1.0;   // PE time jitter

  double gate_ns        = 400.0; // 1 us
  double t_offset_ns    = 30.0;   // (need to check)

  double Qpe_C          = 1.6e-13; // PE charge at anode (assuming gain of the pmt is 1e6)
};

static inline double qlsb_coulomb(const FADCConfig& c)
{
  const double Vlsb = c.v_range_volt / std::pow(2.0, c.nbits);
  const double dt_s = c.dt_ns * 1e-9;
  return (Vlsb * dt_s) / c.R_ohm;
}

static inline double charge_in_sample_from_one_pe(const FADCConfig& c,
                                                  double t0_ns,
                                                  double t1_ns, double t2_ns)
{
  // Pulse current i(t) = (Qpe/tau)*exp(-(t-t0)/tau) for t>=t0
  // Charge in sample [t1,t2] is Qpe*(exp(-(a-t0)/tau) - exp(-(b-t0)/tau)), with a=max(t1,t0)
  if (t2_ns <= t0_ns) return 0.0;
  const double tau = c.tau_ns;
  const double a = std::max(t1_ns, t0_ns);
  const double b = t2_ns;
  return c.Qpe_C * (std::exp(-(a - t0_ns) / tau) - std::exp(-(b - t0_ns) / tau));
}

int main(int argc, char** argv)
{
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " input.root output.root treeName\n";
    return 1;
  }

  const char* inname  = argv[1];
  const char* outname = argv[2];
  const char* tname   = argv[3];

  TFile fin(inname, "READ");
  if (fin.IsZombie()) {
    std::cerr << "ERROR: cannot open input " << inname << "\n";
    return 1;
  }

  TTree* Tin = dynamic_cast<TTree*>(fin.Get(tname));
  if (!Tin) {
    std::cerr << "ERROR: cannot find tree '" << tname << "' in " << inname << "\n";
    return 1;
  }

  // sum branches exist?
  bool has_sum = (Tin->GetBranch("sum.det") && Tin->GetBranch("sum.edep"));

  // Reader
  TTreeReader reader(Tin);

  // hit branches (timing + position, and fallback edep if sum missing)
  TTreeReaderArray<int>    hit_det (reader, "hit.det");
  TTreeReaderArray<double> hit_t   (reader, "hit.t");
  TTreeReaderArray<double> hit_x   (reader, "hit.x");
  TTreeReaderArray<double> hit_y   (reader, "hit.y");
  TTreeReaderArray<double> hit_z   (reader, "hit.z");
  TTreeReaderArray<double> hit_edep(reader, "hit.edep");

  // sum branches (preferred for edep)
  std::unique_ptr<TTreeReaderArray<int>>    sum_det;
  std::unique_ptr<TTreeReaderArray<double>> sum_edep;
  if (has_sum) {
    sum_det  = std::make_unique<TTreeReaderArray<int>>(reader, "sum.det");
    sum_edep = std::make_unique<TTreeReaderArray<double>>(reader, "sum.edep");
  }

  // Output
  TFile fout(outname, "RECREATE");
  TTree Tout("digi", "remoll -> libsbsdig-like PMT digitization (FADC charge-integrating)");

  int evnum = 0;

  // Per-detid (tile) quantities (only for detids that appear in the event)
  std::vector<int>    detid;
  std::vector<double> edep_mev;
  std::vector<double> leff_cm;
  std::vector<double> meanpe;
  std::vector<int>    npe_poiss;

  // Integrated ADC (raw sum of ADC samples including pedestal) and pedestal-subtracted
  std::vector<int> adc_int;
  std::vector<int> adc_int_pedsub;

  // Sample-wise output (libsbsdig style): one entry per (detid, sample)
  std::vector<int>            samp_detid;
  std::vector<int>            samp;
  std::vector<unsigned short> adc; // 0..4095

  // Position estimate per detid (energy-weighted from hits)
  std::vector<int>    pos_detid;
  std::vector<double> pos_x;
  std::vector<double> pos_y;
  std::vector<double> pos_z;

  // Hit-level timing per detid: store all hit times and earliest hit time
  std::vector<int>    hit_detid;
  std::vector<double> hit_time;       // parallel to hit_detid, one entry per hit we keep
  std::vector<int>    t0_detid;
  std::vector<double> t0_time;        // earliest hit time per detid (from hit branch)

  // Branches
  Tout.Branch("evnum", &evnum, "evnum/I");

  Tout.Branch("detid",    &detid);
  Tout.Branch("edep_mev", &edep_mev);
  Tout.Branch("leff_cm",  &leff_cm);
  Tout.Branch("meanpe",   &meanpe);
  Tout.Branch("npe_poiss",&npe_poiss);

  Tout.Branch("adc_int",        &adc_int);
  Tout.Branch("adc_int_pedsub", &adc_int_pedsub);

  Tout.Branch("samp_detid", &samp_detid);
  Tout.Branch("samp",       &samp);
  Tout.Branch("adc",        &adc);

  Tout.Branch("pos_detid", &pos_detid);
  Tout.Branch("pos_x",     &pos_x);
  Tout.Branch("pos_y",     &pos_y);
  Tout.Branch("pos_z",     &pos_z);

  Tout.Branch("hit_detid", &hit_detid);
  Tout.Branch("hit_time",  &hit_time);

  Tout.Branch("t0_detid",  &t0_detid);
  Tout.Branch("t0_time",   &t0_time);

  // Physics / optics knobs
  const double n_quartz = 1.47;
  const double beta = 1.0;
  const double lambda1 = 200e-9;
  const double lambda2 = 600e-9;
  const double eps_photon_to_pe = 0.03;

  // Edep -> Leff proxy via dE/dx (tune!)
  const double dEdx_mev_per_cm = 4.0;

  const double yield_per_m  = frank_tamm_photons_per_meter(n_quartz, beta, lambda1, lambda2);
  const double yield_per_cm = yield_per_m * 0.01;

  FADCConfig fadc;
  const int nsamp = int(std::lround(fadc.gate_ns / fadc.dt_ns)); // 250 for 1us/4ns
  const double qlsb = qlsb_coulomb(fadc);
  const int adc_max = (1 << fadc.nbits) - 1;

  TRandom3 R(0);

  std::cout << "Input entries: " << Tin->GetEntries() << "\n";
  std::cout << "has_sum = " << (has_sum ? "true" : "false") << "\n";
  std::cout << "Frank–Tamm yield ~ " << yield_per_cm
            << " photons/cm in 200–450 nm (n=" << n_quartz << ", beta~1)\n";
  std::cout << "ADC: " << fadc.nbits << " bits, " << fadc.v_range_volt << " V range, "
            << fadc.dt_ns << " ns bins, R=" << fadc.R_ohm
            << " ohm => Q_LSB=" << qlsb * 1e15 << " fC\n";

  // Event loop
  while (reader.Next()) {
    evnum++;

    // Clear outputs
    detid.clear();
    edep_mev.clear();
    leff_cm.clear();
    meanpe.clear();
    npe_poiss.clear();
    adc_int.clear();
    adc_int_pedsub.clear();
    samp_detid.clear();
    samp.clear();
    adc.clear();
    pos_detid.clear();
    pos_x.clear();
    pos_y.clear();
    pos_z.clear();
    hit_detid.clear();
    hit_time.clear();
    t0_detid.clear();
    t0_time.clear();
    
    if (evnum%1000 == 0) std::cout << " event no " << evnum << " out of " << Tin->GetEntries() << std::endl;

    // --- Build per-detid Edep (MeV) using sum branch (preferred) or hits fallback
    std::unordered_map<int, double> edep_by_det;

    if (has_sum) {
      for (auto i = 0u; i < sum_det->GetSize(); i++) {
        int did = (*sum_det)[i];
        double e = (*sum_edep)[i];
        if (e <= 0) continue;
        edep_by_det[did] += e;
      }
    } else {
      for (auto i = 0u; i < hit_det.GetSize(); i++) {
        int did = hit_det[i];
        double e = hit_edep[i];
        if (e <= 0) continue;
        edep_by_det[did] += e;
      }
    }

    // --- Hit-level timing and position estimates from hit branch
    // Store all hit times (as requested) and also build earliest time + pos weights
    std::unordered_map<int, double> t0_by_det;
    std::unordered_map<int, double> wx_by_det, wy_by_det, wz_by_det, wsum_by_det;

    // initialize t0 to +inf when first seen
    for (auto i = 0u; i < hit_det.GetSize(); i++) {
      int did = hit_det[i];
      double t = hit_t[i];

      // store hit-level timing (all hits)
      hit_detid.push_back(did);
      hit_time.push_back(t);

      // earliest time
      auto it = t0_by_det.find(did);
      if (it == t0_by_det.end()) {
        t0_by_det[did] = t;
      } else {
        if (t < it->second) it->second = t;
      }

      // energy-weighted position using hit.edep as weight (fallback tiny if 0)
      double w = std::max(hit_edep[i], 0.0);
      if (w == 0.0) w = 1e-12;
      wx_by_det[did] += w * hit_x[i];
      wy_by_det[did] += w * hit_y[i];
      wz_by_det[did] += w * hit_z[i];
      wsum_by_det[did] += w;
    }

    // Build compact earliest-time vectors
    t0_detid.reserve(t0_by_det.size());
    t0_time.reserve(t0_by_det.size());
    for (const auto& kv : t0_by_det) {
      t0_detid.push_back(kv.first);
      t0_time.push_back(kv.second);
    }

    // --- For each detid in edep map: compute meanpe, Poisson npe, waveform, ADC
    // We only digitize detids with edep > 0 (from sum/hit fallback)
    detid.reserve(edep_by_det.size());
    edep_mev.reserve(edep_by_det.size());
    leff_cm.reserve(edep_by_det.size());
    meanpe.reserve(edep_by_det.size());
    npe_poiss.reserve(edep_by_det.size());
    adc_int.reserve(edep_by_det.size());
    adc_int_pedsub.reserve(edep_by_det.size());

    for (const auto& kv : edep_by_det) {
      const int did = kv.first;
      const double e_mev = kv.second;
      if (e_mev <= 0) continue;

      // L_eff from Edep proxy
      double L_eff_cm = e_mev / dEdx_mev_per_cm;
      if (L_eff_cm < 0) L_eff_cm = 0;

      // Frank–Tamm photons and PE expectation
      const double mean_phot = yield_per_cm * L_eff_cm;
      const double mean_pe   = eps_photon_to_pe * mean_phot;

      // Poisson smear
      const int npe = R.Poisson(mean_pe);
      if (npe <= 0) {
        // Still store bookkeeping (so you can see meanpe/leff even if npe==0)
        detid.push_back(did);
        edep_mev.push_back(e_mev);
        leff_cm.push_back(L_eff_cm);
        meanpe.push_back(mean_pe);
        npe_poiss.push_back(0);
        adc_int.push_back(0);
        adc_int_pedsub.push_back(0);
        continue;
      }

      // Reference time = earliest hit time if available, else 0
      double t_ref = 0.0;
      auto it0 = t0_by_det.find(did);
      if (it0 != t0_by_det.end() && std::isfinite(it0->second)) t_ref = it0->second;

      // Build waveform as charge per 4ns sample
      std::vector<double> q_samp(nsamp, 0.0);

      // Gate window [tmin, tmin+gate]
      const double tmin = -0.5 * fadc.gate_ns + fadc.t_offset_ns;

      // Generate npe photoelectrons with time jitter and add contributions
      for (int k = 0; k < npe; k++) {
        const double tpe = t_ref + R.Gaus(0.0, fadc.sigma_time_ns);
        for (int is = 0; is < nsamp; is++) {
          const double t1 = tmin + is * fadc.dt_ns;
          const double t2 = t1 + fadc.dt_ns;
          q_samp[is] += charge_in_sample_from_one_pe(fadc, tpe, t1, t2);
        }
      }

      // Convert to ADC counts (charge-integrating); store samples and integrated sums
      int adc_sum_raw = 0;
      int adc_sum_pedsub = 0;

      for (int is = 0; is < nsamp; is++) {
        double adc_f = q_samp[is] / qlsb;

        // add pedestal/noise
        double ped = double(fadc.pedestal_mean);
        if (fadc.pedestal_sigma > 0.0) ped += R.Gaus(0.0, fadc.pedestal_sigma);
        adc_f += ped;

        int adc_i = (int)std::llround(adc_f);
        if (adc_i < 0) adc_i = 0;
        if (adc_i > adc_max) adc_i = adc_max;

        adc_sum_raw += adc_i;

        // pedestal-subtracted sample (use mean pedestal only; noise stays in raw)
        int adc_pedsub_i = adc_i - fadc.pedestal_mean;
        adc_sum_pedsub += adc_pedsub_i;

        // sample-wise output
        samp_detid.push_back(did);
        samp.push_back(is);
        adc.push_back((unsigned short)adc_i);
      }

      // Store per-detid outputs
      detid.push_back(did);
      edep_mev.push_back(e_mev);
      leff_cm.push_back(L_eff_cm);
      meanpe.push_back(mean_pe);
      npe_poiss.push_back(npe);
      adc_int.push_back(adc_sum_raw);
      adc_int_pedsub.push_back(adc_sum_pedsub);

      // Position estimate if available
      auto itw = wsum_by_det.find(did);
      if (itw != wsum_by_det.end() && itw->second > 0.0) {
        pos_detid.push_back(did);
        pos_x.push_back(wx_by_det[did] / itw->second);
        pos_y.push_back(wy_by_det[did] / itw->second);
        pos_z.push_back(wz_by_det[did] / itw->second);
      }
    }

    Tout.Fill();
  }

  fout.Write();
  fout.Close();
  fin.Close();

  std::cout << "Wrote output: " << outname << "\n";
  return 0;
}


#include "MOLLERPMTDigitizer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <set>
#include <tuple>

MOLLERPMTDigitizer::MOLLERPMTDigitizer()
{
  fRand.SetSeed(0);
}

void MOLLERPMTDigitizer::SetConfig(const Config& cfg)
{
  if(cfg.pe_method != 0 && cfg.pe_method != 1)
    throw std::runtime_error("pe_method_mollerpmt must be 0 or 1");
  if(cfg.pe_method == 1) fLookup.Load(cfg.lut_directory, cfg.lut_ring);
  LoadDetectorMap(cfg.detector_map_file);
  fCfg = cfg;
}

void MOLLERPMTDigitizer::SetRandomSeed(unsigned int seed)
{
  fRand.SetSeed(seed);
}

void MOLLERPMTDigitizer::SetEventNumber(int iev)
{
  fEvnum = iev;
}

double MOLLERPMTDigitizer::FrankTammPhotonsPerCm() const
{
  const double alpha = 1.0 / 137.036;

  const double l1 = fCfg.lambda_min_nm * 1.0e-9;
  const double l2 = fCfg.lambda_max_nm * 1.0e-9;

  if(fCfg.beta <= 0.0) return 0.0;

  double term =
    1.0 - 1.0 / (fCfg.beta * fCfg.beta * fCfg.n_quartz * fCfg.n_quartz);

  if(term <= 0.0) return 0.0;

  double yield_per_m =
    2.0 * M_PI * alpha * term * (1.0 / l1 - 1.0 / l2);

  return yield_per_m * 0.01;
}

double MOLLERPMTDigitizer::QLSB() const
{
  const double Vlsb = fCfg.v_range_volt / std::pow(2.0, fCfg.nbits);
  const double dt_s = fCfg.dt_ns * 1.0e-9;
  return (Vlsb * dt_s) / fCfg.R_ohm;
}

double MOLLERPMTDigitizer::ChargeInSampleFromOnePE(
    double t0_ns, double t1_ns, double t2_ns) const
{
  if(t2_ns <= t0_ns) return 0.0;

  const double tau = fCfg.tau_ns;
  const double a = std::max(t1_ns, t0_ns);
  const double b = t2_ns;

  return fCfg.Qpe_C *
         (std::exp(-(a - t0_ns) / tau) -
          std::exp(-(b - t0_ns) / tau));
}

void MOLLERPMTDigitizer::BookBranches(TTree* tree)
{
  if(!tree) return;

  tree->Branch("evnum", &fEvnum, "evnum/I");

  tree->Branch("rocid", &rocid);
  tree->Branch("slot", &slot);
  tree->Branch("channel", &channel);
  tree->Branch("segment", &segment);
  tree->Branch("ring", &ring);
  tree->Branch("segment_group", &segment_group);
  tree->Branch("orientation", &orientation);
  tree->Branch("subdivision", &subdivision);
  tree->Branch("mapping_valid", &mapping_valid);
  tree->Branch("channel_enabled", &channel_enabled);
  tree->Branch("detid",     &detid);
  tree->Branch("edep_mev",  &edep_mev);
  tree->Branch("leff_cm",   &leff_cm);
  tree->Branch("meanpe",    &meanpe);
  tree->Branch("npe_poiss", &npe_poiss);

  tree->Branch("adc_int",        &adc_int);
  tree->Branch("adc_int_pedsub", &adc_int_pedsub);

  tree->Branch("hit_detid", &hit_detid);
  tree->Branch("hit_time",  &hit_time);

  tree->Branch("t0_detid", &t0_detid);
  tree->Branch("t0_time",  &t0_time);

  if(fCfg.store_waveforms) {
    tree->Branch("wf_rocid", &wf_rocid);
    tree->Branch("wf_slot", &wf_slot);
    tree->Branch("wf_channel", &wf_channel);
    tree->Branch("wf_mapping_valid", &wf_mapping_valid);
    tree->Branch("wf_detid", &wf_detid);
    tree->Branch("wf_samp",  &wf_samp);
    tree->Branch("wf_adc",   &wf_adc);
  }
}

void MOLLERPMTDigitizer::Clear()
{
  rocid.clear(); slot.clear(); channel.clear(); segment.clear(); ring.clear();
  segment_group.clear(); orientation.clear(); subdivision.clear();
  mapping_valid.clear(); channel_enabled.clear();
  wf_rocid.clear(); wf_slot.clear(); wf_channel.clear(); wf_mapping_valid.clear();
  detid.clear();
  edep_mev.clear();
  leff_cm.clear();
  meanpe.clear();
  npe_poiss.clear();

  adc_int.clear();
  adc_int_pedsub.clear();

  hit_detid.clear();
  hit_time.clear();

  t0_detid.clear();
  t0_time.clear();

  wf_detid.clear();
  wf_samp.clear();
  wf_adc.clear();
}

void MOLLERPMTDigitizer::DigitizeEvent(
    const std::vector<int>& hit_det,
    const std::vector<double>& hit_t,
    const std::vector<double>& hit_x,
    const std::vector<double>& hit_y,
    const std::vector<double>& hit_z,
    const std::vector<double>& hit_edep,
    const std::vector<int>& sum_det,
    const std::vector<double>& sum_edep,
    bool has_sum)
{
  Clear();

  std::unordered_map<int, double> edep_by_det;

  if(has_sum) {
    for(size_t i = 0; i < sum_det.size(); i++) {
      if(sum_edep[i] <= 0.0) continue;
      edep_by_det[sum_det[i]] += sum_edep[i];
    }
  } else {
    for(size_t i = 0; i < hit_det.size(); i++) {
      if(hit_edep[i] <= 0.0) continue;
      edep_by_det[hit_det[i]] += hit_edep[i];
    }
  }

  std::unordered_map<int, double> t0_by_det;

  for(size_t i = 0; i < hit_det.size(); i++) {
    int did = hit_det[i];

    hit_detid.push_back(did);
    hit_time.push_back(hit_t[i]);

    if(!t0_by_det.count(did) || hit_t[i] < t0_by_det[did])
      t0_by_det[did] = hit_t[i];
  }

  for(const auto& kv : t0_by_det) {
    t0_detid.push_back(kv.first);
    t0_time.push_back(kv.second);
  }

  DigitizePrepared(edep_by_det, t0_by_det);
}

void MOLLERPMTDigitizer::DigitizePrepared(
    const std::unordered_map<int,double>& edep_by_det,
    const std::unordered_map<int,double>& t0_by_det,
    const std::unordered_map<int,double>* lut_means)
{
  const double yield_per_cm = FrankTammPhotonsPerCm();
  const int nsamp = static_cast<int>(std::lround(fCfg.gate_ns / fCfg.dt_ns));
  const double qlsb = QLSB();
  const int adc_max = (1 << fCfg.nbits) - 1;

  for(const auto& kv : edep_by_det) {
    const int did = kv.first;
    const double e_mev = kv.second;

    double L_eff_cm = e_mev / fCfg.dEdx_mev_per_cm;
    if(L_eff_cm < 0.0) L_eff_cm = 0.0;

    const double mean_phot = yield_per_cm * L_eff_cm;
    const double mean_pe = lut_means ? lut_means->at(did)
                                    : fCfg.eps_photon_to_pe * mean_phot;
    const int npe = fRand.Poisson(mean_pe);

    detid.push_back(did);
    AppendMapping(did);
    edep_mev.push_back(e_mev);
    leff_cm.push_back(L_eff_cm);
    meanpe.push_back(mean_pe);
    npe_poiss.push_back(npe);

    if(npe <= 0) {
      adc_int.push_back(0);
      adc_int_pedsub.push_back(0);
      continue;
    }

    double t_ref = 0.0;
    if(t0_by_det.count(did) && std::isfinite(t0_by_det.at(did)))
      t_ref = t0_by_det.at(did);

    std::vector<double> q_samp(nsamp, 0.0);
    const double tmin = -0.5 * fCfg.gate_ns + fCfg.t_offset_ns;

    for(int k = 0; k < npe; k++) {
      const double tpe = t_ref + fRand.Gaus(0.0, fCfg.sigma_time_ns);

      for(int is = 0; is < nsamp; is++) {
        const double t1 = tmin + is * fCfg.dt_ns;
        const double t2 = t1 + fCfg.dt_ns;
        q_samp[is] += ChargeInSampleFromOnePE(tpe, t1, t2);
      }
    }

    int adc_sum_raw = 0;
    int adc_sum_pedsub = 0;

    for(int is = 0; is < nsamp; is++) {
      double adc_f = q_samp[is] / qlsb;

      double ped = static_cast<double>(fCfg.pedestal_mean);
      if(fCfg.pedestal_sigma > 0.0)
        ped += fRand.Gaus(0.0, fCfg.pedestal_sigma);

      adc_f += ped;

      int adc_i = static_cast<int>(std::llround(adc_f));
      if(adc_i < 0) adc_i = 0;
      if(adc_i > adc_max) adc_i = adc_max;

      adc_sum_raw += adc_i;
      adc_sum_pedsub += adc_i - fCfg.pedestal_mean;

      if(fCfg.store_waveforms) {
        wf_rocid.push_back(rocid.back());
        wf_slot.push_back(slot.back());
        wf_channel.push_back(channel.back());
        wf_mapping_valid.push_back(mapping_valid.back());
        wf_detid.push_back(did);
        wf_samp.push_back(is);
        wf_adc.push_back(static_cast<unsigned short>(adc_i));
      }
    }

    adc_int.push_back(adc_sum_raw);
    adc_int_pedsub.push_back(adc_sum_pedsub);
  }
}

void MOLLERPMTDigitizer::DigitizeEventLUT(
    const std::vector<MOLLERLookupResponse::Hit>& hits,
    const std::vector<int>& sum_det,
    const std::vector<double>& sum_edep, bool has_sum)
{
  if(fCfg.pe_method != 1) throw std::runtime_error("LUT mode is not configured");
  Clear();
  std::unordered_map<int,double> means, times, energy, deposited;
  if(has_sum) {
    if(sum_det.size()!=sum_edep.size()) throw std::runtime_error("Mismatched sum vectors");
    for(size_t i=0;i<sum_det.size();++i) if(sum_edep[i]>0) deposited[sum_det[i]]+=sum_edep[i];
  }
  for(const auto& hit:hits) {
    double pe=fLookup.MeanPE(hit);
    if(pe<0) continue;
    means[hit.det]+=pe;
    if(!times.count(hit.det)||hit.t<times[hit.det]) times[hit.det]=hit.t;
    hit_detid.push_back(hit.det); hit_time.push_back(hit.t);
    if(has_sum) energy[hit.det]=deposited[hit.det];
    else { energy[hit.det]+=std::max(hit.edep,0.); }
  }
  for(const auto& kv:times) {t0_detid.push_back(kv.first);t0_time.push_back(kv.second);}
  DigitizePrepared(energy,times,&means);
  // Effective path length is undefined for LUT PE generation.
  for(auto& length:leff_cm) length=std::numeric_limits<double>::quiet_NaN();
}

void MOLLERPMTDigitizer::LoadDetectorMap(const std::string& filename)
{
  std::unordered_map<int,ChannelInfo> loaded;
  if(filename.empty()) { fDetectorMap.clear(); return; }
  std::ifstream input(filename);
  if(!input) throw std::runtime_error("Cannot open detector map: "+filename);
  std::string line;
  if(!std::getline(input,line)) throw std::runtime_error("Empty detector map: "+filename);
  if(!line.empty() && line.back()=='\r') line.pop_back();
  if(line!="detid,segment,segment_group,orientation,ring,subdivision,rocid,slot,channel,enabled")
    throw std::runtime_error("Unexpected detector-map header: "+filename);
  std::set<std::tuple<int,int,int>> addresses;
  size_t row=1;
  while(std::getline(input,line)) {
    ++row;
    if(line.find_first_not_of(" \t\r")==std::string::npos) continue;
    if(!line.empty() && line.back()=='\r') line.pop_back();
    std::stringstream ss(line); std::string token; std::vector<std::string> v;
    while(std::getline(ss,token,',')) v.push_back(token);
    auto fail=[&]() { return std::runtime_error("Invalid detector map row "+std::to_string(row)+" in "+filename); };
    if(v.size()!=10) throw fail();
    auto number=[&](size_t i) {
      size_t used=0; int n;
      try { n=std::stoi(v[i],&used); } catch(...) {throw fail();}
      if(used!=v[i].size()) throw fail();
      return n;
    };
    int did=number(0), seg=number(1), r=number(4), roc=number(6), sl=number(7), ch=number(8), on=number(9);
    if(did<=0 || seg<1 || seg>14 || r<1 || r>6 || roc<0 || sl<1 || sl>21 || ch<0 || ch>15 || (on!=0 && on!=1)) throw fail();
    if((v[2]!="FF" && v[2]!="BF") || (v[3]!="FF" && v[3]!="BF")) throw fail();
    if(v[5]!="single" && v[5]!="left" && v[5]!="centre" && v[5]!="right") throw fail();
    if(loaded.count(did)) throw std::runtime_error("Duplicate detector ID in "+filename);
    if(on && !addresses.insert({roc,sl,ch}).second) throw std::runtime_error("Duplicate enabled hardware address in "+filename);
    loaded.emplace(did,ChannelInfo{seg,r,roc,sl,ch,v[2],v[3],v[5],on!=0});
  }
  if(loaded.empty()) throw std::runtime_error("No detectors in map: "+filename);
  fDetectorMap.swap(loaded);
}

void MOLLERPMTDigitizer::AppendMapping(int did)
{
  const auto it=fDetectorMap.find(did);
  bool valid=it!=fDetectorMap.end();
  mapping_valid.push_back(valid?1:0);
  channel_enabled.push_back(valid ? (it->second.enabled?1:0) : -1);
  rocid.push_back(valid?it->second.rocid:-1);
  slot.push_back(valid?it->second.slot:-1);
  channel.push_back(valid?it->second.channel:-1);
  segment.push_back(valid?it->second.segment:-1);
  ring.push_back(valid?it->second.ring:-1);
  segment_group.push_back(valid?it->second.segment_group:"unknown");
  orientation.push_back(valid?it->second.orientation:"unknown");
  subdivision.push_back(valid?it->second.subdivision:"unknown");
}

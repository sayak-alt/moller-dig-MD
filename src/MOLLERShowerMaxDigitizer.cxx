#include "MOLLERShowerMaxDigitizer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <set>
#include <tuple>

MOLLERShowerMaxDigitizer::MOLLERShowerMaxDigitizer()
{
  fRand.SetSeed(0);
}

void MOLLERShowerMaxDigitizer::SetConfig(const Config& cfg)
{
  fShowerMaxLookup.Load(cfg.lut_directory,cfg.legacy_pion_names);
  LoadDetectorMap(cfg.detector_map_file);
  fCfg = cfg;
}

void MOLLERShowerMaxDigitizer::SetRandomSeed(unsigned int seed)
{
  fRand.SetSeed(seed);
}

void MOLLERShowerMaxDigitizer::SetEventNumber(int iev)
{
  fEvnum = iev;
}

double MOLLERShowerMaxDigitizer::QLSB() const
{
  const double Vlsb = fCfg.v_range_volt / std::pow(2.0, fCfg.nbits);
  const double dt_s = fCfg.dt_ns * 1.0e-9;
  return (Vlsb * dt_s) / fCfg.R_ohm;
}

double MOLLERShowerMaxDigitizer::ChargeInSampleFromOnePE(
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

void MOLLERShowerMaxDigitizer::BookBranches(TTree* tree)
{
  if(!tree) return;

  tree->Branch("evnum", &fEvnum, "evnum/I");
  tree->Branch("event_rate_hz", &event_rate_hz, "event_rate_hz/D");
  tree->Branch("rate_GHz", &rate_GHz, "rate_GHz/D");
  tree->Branch("rate_x_meanpe", &rate_x_meanpe);
  tree->Branch("rate_x_npe", &rate_x_npe);

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

void MOLLERShowerMaxDigitizer::Clear()
{
  rocid.clear(); slot.clear(); channel.clear(); segment.clear(); ring.clear();
  segment_group.clear(); orientation.clear(); subdivision.clear();
  mapping_valid.clear(); channel_enabled.clear();
  wf_rocid.clear(); wf_slot.clear(); wf_channel.clear(); wf_mapping_valid.clear();
  detid.clear();
  edep_mev.clear();
  leff_cm.clear();
  rate_x_meanpe.clear(); rate_x_npe.clear();
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

void MOLLERShowerMaxDigitizer::DigitizePrepared(
    const std::unordered_map<int,double>& edep_by_det,
    const std::unordered_map<int,double>& t0_by_det,
    const std::unordered_map<int,double>* lut_means)
{
  const int nsamp = static_cast<int>(std::lround(fCfg.gate_ns / fCfg.dt_ns));
  const double qlsb = QLSB();
  const int adc_max = (1 << fCfg.nbits) - 1;

  for(const auto& kv : edep_by_det) {
    const int did = kv.first;
    const double e_mev = kv.second;

    const double L_eff_cm = std::numeric_limits<double>::quiet_NaN();
    const double mean_pe = lut_means->at(did);
    const int npe = fRand.Poisson(mean_pe);

    detid.push_back(did);
    AppendMapping(did);
    edep_mev.push_back(e_mev);
    leff_cm.push_back(L_eff_cm);
    meanpe.push_back(mean_pe);
    npe_poiss.push_back(npe);
    rate_x_meanpe.push_back(rate_GHz * mean_pe);
    rate_x_npe.push_back(rate_GHz * npe);

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

void MOLLERShowerMaxDigitizer::LoadDetectorMap(const std::string& filename)
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
    if(did<=0 || seg<1 || seg>28 || r<1 || r>7 || roc<0 || sl<1 || sl>21 || ch<0 || ch>15 || (on!=0 && on!=1)) throw fail();
    if((v[2]!="FF" && v[2]!="BF" && v[2]!="NA") || (v[3]!="FF" && v[3]!="BF" && v[3]!="NA")) throw fail();
    if(v[5]!="single" && v[5]!="left" && v[5]!="centre" && v[5]!="right") throw fail();
    if(loaded.count(did)) throw std::runtime_error("Duplicate detector ID in "+filename);
    if(on && !addresses.insert({roc,sl,ch}).second) throw std::runtime_error("Duplicate enabled hardware address in "+filename);
    loaded.emplace(did,ChannelInfo{seg,r,roc,sl,ch,v[2],v[3],v[5],on!=0});
  }
  if(loaded.empty()) throw std::runtime_error("No detectors in map: "+filename);
  fDetectorMap.swap(loaded);
}

void MOLLERShowerMaxDigitizer::AppendMapping(int did)
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

void MOLLERShowerMaxDigitizer::DigitizeEvent(const std::vector<MOLLERShowerMaxLookupResponse::Hit>& hits)
{
 Clear();std::unordered_map<int,double> means,times,energy;
 for(const auto& h:hits){double pe=fShowerMaxLookup.MeanPE(h);if(pe<0)continue;
  means[h.det]+=pe;energy[h.det]=0;
  if(!times.count(h.det)||h.t<times[h.det])times[h.det]=h.t;
  hit_detid.push_back(h.det);hit_time.push_back(h.t);
 }
 for(const auto& v:times){t0_detid.push_back(v.first);t0_time.push_back(v.second);}
 DigitizePrepared(energy,times,&means);
 for(auto& e:edep_mev)e=std::numeric_limits<double>::quiet_NaN();
 for(auto& l:leff_cm)l=std::numeric_limits<double>::quiet_NaN();
}

void MOLLERShowerMaxDigitizer::SetEventRate(double rate_hz, double divisor)
{
  if(!std::isfinite(rate_hz) || !std::isfinite(divisor) || divisor <= 0.0)
    throw std::runtime_error("Invalid event rate or rate normalization divisor");
  event_rate_hz = rate_hz;
  rate_GHz = rate_hz / 1.0e9 / divisor;
}

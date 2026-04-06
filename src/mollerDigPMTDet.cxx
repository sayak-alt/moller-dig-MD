#include "mollerDigPMTDet.h"

#include <TTree.h>
#include <TRandom3.h>

#include <algorithm>
#include <stdexcept>

mollerDigPMTDet::mollerDigPMTDet()
  : fUniqueID(0), fEvNum(0)
{
}

mollerDigPMTDet::mollerDigPMTDet(unsigned short uniqueid)
  : fUniqueID(uniqueid), fEvNum(0)
{
}

mollerDigPMTSignal& mollerDigPMTDet::GetOrCreateSignal(int detid)
{
  auto it = fSignals.find(detid);
  if (it == fSignals.end()) {
    it = fSignals.emplace(detid, mollerDigPMTSignal(detid)).first;
  }
  return it->second;
}

void mollerDigPMTDet::ClearEvent(bool clearSignals)
{
  if (clearSignals) {
    fSignals.clear();
  } else {
    for (auto& kv : fSignals) {
      kv.second.Clear(true);
    }
  }

  fOutData.ClearBranches();
}

void mollerDigPMTDet::AddHit(int detid, double t_ns, double x, double y, double z, double edep_mev)
{
  GetOrCreateSignal(detid).AddHit(t_ns, x, y, z, edep_mev);
}

void mollerDigPMTDet::AddSum(int detid, double edep_mev)
{
  GetOrCreateSignal(detid).AddSumEdep(edep_mev);
}

void mollerDigPMTDet::LoadEvent(const std::vector<int>& hit_det,
                                const std::vector<double>& hit_t,
                                const std::vector<double>& hit_x,
                                const std::vector<double>& hit_y,
                                const std::vector<double>& hit_z,
                                const std::vector<double>& hit_edep,
                                const std::vector<int>* sum_det,
                                const std::vector<double>* sum_edep)
{
  ClearEvent(true);

  const std::size_t nhit = hit_det.size();
  if (hit_t.size() != nhit || hit_x.size() != nhit || hit_y.size() != nhit ||
      hit_z.size() != nhit || hit_edep.size() != nhit) {
    throw std::runtime_error("mollerDigPMTDet::LoadEvent: hit arrays have inconsistent sizes");
  }

  for (std::size_t i = 0; i < nhit; ++i) {
    AddHit(hit_det[i], hit_t[i], hit_x[i], hit_y[i], hit_z[i], hit_edep[i]);
  }

  if (sum_det && sum_edep) {
    if (sum_det->size() != sum_edep->size()) {
      throw std::runtime_error("mollerDigPMTDet::LoadEvent: sum arrays have inconsistent sizes");
    }
    for (std::size_t i = 0; i < sum_det->size(); ++i) {
      if ((*sum_edep)[i] > 0.0) {
        AddSum((*sum_det)[i], (*sum_edep)[i]);
      }
    }
  } else {
    // fallback: sum hit.edep exactly as in your standalone code
    for (std::size_t i = 0; i < nhit; ++i) {
      if (hit_edep[i] > 0.0) {
        AddSum(hit_det[i], hit_edep[i]);
      }
    }
  }
}

void mollerDigPMTDet::Digitize(TRandom3& rng)
{
  for (auto& kv : fSignals) {
    kv.second.Digitize(rng, fOptics, fFADC);
  }
}

bool mollerDigPMTDet::SetupOutputTree(TTree* tree, const char* prefix)
{
  if (!tree) {
    throw std::runtime_error("mollerDigPMTDet::SetupOutputTree: null TTree pointer");
  }
  return fOutData.SetupBranches(tree, prefix);
}

void mollerDigPMTDet::FillOutputBranches()
{
  fOutData.ClearBranches();

  for (const auto& kv : fSignals) {
    const int detid = kv.first;
    const mollerDigPMTSignal& sig = kv.second;

    for (double t : sig.GetHitTimes()) {
      fOutData.hit_detid->push_back(detid);
      fOutData.hit_time->push_back(t);
    }

    if (sig.HasEarliestTime()) {
      fOutData.t0_detid->push_back(detid);
      fOutData.t0_time->push_back(sig.GetEarliestTime());
    }

    if (sig.GetEdepMeV() <= 0.0) {
      continue;
    }

    fOutData.detid->push_back(detid);
    fOutData.edep_mev->push_back(sig.GetEdepMeV());
    fOutData.leff_cm->push_back(sig.GetLeffCM());
    fOutData.meanpe->push_back(sig.GetMeanPE());
    fOutData.npe_poiss->push_back(sig.GetNPEPoiss());
    fOutData.adc_int->push_back(sig.GetADCInt());
    fOutData.adc_int_pedsub->push_back(sig.GetADCIntPedSub());

    if (sig.HasPosition()) {
      fOutData.pos_detid->push_back(detid);
      fOutData.pos_x->push_back(sig.GetPosX());
      fOutData.pos_y->push_back(sig.GetPosY());
      fOutData.pos_z->push_back(sig.GetPosZ());
    }

    const auto& samples = sig.GetADCSamples();
    for (std::size_t is = 0; is < samples.size(); ++is) {
      fOutData.samp_detid->push_back(detid);
      fOutData.samp->push_back(static_cast<int>(is));
      fOutData.adc->push_back(static_cast<int>(samples[is]));
    }
  }

  fOutData.FillBranches();
}

#ifndef MollerDigPMTDet_H
#define MollerDigPMTDet_H

#include "mollerDigPMTSignal.h"
#include "mollerDigOutDataPMT.h"

#include <map>
#include <memory>
#include <vector>

class TTree;
class TRandom3;

class mollerDigPMTDet {
public:
  mollerDigPMTDet();
  explicit mollerDigPMTDet(unsigned short uniqueid);
  virtual ~mollerDigPMTDet() = default;

  void SetUniqueID(unsigned short uniqueid) { fUniqueID = uniqueid; }
  unsigned short GetUniqueID() const { return fUniqueID; }

  MollerDigOpticsConfig& Optics() { return fOptics; }
  const MollerDigOpticsConfig& Optics() const { return fOptics; }

  MollerDigFADCConfig& FADC() { return fFADC; }
  const MollerDigFADCConfig& FADC() const { return fFADC; }

  void ClearEvent(bool clearSignals = true);

  void AddHit(int detid, double t_ns, double x, double y, double z, double edep_mev);
  void AddSum(int detid, double edep_mev);

  // Convenience helper for remoll-style event content.
  void LoadEvent(const std::vector<int>& hit_det,
                 const std::vector<double>& hit_t,
                 const std::vector<double>& hit_x,
                 const std::vector<double>& hit_y,
                 const std::vector<double>& hit_z,
                 const std::vector<double>& hit_edep,
                 const std::vector<int>* sum_det = nullptr,
                 const std::vector<double>* sum_edep = nullptr);

  void Digitize(TRandom3& rng);

  // ROOT tree interface in the same prefix-based style as mollerDigOutData.h
  bool SetupOutputTree(TTree* tree, const char* prefix = "md");
  void FillOutputBranches();
  DigPMTData_t& GetOutData() { return fOutData; }
  const DigPMTData_t& GetOutData() const { return fOutData; }

  int GetEventNumber() const { return fEvNum; }
  void SetEventNumber(int evnum) { fEvNum = evnum; }

private:
  mollerDigPMTSignal& GetOrCreateSignal(int detid);

private:
  unsigned short fUniqueID;
  int fEvNum;

  MollerDigOpticsConfig fOptics;
  MollerDigFADCConfig   fFADC;

  std::map<int, mollerDigPMTSignal> fSignals;

  DigPMTData_t fOutData;
};

#endif

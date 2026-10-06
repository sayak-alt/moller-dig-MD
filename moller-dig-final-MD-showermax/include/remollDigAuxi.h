#ifndef REMOLL_DIG_AUXI_H
#define REMOLL_DIG_AUXI_H

#include <map>
#include <vector>
#include <TRandom3.h>

#include "GEMDet.h"

// Forward declare
class remollGenericDetectorHit_t;
/*
struct remollSimHit {
    int det;
    int pid;
    int trid;
    int mtrid;

    double x, y, z;
    double px, py, pz;
    double edep;
};
*/
//class remollDigAuxi {
//public:
    bool UnfoldEvent(const std::vector<remollGenericDetectorHit_t>* hits, TRandom3* R,
	std::vector<GEMDet*> gemdets, std::vector<int> gemmap,
	double tzero, int signal);

 //   const std::vector<remollSimHit>& GetHits(int detID) const;
   // std::vector<int> GetActiveDetectors() const;

//    void Clear();
/*
private:
    std::map<int, std::vector<remollSimHit>> fHitsByDet;
    static const std::vector<remollSimHit> fEmpty;
*/
//};

#endif

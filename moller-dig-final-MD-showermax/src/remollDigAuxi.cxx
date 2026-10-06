#include "remollDigAuxi.h"
#include "remollGenericDetectorHit.hh"
#include <iostream>

//const std::vector<remollSimHit> remollDigAuxi::fEmpty;

//void remollDigAuxi::Clear() {
 //   fHitsByDet.clear();
//}

bool UnfoldEvent(
    const std::vector<remollGenericDetectorHit_t>* hits, TRandom3* R,
	std::vector<GEMDet*> gemdets, std::vector<int> gemmap, 
	double tzero, int signal) {

    //Clear();
    bool has_data = false;
    if (!hits) return has_data;
     
    int chan, mod;

    //for (const auto& h : *hits) {
    for (size_t ihit=0; ihit<hits->size(); ihit++) {
	//has_data=false;
        //remollSimHit remollh; //remoll hit
        const auto&  h = hits->at(ihit);
        if(h.det>=300 && h.det<=327){
	GEMDet::gemhit simdh; 

	mod=0;
	simdh.source = signal;
	simdh.module = h.det-300;
	simdh.edep = h.edep*1.e6;
	simdh.t = h.t;
        // remoll stores local coords in mm; digitizer expects meters (*1e3 → mm)
        // yl is the radial distance from the beam axis (not module-centered).
        // Subtract the module y-center so the digitizer sees coordinates near 0.
        // plane12 solid (det 300-313): y range ~583-1065mm, center ~824mm
        // plane34 solid (det 314-327): y range ~612-1096mm, center ~854mm
        const double yl_center = (h.det < 314) ? 824.0 : 854.0; // mm
        simdh.xin = h.xl * 1.e-3;
        simdh.yin = (h.yl - yl_center) * 1.e-3;
        simdh.zin = h.zl * 1.e-3;
        // remoll provides one hit position; extrapolate exit using local momentum
        // direction through the GEM drift gap (3 mm = 0.003 m)
        const double gem_dz = 0.003; // m
        if(fabs(h.pzl) > 0){ 
          simdh.xout = simdh.xin + (h.pxl / fabs(h.pzl)) * gem_dz;
          simdh.yout = simdh.yin + (h.pyl / fabs(h.pzl)) * gem_dz;
          simdh.zout = simdh.zin + gem_dz;
        } else {
          simdh.xout = simdh.xin;
          simdh.yout = simdh.yin;
          simdh.zout = simdh.zin;
        }
	//std::cout<<"Has_data = "<<has_data<<" ; Det = "<<h.det<<" , pid = "<<h.pid<<" ; xl = "<<h.xl<<std::endl;
	gemdets[gemmap[0]]->fGEMhits.push_back(simdh);
	
	has_data = true;
	//std::cout<<"Has_data = "<<has_data<<" ; Det = "<<h.det<<" , source = "<<simdh.source<<" ; module = "<<simdh.module<<" edep = "<<simdh.edep<<" ; t = "<<simdh.t<<
//	" ; xin = "<<simdh.xin	<<" ;yin = "<<simdh.yin<<" ;zin = "<<simdh.zin<<
//	" ; xout = "<<simdh.xout	<<" ;yout = "<<simdh.yout<<" ;zout = "<<simdh.zout<<
//	std::endl;
	}
	else continue;
	/*
        remollh.det   = h.det;
        remollh.pid   = h.pid;
        remollh.trid  = h.trid;
        remollh.mtrid = h.mtrid;

        remollh.x = h.x;
        remollh.y = h.y;
        remollh.z = h.z;

        remollh.px = h.px;
        remollh.py = h.py;
        remollh.pz = h.pz;

        remollh.edep = h.edep;

        fHitsByDet[remollh.det].push_back(remollh);
*/
    }
    return has_data;
}
/*
const std::vector<remollSimHit>&
remollDigAuxi::GetHits(int detID) const {
    auto it = fHitsByDet.find(detID);
    if (it == fHitsByDet.end()) return fEmpty;
    return it->second;
}

std::vector<int> remollDigAuxi::GetActiveDetectors() const {
    std::vector<int> ids;
    for (const auto& kv : fHitsByDet) {
        ids.push_back(kv.first);
    }
    return ids;
}
*/

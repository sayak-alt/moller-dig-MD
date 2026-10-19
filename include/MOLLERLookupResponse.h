#pragma once
#include <map>
#include <string>
#include <vector>
class MOLLERLookupResponse {
public:
 struct Hit { int det=0, pid=0; double t=0, xl=0, yl=0, zl=0, px=0, py=0, pz=0, k=0, edep=0; };
 void Load(const std::string& directory, int ring=0);
 // Negative means rejected or outside position coverage.
 double MeanPE(const Hit& hit) const;
private:
 struct Position { double hmin,hmax,vmin,vmax,mean,langau,rms,res,p0,p1,p2,p3,p4; };
 struct Energy { double energy,scale; };
 std::map<int,std::vector<Position>> positions;
 std::map<int,std::vector<Energy>> energies;
};

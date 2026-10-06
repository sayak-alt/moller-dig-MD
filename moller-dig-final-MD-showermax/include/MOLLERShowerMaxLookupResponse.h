#pragma once
#include <map>
#include <string>
#include <vector>
class MOLLERShowerMaxLookupResponse {
public:
 struct Hit {int det=0,pid=0; double t=0,e=0,x=0,y=0,pz=0;};
 void Load(const std::string& directory, bool legacy_pion_names=true);
 double MeanPE(const Hit& hit) const;
private:
 struct Fit {double energy,x0,x1,y0,y1,y2;};
 std::map<int,std::vector<Fit>> tables;
};

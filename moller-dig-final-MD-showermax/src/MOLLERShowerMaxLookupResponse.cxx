#include "MOLLERShowerMaxLookupResponse.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <array>
void MOLLERShowerMaxLookupResponse::Load(const std::string& directory,bool legacy) {
 std::map<int,std::vector<Fit>> loaded;
 const std::pair<int,std::string> names[]={ {11,"e-"},{-11,"e+"},{22,"gamma"},{13,"mu-"},{-13,"mu+"},{211,legacy?"pi-":"pi+"},{-211,legacy?"pi+":"pi-"},{2112,"neutron"} };
 for(const auto& name:names) {
  auto path=directory+"/fit_param_xy_"+name.second+"_ifarm.csv";
  std::ifstream in(path);if(!in) throw std::runtime_error("Cannot open shower-max LUT: "+path);
  std::string line;std::getline(in,line);
  if(!line.empty()&&line.back()=='\r')line.pop_back();
  if(line!="energy,xp0,xp1,yp0,yp1,yp2")throw std::runtime_error("Unexpected shower-max CSV header: "+path);
  std::map<double,std::array<double,6>> unique;
  while(std::getline(in,line)) {
   if(line.find_first_not_of(" \t\r")==std::string::npos || line[0]=='#')continue;
   std::stringstream ss(line);std::string token;std::vector<double> v;
   while(std::getline(ss,token,',')) {
    size_t used=0;double x=std::stod(token,&used);
    if(!std::isfinite(x)||token.find_first_not_of(" \t\r",used)!=std::string::npos)throw std::runtime_error("Invalid shower-max value: "+path);
    v.push_back(x);
   }
   if(v.size()!=6||v[0]<=0)throw std::runtime_error("Invalid shower-max row: "+path);
   std::array<double,6>a;std::copy(v.begin(),v.end(),a.begin());
   auto found=unique.find(v[0]);if(found!=unique.end()&&found->second!=a)throw std::runtime_error("Conflicting duplicate energy: "+path);
   unique[v[0]]=a;
  }
  if(unique.size()<2)throw std::runtime_error("Too few shower-max energies: "+path);
  for(const auto& row:unique){const auto& v=row.second;loaded[name.first].push_back({v[0],v[1],v[2],v[3],v[4],v[5]});}
 }
 tables.swap(loaded);
}
double MOLLERShowerMaxLookupResponse::MeanPE(const Hit& h) const {
 if(h.det<73001||h.det>73028||!std::isfinite(h.e)||!std::isfinite(h.x)||!std::isfinite(h.y)||!std::isfinite(h.pz)||!std::isfinite(h.t)||h.e<=10||h.pz<=0)return -1;
 double r=std::hypot(h.x,h.y);if(r<=1020||r>=1180)return -1;
 auto found=tables.find(h.pid);if(found==tables.end())return -1;
 const auto& rows=found->second;
 // Match the original helper's fixed energy grid and fallback bounds.
 const double grid[]={5,10,50,100,500,1000,2000,3000,4000,5000,6000,7000,8000,9000};
 double low=5,high=9000;
 for(int i=0;i<13;++i)if(h.e>=grid[i]&&h.e<grid[i+1]){low=grid[i];high=grid[i+1];break;}
 constexpr double pi=3.14159265358979323846;
 double phi=std::atan2(h.y,h.x);if(phi<0)phi+=2*pi;
 int nearest=0;double best=pi;
 for(int i=0;i<28;++i){double delta=std::abs(phi-i*2*pi/28);delta=std::min(delta,2*pi-delta);if(delta<best){best=delta;nearest=i;}}
 double angle=nearest*2*pi/28;
 double x=h.x*std::cos(angle)+h.y*std::sin(angle)-1100;
 double y=-h.x*std::sin(angle)+h.y*std::cos(angle);
 auto eval=[&](const Fit& f){return .5*(f.x0+x*f.x1)+.5*(f.y0+y*f.y1+y*y*f.y2);};
 // A missing energy row gives zero, matching the original default TF2 parameters.
 auto atEnergy=[&](double energy){
   for(const auto& row:rows)if(row.energy==energy)return eval(row);
   return 0.;
 };
 double vlo=atEnergy(low),vhi=atEnergy(high);
 double value=vlo+(vhi-vlo)*(h.e-low)/(high-low);
 double preserve=((h.det-73000)%4==3)?.22:.33;
 value*=preserve;
 return std::isfinite(value)?std::max(value,0.):-1;
}

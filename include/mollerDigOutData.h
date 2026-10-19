#ifndef MOLLERDigOutData_H
#define MOLLERDigOutData_H

#include <Rtypes.h>
class TTree;
class TBranch;

struct VDigOutData_t {
	VDigOutData_t() {};
	virtual ~VDigOutData_t() {};
	virtual bool SetupBranches(TTree *t, const char* prefix) = 0;
	template<typename T> int SetupBranch(TTree* tree, const char* prefix, const char* varname, T &var);
};

struct DigGEMData_t : public VDigOutData_t {
	Int_t nstrips;
	std::vector<Int_t> *module;
	std::vector<Int_t> *strip;
	std::vector<Int_t> *adc;
	std::vector<Int_t> *samp;

	TBranch *b_nstrips;
	TBranch *b_module;
	TBranch *b_strip;
	TBranch *b_adc;
	TBranch *b_samp;
	
	DigGEMData_t() : nstrips(0), module(0), strip(0), adc(0), samp(0), b_nstrips(0), b_module(0), b_strip(0), b_adc(0), b_samp(0) {};

	virtual ~DigGEMData_t() {};
	virtual bool SetupBranches(TTree *t, const char *prefix);
	void ClearBranches();
	void FillBranches();
};

#endif

#include "mollerDigOutData.h"
#include <TTree.h>
#include <TBranch.h>
#include <iostream>

template<typename T> int VDigOutData_t::SetupBranch(TTree* tree, const char* prefix, const char* varname, T &var){
	TString branchname = TString::Format("%s",varname);
	if(!tree)
	return 1;


	return 0;
}

bool DigGEMData_t::SetupBranches(TTree* tree, const char *prefix)
{
	if(!tree)return(false);
	
	module = new std::vector<int>;
	strip = new std::vector<int>;
	adc = new std::vector<int>;
	samp = new std::vector<int>;
	
	b_nstrips = tree->Branch(Form("%s.nstrips", prefix), &nstrips);
	b_module = tree->Branch(Form("%s.module", prefix), &module);
	b_strip = tree->Branch(Form("%s.strip", prefix), &strip);
	b_adc = tree->Branch(Form("%s.adc", prefix), &adc);
	b_samp = tree->Branch(Form("%s.samp", prefix), &samp);
	return true;
}
  
void DigGEMData_t::ClearBranches()
{
	if(strip){
	nstrips = 0;
	strip->clear();
	module->clear();
	adc->clear();
	samp->clear();
	}   
}
  
void DigGEMData_t::FillBranches()
{
/*	if(b_nstrips){
	b_nstrips->Fill();
	b_module->Fill();
	b_strip->Fill();
	b_adc->Fill();
	b_samp->Fill();
	}
*/
}

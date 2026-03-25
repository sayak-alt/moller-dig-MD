#include "remoll_tree.h"
#include <iostream>

// IMPORTANT: include remoll class
#include "remollGenericDetectorHit.hh"
//#include "mollerDigOutData.h"

void remoll_tree::Init(TTree* tree) {
    if (!tree) {
        std::cerr << "Error: null tree!" << std::endl;
        return;
    }

    fChain = tree;

    // Bind the full object branch
    fChain->SetBranchAddress("hit", &hit);

    std::cout << "Initialized tree with "
              << fChain->GetEntries() << " entries\n";

//    SetupDetBranch(mollergem_dig,"moller");

}

void remoll_tree::SetupOutputBranches(TTree* out_tree) {
    SetupDetBranch(mollergem_dig, "moller", out_tree);
}

Long64_t remoll_tree::GetEntry(Long64_t entry) {
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

void remoll_tree::ClearDigBranches()
{   
    mollergem_dig.ClearBranches();
}

void remoll_tree::FillDigBranches()
{   
    mollergem_dig.FillBranches();
}

void remoll_tree::SetupDetBranch(VDigOutData_t &det, const char* prefix, TTree* tree)
{
    TTree* target = tree ? tree : fChain;   
    det.SetupBranches(target,prefix);
}

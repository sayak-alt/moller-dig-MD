#ifndef REMOLL_TREE_H
#define REMOLL_TREE_H

#include "TTree.h"
#include <vector>
#include "mollerDigOutData.h"

// Forward declaration (comes from remoll)
class remollGenericDetectorHit_t;
class mollerDigOutData;

class remoll_tree {
public:
    remoll_tree() = default;
    ~remoll_tree() = default;

    DigGEMData_t mollergem_dig;
     
    void Init(TTree* tree);
    void SetupOutputBranches(TTree* out_tree);

    Long64_t GetEntry(Long64_t entry);
    Long64_t GetEntries() const { return fChain ? fChain->GetEntries() : 0; }

    // --- Data members ---
    std::vector<remollGenericDetectorHit_t>* hit = nullptr;

//private:
    TTree* fChain = nullptr;

    void ClearDigBranches();
    void FillDigBranches();

private:
    void SetupDetBranch(VDigOutData_t &det, const char* prefix, TTree* tree  =nullptr);
};

#endif

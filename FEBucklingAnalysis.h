#ifndef _FEBUCKLINGANALYSIS_
#define _FEBUCKLINGANALYSIS_

#include "FEAnalysis.h"

class FEBucklingAnalysis : public FEModeOperator
{
public:
    std::shared_ptr<FEDeformOperator> InitailDeformOp = nullptr;
    std::shared_ptr<FEDeformOperator> deform_case;
    FEBucklingAnalysis(std::shared_ptr<FEDeformOperator> deform_op) : deform_case(deform_op), FEModeOperator(deform_op->model) {};

    int ModeNum() override { return mode_num; }
    std::vector<std::vector<Displacement>> ModeVectors() override { return mode_vectors; }
    std::vector<double> EigenValues() override { return eigs; }

    int SolveBuckling();
    int mode_num = 1;
    std::vector<double> eigs;
    std::vector<std::vector<Displacement>> mode_vectors;
};

#endif
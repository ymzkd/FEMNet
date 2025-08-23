#ifndef _FEANALYSIS_
#define _FEANALYSIS_

#include "Model.h"
#include "Elements/Elements.h"

// 変形ケースを表す基底クラス
class FEDeformOperator
{
public:
    std::shared_ptr<FEModel> model;

    FEDeformOperator() {};
    FEDeformOperator(std::shared_ptr<FEModel> model) : model(model) {};

    virtual BeamStressData GetBeamStress(int eid, double p) = 0;
    virtual PlateStressData GetPlateStressData(int eid, double xi, double eta) = 0;
    virtual Displacement GetBeamDisplace(int eid, double p) = 0;

    virtual std::vector<Displacement> GetDisplacements() = 0;
    virtual std::vector<Displacement> GetVelocities(){return std::vector<Displacement>();};
    virtual std::vector<Displacement> GetAccelerations(){return std::vector<Displacement>();};

    virtual std::vector<NodeLoad> GetReactForces()
    {
        return std::vector<NodeLoad>();
    };
};

class FEModeOperator
{
public:
    std::shared_ptr<FEModel> model;
    FEModeOperator() {};
    FEModeOperator(std::shared_ptr<FEModel> model) : model(model) {};

    virtual std::vector<std::vector<Displacement>> ModeVectors() = 0;
    virtual std::vector<double> EigenValues() = 0;
    virtual int ModeNum() = 0;
};

#endif
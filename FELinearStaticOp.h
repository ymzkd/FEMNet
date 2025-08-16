#ifndef _FELINEARSTATICOP_
#define _FELINEARSTATICOP_

#include "FEAnalysis.h"

// Static Solver Package
class FELinearStaticOp : public FEDeformOperator
{
private:
    bool m_computed = false;

public:
    std::vector<std::shared_ptr<LoadBase>> loads;
    std::vector<Displacement> displace;
    std::vector<NodeLoad> react_force;

    FELinearStaticOp(
        std::shared_ptr<FEModel> model,
        std::vector<std::shared_ptr<LoadBase>> loads)
        : FEDeformOperator(model), loads(loads) {
          };

    bool Computed() { return m_computed; }
    void Compute();

    BeamStressData GetBeamStress(int eid, double p) override;

    /// <summary>
    /// Obtain stress data for plate elements
    /// </summary>
    /// <param name="eid">element index</param>
    /// <param name="xi">
    /// xi for square isoparametric elements and L1 for triangular element
    /// area coordinate system in the first parameter.
    /// </param>
    /// <param name="eta">
    /// eta for square isoparametric elements and L2 for triangular element
    /// area coordinate system in the second parameter.
    /// </param>
    /// <returns>Plate element stress data</returns>
    PlateStressData GetPlateStressData(int eid, double xi, double eta) override;

    Displacement GetBeamDisplace(int eid, double p) override;

    // FEDeformCase を介して継承されました
    std::vector<Displacement> GetDisplacements() override;

    std::vector<NodeLoad> GetReactForces() override;
}; // FEStaticResult

struct LinearStaticDeformFactor
{
public:
    std::shared_ptr<FELinearStaticOp> op;
    double factor;
    LinearStaticDeformFactor(std::shared_ptr<FELinearStaticOp> op, double factor)
        : op(op), factor(factor) {
          };
};

class LinearStaticCombinationOperator : public FEDeformOperator
{
public:
    std::vector<LinearStaticDeformFactor> cases;
    LinearStaticCombinationOperator() {};
    LinearStaticCombinationOperator(std::vector<LinearStaticDeformFactor> cases)
        : cases(cases) {
          };

    // FEDeformCase を介して継承されました
    BeamStressData GetBeamStress(int eid, double p) override;

    PlateStressData GetPlateStressData(int eid, double xi, double eta) override;

    Displacement GetBeamDisplace(int eid, double p) override;

    std::vector<Displacement> GetDisplacements() override;

    std::vector<NodeLoad> GetReactForces() override;

}; // StaticCombinationOperator

#endif
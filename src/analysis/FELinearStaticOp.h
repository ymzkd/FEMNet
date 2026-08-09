#ifndef _FELINEARSTATICOP_
#define _FELINEARSTATICOP_

#include "FEAnalysis.h"

struct LinearStaticDeformFactor;

// Static Solver Package
class FELinearStaticOp : public FEDeformOperator
{
private:
    bool m_computed = false;

    // 状態依存要素(非抗圧トラス等)の解析状態。要素内の可変状態ではなく
    // Operatorが唯一の保存場所として所有する。Compute()の反復で更新され、
    // 別ケースの解析後もこのケースの応力を正しく復元できる。
    ElementStates m_states;

public:
    std::vector<std::shared_ptr<LoadBase>> loads;
    std::vector<Displacement> displace;
    std::vector<NodeLoad> react_force;

    // 直近のCompute()における状態依存要素の反復回数
    int Iterations = 0;
    // 状態依存要素の反復が収束したか(最大反復数到達でfalse)
    bool Converged = true;

    FELinearStaticOp(
        std::shared_ptr<FEModel> model,
        std::vector<std::shared_ptr<LoadBase>> loads)
        : FEDeformOperator(model), loads(loads) {
          };

    // 荷重組み合わせ(ケースOp+係数のリスト)を、係数倍した荷重を連結した
    // 単一の静的荷重ケースとして構築する。
    // 状態依存要素(引張専用トラス等)を含むモデルでは、結果の重ね合わせ
    // (LinearStaticCombinationOperator)の代わりにこちらを使用する。
    static std::shared_ptr<FELinearStaticOp> FromCombination(
        std::shared_ptr<FEModel> model,
        const std::vector<LinearStaticDeformFactor>& cases);

    bool Computed() { return m_computed; }

    // 状態依存要素の反復解法による線形静的解析。
    // 要素状態(m_states)はOperatorが所有し、組立時にFEModelへ明示的に渡す。
    void Compute();

    // 状態依存要素の収束時状態を返す(未登録・非状態依存要素はtrue)
    bool GetElementState(int eid) const { return m_states.Get(eid); }

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

    std::vector<NodeLoadData> GetPlateNodalForces(int eid, bool local) override;

    Displacement GetBeamDisplace(int eid, double p) override;

    // FEDeformCase を介して継承されました
    std::vector<Displacement> GetDisplacements() override;

    std::vector<NodeLoad> GetReactForces() override;
}; // FELinearStaticOp

struct LinearStaticDeformFactor
{
public:
    std::shared_ptr<FELinearStaticOp> op;
    double factor;
    LinearStaticDeformFactor() : op(nullptr), factor(0.0) {};
    LinearStaticDeformFactor(std::shared_ptr<FELinearStaticOp> op, double factor)
        : op(op), factor(factor) {
          };
};

class LinearStaticCombinationOperator : public FEDeformOperator
{
public:
    std::vector<LinearStaticDeformFactor> cases;
    LinearStaticCombinationOperator() {};

    LinearStaticCombinationOperator(
        std::shared_ptr<FEModel> model,
        std::vector<LinearStaticDeformFactor> cases)
        : FEDeformOperator(model), cases(cases) {
          };

    // FEDeformCase を介して継承されました
    BeamStressData GetBeamStress(int eid, double p) override;

    // 各LinearStaticDeformFactorのインデックスに対応する応力を返す
    std::vector<BeamStressData> GetBeamStressComponents(int eid, double p);

    PlateStressData GetPlateStressData(int eid, double xi, double eta) override;

    // 各LinearStaticDeformFactorのインデックスに対応する応力を返す
    std::vector<PlateStressData> GetPlateStressDataComponents(int eid, double xi, double eta);

    std::vector<NodeLoadData> GetPlateNodalForces(int eid, bool local) override;

    Displacement GetBeamDisplace(int eid, double p) override;

    std::vector<Displacement> GetDisplacements() override;

    std::vector<NodeLoad> GetReactForces() override;

}; // StaticCombinationOperator

#endif

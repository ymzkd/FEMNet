#ifndef _FEVIBRATIONANALYSIS_
#define _FEVIBRATIONANALYSIS_

#include "FEAnalysis.h"

// 固有値(振動モード)解析Operator。
// Compute(nev) で一般化固有値問題 K φ = ω^2 M φ を解き、
// 固有値(角振動数)とモード形状を保持する。
// 解析済みの結果を直接与えて結果コンテナとして使うこともできる。
class FEVibrationAnalysis : public FEModeOperator
{
private:
    std::vector<double> eigs; // Omegas
    std::vector<std::vector<Displacement>> mode_vectors;
    bool m_computed = false;

public:
    int ModeNum() override { return eigs.size(); };
    const std::vector<std::vector<Displacement>>& ModeVectors() override { return mode_vectors; };
    std::vector<double> EigenValues() const override { return eigs; };

    FEVibrationAnalysis() {};

    FEVibrationAnalysis(std::shared_ptr<FEModel> model) : FEModeOperator(model) {};

    // 解析済みの結果から構築する(結果コンテナとしての利用)
    FEVibrationAnalysis(
        std::shared_ptr<FEModel> model,
        std::vector<std::vector<Displacement>> mode_vectors,
        std::vector<double> eigs)
        : mode_vectors(mode_vectors),
          eigs(eigs), FEModeOperator(model) { m_computed = true; };

    bool Computed() { return m_computed; }

    // 固有値解析を実行する。
    // nev: 要求モード数
    // 戻り値: 収束したモード数(失敗時は負値)
    int Compute(const int nev);

    // モードごとの刺激係数計算
    std::vector<double> ParticipationFactors();
    std::vector<double> ParticipationFactors(Vector direction);

    // 方向別のモードごとの刺激係数計算
    std::vector<Displacement> ParticipationDirectedFactors();

    // モードごとの有効質量比計算
    std::vector<double> EffectiveMassRates();

    // 方向別のモードごとの有効質量比計算
    std::vector<Displacement> EffectiveDirectedMassRates();

    std::vector<double> NaturalPeriods()
    {
        std::vector<double> periods;

        for (double eig : eigs)
            periods.push_back(2.0 * PI / eig);

        return periods;
    }
};

#endif

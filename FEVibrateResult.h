#ifndef _FEVIBRATERESULT_
#define _FEVIBRATERESULT_

#include "FEAnalysis.h"

class FEVibrateResult : public FEModeOperator
{
private:
    std::vector<double> eigs; // Omegas
    std::vector<std::vector<Displacement>> mode_vectors;

public:
    int ModeNum() override { return eigs.size(); };
    std::vector<std::vector<Displacement>> ModeVectors() override { return mode_vectors; };
    std::vector<double> EigenValues() override { return eigs; };

    FEVibrateResult() {};

    FEVibrateResult(
        std::shared_ptr<FEModel> model,
        std::vector<std::vector<Displacement>> mode_vectors,
        std::vector<double> eigs)
        : mode_vectors(mode_vectors),
          eigs(eigs), FEModeOperator(model) {};

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
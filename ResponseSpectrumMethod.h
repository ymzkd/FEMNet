#ifndef _RESPONSE_SPECTRUM_METHOD_
#define _RESPONSE_SPECTRUM_METHOD_

#include "FEAnalysis.h"
#include "FELinearStaticOp.h"
#include "FEVibrateResult.h"


enum ResponseSpectrumMethodType
{
    ABS,
    SRSS,
    CQC
};

enum class ResponseValueType
{
    Displacement,
    Velocity,
    Acceleration
};

class ResponseSpectrumMethod : public FEDeformOperator
{
private:
    std::vector<Displacement> calculate_responseCQC(ResponseValueType vt);
    std::vector<Displacement> calculate_responseSRSS(ResponseValueType vt);
    std::vector<Displacement> calculate_responseABS(ResponseValueType vt);
    std::vector<Displacement> calculate_response(ResponseValueType vt);

    std::vector<Displacement> displacements; // 解析結果の変位ベクトル
    std::vector<Displacement> velocities;    // 解析結果の速度ベクトル
    std::vector<Displacement> accelerations; // 解析結果の加速度ベクトル

    bool m_computed = false;

public:
    double damping_rate = 0.05;               // 減衰比(CQC法の場合のみ計算に影響)
    Vector Direction = Vector(1.0, 1.0, 1.0); // 応答スペクトルの方向

    FEVibrateResult VibrateResult;
    IResponseSpectrum *SpectrumFunction;
    ResponseSpectrumMethodType MethodType = ResponseSpectrumMethodType::ABS;

    ResponseSpectrumMethod() {}
    ResponseSpectrumMethod(std::shared_ptr<FEModel> model, FEVibrateResult vibrate_result, Vector direction,
                           IResponseSpectrum *spectrum_function, ResponseSpectrumMethodType type);

    bool Computed() { return m_computed; }
    void Compute();
    std::vector<Displacement> GetDisplacements() override;
    std::vector<Displacement> GetVelocities();
    std::vector<Displacement> GetAccelerations();

    // FEDeformCase を介して継承されました
    BeamStressData GetBeamStress(int eid, double p) override;
    PlateStressData GetPlateStressData(int eid, double xi, double eta) override;
    Displacement GetBeamDisplace(int eid, double p) override;

    FELinearStaticOp GetLinearStaticCase();
};

#endif
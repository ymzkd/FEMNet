#ifndef _RESPONSE_SPECTRUM_METHOD_
#define _RESPONSE_SPECTRUM_METHOD_

#ifndef SWIG
#include <cmath>
#endif

#include "FEAnalysis.h"
#include "FELinearStaticOp.h"
#include "FEVibrationAnalysis.h"


// 応答スペクトルのモード合成方式。
//   ABS     : 各モード応答の絶対値和
//   SRSS    : 二乗和平方根
//   CQC     : 完全二次結合(Complete Quadratic Combination)
//   AWA_ABS : AWA法(最悪モード W + 残差項)。残差 = calculate_responseABS (絶対値和 0.5·Σ_j|β_j S_j φ_j|)
//   AWA_CQC : AWA法(最悪モード W + 残差項)。残差 = calculate_responseCQC (0.5·CQC結合場)
// AWA_ABS/AWA_CQC はいずれも残差に AWA の平均化 0.5 を掛け、合成則(成分ごとローカル符号
// reinforce: R = W + sign(W)·|残差|)は共通。違いは残差に使う関数のみ。
enum ResponseSpectrumMethodType
{
    ABS,
    SRSS,
    CQC,
    AWA_ABS,
    AWA_CQC,
};

// 応答成分の符号調整方式。CQC/SRSS/ABS や剛応答のSRSS合成は符号を失い正値になるため、
// 各成分の絶対値は保持したまま符号を指定モードの形状の符号に合わせる。
//   SIGN_NONE                : 符号調整なし(そのまま。CQC等は正値, AWAは自身の符号)
//   SIGN_SPECIFIED_MODE      : 指定次数(sign_mode_index, 0始まり)のモード形状の符号
//   SIGN_WORST_MODE          : AWAの最悪モード(|β|/ω 最大)の符号(採用次数は sign_mode_index に記録)
//   SIGN_STRAIN_ENERGY_MODE  : モードひずみエネルギー sE∝β²·Sv² 最大のモードの符号(採用次数を記録)
enum ResponseSignType
{
    SIGN_NONE,
    SIGN_SPECIFIED_MODE,
    SIGN_WORST_MODE,
    SIGN_STRAIN_ENERGY_MODE,
};

enum class ResponseValueType
{
    Displacement,
    Velocity,
    Acceleration
};

struct RigidResponseComposition{
public:
    double f1, f2;
    RigidResponseComposition() : f1(0.0), f2(0.0){}
    RigidResponseComposition(double f1_, double f2_) : f1(f1_), f2(f2_){}

    // ti: 固有周期(s)。式は周波数ベースのため fi=1/ti に変換して評価する。
    double RigidResponseFactor(const double ti){
        double fi = 1.0 / ti;
        if (fi <= f1) return 0.0;
        if (fi >= f2) return 1.0;
        return std::log(fi / f1) / std::log(f2 / f1);
    }
};

class ResponseSpectrumMethod : public FEDeformOperator
{
private:
    std::vector<Displacement> calculate_responseAWA(ResponseValueType vt);
    std::vector<Displacement> calculate_responseCQC(ResponseValueType vt);
    std::vector<Displacement> calculate_responseSRSS(ResponseValueType vt);
    std::vector<Displacement> calculate_responseABS(ResponseValueType vt);
    std::vector<Displacement> calculate_rigid_response(ResponseValueType vt);

    size_t worst_mode_index();              // AWAの最悪モード(エネルギー寄与|β|/ω が最大)の添字
    size_t max_strain_energy_mode_index();  // モードひずみエネルギー sE∝β²·Sv² が最大の添字(符号採用用)

    std::vector<Displacement> calculate_response(ResponseValueType vt);
    std::vector<NodeLoad> calculate_react_forces(const std::vector<Displacement> &disp);

    std::vector<Displacement> displacements; // 解析結果の変位ベクトル
    std::vector<Displacement> velocities;    // 解析結果の速度ベクトル
    std::vector<Displacement> accelerations; // 解析結果の加速度ベクトル
    std::vector<NodeLoad> react_forces;      // 合成変位分布に基づく支点反力

    bool m_computed = false;

public:
    // CQCの減衰比はSpectrumFunction->effective_damping_rate()から取得する
    // (旧damping_rateは廃止。スペクトル側と二重管理になっていたため)
    Vector Direction = Vector(1.0, 1.0, 1.0); // 応答スペクトルの方向

    FEVibrationAnalysis VibrateResult;
    IResponseSpectrum *SpectrumFunction;
    ResponseSpectrumMethodType MethodType = ResponseSpectrumMethodType::ABS;
    ResponseSignType sign_type = ResponseSignType::SIGN_NONE; // 応答成分の符号調整(既定=なし)。Compute()前に設定する。
    int sign_mode_index = -1; // 符号付与に用いるモード添字(0始まり)。SIGN_SPECIFIED_MODEでは入力、
                              // SIGN_WORST_MODEではCompute時に採用次数を記録。既定-1(モード非依存)。

    // 剛応答・欠落質量の考慮オプション
    bool EnableRigidResponse = false;
    RigidResponseComposition RigidResponse; // 剛応答と欠落質量補正

    ResponseSpectrumMethod() {}
    ResponseSpectrumMethod(std::shared_ptr<FEModel> model, FEVibrationAnalysis vibrate_result, Vector direction,
                           IResponseSpectrum *spectrum_function, ResponseSpectrumMethodType type);

    bool Computed() { return m_computed; }
    void Compute();
    std::vector<Displacement> GetDisplacements() override;
    std::vector<Displacement> GetVelocities() override;
    std::vector<Displacement> GetAccelerations() override;

    // 合成後の変位分布から K・u で算出した反力を返す
    std::vector<NodeLoad> GetReactForces() override;

    // FEDeformCase を介して継承されました
    BeamStressData GetBeamStress(int eid, double p) override;
    PlateStressData GetPlateStressData(int eid, double xi, double eta) override;
    Displacement GetBeamDisplace(int eid, double p) override;

    FELinearStaticOp GetLinearStaticCase();
};

#endif
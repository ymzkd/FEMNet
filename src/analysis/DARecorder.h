#ifndef _DARECORDER_
#define _DARECORDER_

#ifndef SWIG
#include <vector>
#include <string>
#endif

// 前方宣言
class DynamicAnalysis;

/// <summary>
/// 時刻歴応答解析の各ステップで任意の量を記録するレコーダの基底クラス。
/// DASampler と同様に DynamicAnalysis へ複数登録でき、C#/Python 側でも派生できる。
/// </summary>
class DARecorder
{
public:
    virtual ~DARecorder() = default;

    std::string Name;
    std::string Description;
    std::vector<double> Values;

    DARecorder() = default;
    DARecorder(std::string name, std::string description = "")
        : Name(name), Description(description) {}

    /// 解析の初期化時(DynamicAnalysis::Initialize)に呼ばれる。記録バッファの初期化に使う。
    virtual void Initialize(DynamicAnalysis &da) { Values.clear(); }
    /// 各ステップの計算後に呼ばれる。
    virtual void Record(DynamicAnalysis &da) = 0;
};

/// <summary>
/// 運動エネルギー 1/2·vᵀ·M·v を記録する。
/// </summary>
class DARecorder_KineticEnergy : public DARecorder
{
public:
    DARecorder_KineticEnergy()
        : DARecorder("KineticEnergy", "Kinetic energy 1/2*v^T*M*v") {}

    /// step 0 の値として 0 を積んでから記録を開始する
    void Initialize(DynamicAnalysis &da) override { DARecorder::Initialize(da); Values.push_back(0.0); }
    void Record(DynamicAnalysis &da) override;
};

/// <summary>
/// ポテンシャル(ひずみ)エネルギー 1/2·dᵀ·K·d を記録する。
/// </summary>
class DARecorder_PotentialEnergy : public DARecorder
{
public:
    DARecorder_PotentialEnergy()
        : DARecorder("PotentialEnergy", "Potential (strain) energy 1/2*d^T*K*d") {}

    void Initialize(DynamicAnalysis &da) override { DARecorder::Initialize(da); Values.push_back(0.0); }
    void Record(DynamicAnalysis &da) override;
};

/// <summary>
/// 減衰による消散量 vᵀ·C·v を記録する(各ステップの瞬時値。時間積分は利用側で行う)。
/// </summary>
class DARecorder_DampingEnergy : public DARecorder
{
public:
    DARecorder_DampingEnergy()
        : DARecorder("DampingEnergy", "Damping dissipation v^T*C*v") {}

    void Initialize(DynamicAnalysis &da) override { DARecorder::Initialize(da); Values.push_back(0.0); }
    void Record(DynamicAnalysis &da) override;
};

/// <summary>
/// 外力による入力量 -F_ext·v を記録する(各ステップの瞬時値。時間積分は利用側で行う)。
/// </summary>
class DARecorder_InputEnergy : public DARecorder
{
public:
    DARecorder_InputEnergy()
        : DARecorder("InputEnergy", "Input energy -F_ext*v") {}

    void Initialize(DynamicAnalysis &da) override { DARecorder::Initialize(da); Values.push_back(0.0); }
    void Record(DynamicAnalysis &da) override;
};

#endif

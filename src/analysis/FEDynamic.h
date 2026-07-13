#ifndef _FEDYNAMIC_
#define _FEDYNAMIC_

#ifndef SWIG
#include <vector>
#include <map>
#include <array>

#include "SparseSolver.h"

#endif

#include "FEAnalysis.h"
#include "LoadComponent.h"

// 前方宣言
class DynamicAnalysis;

class DASampler
{

public:
    int step;
    std::string Name;
    std::vector<Displacement> velocity, displacement, acceleration;

    DASampler() : step(0), Name("") {}
    DASampler(std::string name) : step(0), Name(name) {}

    virtual void Sampling(DynamicAnalysis &analysis) = 0;
};

class DASampler_MaxDisplacement : public DASampler
{
public:
    double max_displacement = 0.0;

    DASampler_MaxDisplacement() : DASampler("MaxDisplacement") {}

    void Sampling(DynamicAnalysis &da) override;
};

class DARecorder
{
public:
    virtual void Record(DynamicAnalysis &da) = 0;
};

class DAEnergyRecorder : public DARecorder
{
public:
    std::vector<double> kinetic_energy, potential_energy, damping_energy, input_energy;

    void Initialize();
    void RecordKineticEnergy(DynamicAnalysis &da);
    void RecordPotentialEnergy(DynamicAnalysis &da);
    void RecordDampingEnergy(DynamicAnalysis &da);
    void RecordInputEnergy(DynamicAnalysis &da);

    void Record(DynamicAnalysis &da) override;
};

/// <summary>
/// 時刻歴応答解析クラス
/// </summary>
class DynamicAnalysis : public FEDeformOperator
{
private:
    std::unique_ptr<ISparseSolver> solver;
    Eigen::SparseMatrix<double> matM_aa, matM_ab, matM_bb;
    Eigen::SparseMatrix<double> matK_aa, matK_ab, matK_bb;
    Eigen::SparseMatrix<double> matC_aa, matC_ab, matC_bb;
    std::vector<int> free_indices, fixed_indices;
    std::vector<int> slave_indices;           // RigidLink: スレーブDOFインデックス
    Eigen::SparseMatrix<double> linkTransMat; // RigidLink: 変換行列
    int master_dof_num = 0;                   // RigidLink: マスターDOF数

    Eigen::VectorXd current_disp, current_vel, current_accel;
    std::vector<NodeLoad> current_react_force;

    friend class FEDynamicDampInitializer;
    friend class FEDynamicStiffDampInitializer;
    friend class FEDynamicMassDampInitializer;
    friend class FEDynamicRayleighDampInitializer;
    friend class DAEnergyRecorder;

public:
    DynamicAccelLoad accel_load;
    FEDynamicDampInitializer *damp_initializer = nullptr;

    int current_step = 0;

    std::vector<std::shared_ptr<DASampler>> samplers;
    DAEnergyRecorder energy_recorder;
    bool RecordEnabled = true;

    double beta = 0.25; // 平均加速度法
                        // double beta = 1.0/6.0; // 線形加速度法(発散しがち)

    DynamicAnalysis(std::shared_ptr<FEModel> model,
                    const DynamicAccelLoad& accel_load, FEDynamicDampInitializer *damp = nullptr);

    bool Initialize();

    void Clear()
    {
        current_step = 0;
        solver.reset();
        matM_aa.resize(0, 0);
        matK_aa.resize(0, 0);
        matC_aa.resize(0, 0);
        linkTransMat.resize(0, 0);
        master_dof_num = 0;
        slave_indices.clear();
    }

    // Newmarkのβ法による動的解析
    void ComputeStep();
    void ComputeSteps(int steps);

    bool SetDisplacements(std::vector<Displacement> disps);
    bool SetVelocities(std::vector<Displacement> vels);
    bool SetAccelerations(std::vector<Displacement> accs);

    std::vector<Displacement> GetDisplacements() override;
    std::vector<Displacement> GetVelocities() override;
    std::vector<Displacement> GetAccelerations() override;

    std::vector<BarElementBase *> GetBarElements()
    {
        std::vector<BarElementBase *> bar_elements;
        for (const auto &elem : model->Elements)
        {
            if (IsBarType(elem->Type()))
            {
                if (auto *bar_elem = dynamic_cast<BarElementBase *>(elem.get()))
                {
                    bar_elements.push_back(bar_elem);
                }
            }
        }
        return bar_elements;
    }

    std::vector<PlaneElementBase *> GetPlaneElements()
    {
        std::vector<PlaneElementBase *> plane_elements;
        for (const auto &elem : model->Elements)
        {
            if (IsPlaneType(elem->Type()))
            {
                if (auto *plane_elem = dynamic_cast<PlaneElementBase *>(elem.get()))
                {
                    plane_elements.push_back(plane_elem);
                }
            }
        }
        return plane_elements;
    }

    // FEDeformCase を介して継承されました
    BeamStressData GetBeamStress(int eid, double p) override;
    PlateStressData GetPlateStressData(int eid, double xi, double eta) override;
    Displacement GetBeamDisplace(int eid, double p) override;

    std::vector<NodeLoad> GetReactForces() override { return current_react_force; };
};

/// <summary>
/// 時刻歴応答解析における比例減衰マトリクスを初期化する抽象基底クラス
/// </summary>
class FEDynamicDampInitializer
{
private:
    DynamicAnalysis *analysis;

public:
    virtual bool Initialize(DynamicAnalysis *analysis) = 0;

    /// <summary>
    /// 固有周期 t における減衰比を返す(算定不能な場合は -1)
    /// </summary>
    virtual double DampRateAtPeriod(double t) = 0;
};

/// <summary>
/// 時刻歴応答解析において比例減衰マトリクスを剛性比例として初期化するクラス
/// </summary>
class FEDynamicStiffDampInitializer : public FEDynamicDampInitializer
{

public:
    double natural_angle_velocity = 0.0; // 自然角速度
    double damp_rate = 0.05;             // 減衰比

    FEDynamicStiffDampInitializer(double damp_rate = 0.05)
        : damp_rate(damp_rate) {}

    bool Initialize(DynamicAnalysis *analysis) override;
    double DampRateAtPeriod(double t) override;
};

/// <summary>
/// 時刻歴応答解析において比例減衰マトリクスを質量比例として初期化するクラス
/// </summary>
class FEDynamicMassDampInitializer : public FEDynamicDampInitializer
{

public:
    double natural_angle_velocity = 0.0; // 自然角速度
    double damp_rate = 0.05;             // 減衰比

    FEDynamicMassDampInitializer(double damp_rate = 0.05)
        : damp_rate(damp_rate) {}

    bool Initialize(DynamicAnalysis *analysis) override;
    double DampRateAtPeriod(double t) override;
};

/// <summary>
/// 時刻歴応答解析において比例減衰マトリクスをレイリー減衰(C = αM + βK)として初期化するクラス
/// </summary>
class FEDynamicRayleighDampInitializer : public FEDynamicDampInitializer
{

public:
    double alpha = 0.0;                          // 質量比例係数
    double beta = 0.0;                           // 剛性比例係数
    double damp_rate1 = 0.05, damp_rate2 = 0.05; // 対象モードの減衰比
    int mode1 = 1, mode2 = 2;                    // 対象モード次数(1始まり)
    double natural_angle_velocity1 = 0.0, natural_angle_velocity2 = 0.0; // 対象モードの自然角速度
    bool direct_coefficients = false;            // α, βを直接指定する場合true

    // α, βを直接指定
    FEDynamicRayleighDampInitializer(double alpha, double beta)
        : alpha(alpha), beta(beta), direct_coefficients(true) {}

    // 2つのモードの減衰比を指定(Initialize時の固有値解析でα, βを算出)
    FEDynamicRayleighDampInitializer(double damp_rate1, double damp_rate2, int mode1, int mode2)
        : damp_rate1(damp_rate1), damp_rate2(damp_rate2), mode1(mode1), mode2(mode2) {}

    bool Initialize(DynamicAnalysis *analysis) override;
    double DampRateAtPeriod(double t) override;
};
#endif
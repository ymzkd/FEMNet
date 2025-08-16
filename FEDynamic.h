#ifndef _FEDYNAMIC_
#define _FEDYNAMIC_

#ifndef SWIGCSHARP
#include <vector>
#include <map>
#include <array>

// #include <Eigen/Dense>

#ifdef EIGEN_USE_MKL_ALL
// #define EIGEN_USE_MKL_ALL
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>
#else
#include <Eigen/Sparse>
#endif

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
#ifdef EIGEN_USE_MKL_ALL
    Eigen::PardisoLLT<Eigen::SparseMatrix<double>> solver;
#else
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
#endif
    Eigen::SparseMatrix<double> matM_aa, matM_ab, matM_bb;
    Eigen::SparseMatrix<double> matK_aa, matK_ab, matK_bb;
    Eigen::SparseMatrix<double> matC_aa, matC_ab, matC_bb;
    std::vector<int> free_indices, fixed_indices;

    Eigen::VectorXd current_disp, current_vel, current_accel;
    std::vector<NodeLoad> current_react_force;

    friend class FEDynamicDampInitializer;
    friend class FEDynamicStiffDampInitializer;
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
                    DynamicAccelLoad accel_load, FEDynamicDampInitializer *damp = nullptr);

    bool Initialize();

    void Clear()
    {
        current_step = 0;
        matM_aa.resize(0, 0);
        matK_aa.resize(0, 0);
        matC_aa.resize(0, 0);
    }

    // Newmarkのβ法による動的解析
    void ComputeStep();
    void ComputeSteps(int steps);

    bool SetDisplacements(std::vector<Displacement> disps);
    bool SetVelocities(std::vector<Displacement> vels);
    bool SetAccelerations(std::vector<Displacement> accs);

    std::vector<Displacement> GetDisplacements();
    std::vector<Displacement> GetVelocities();
    std::vector<Displacement> GetAccelerations();

    std::vector<BarElementBase *> GetBarElements()
    {
        std::vector<BarElementBase *> bar_elements;
        for (const auto &elem : model->Elements)
        {
            if (elem->Type() == ElementType::Beam || elem->Type() == ElementType::Truss)
            {
                if (auto *bar_elem = dynamic_cast<BarElementBase *>(elem.get()))
                {
                    bar_elements.push_back(bar_elem);
                }
            }
        }
        return bar_elements;
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
};
#endif
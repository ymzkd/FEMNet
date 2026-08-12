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
#include "FEVibrationAnalysis.h"

// 前方宣言
class DynamicAnalysis;

class DASampler
{

public:
    virtual ~DASampler() = default;

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

class DASampler_MaxVelocity : public DASampler
{
public:
    double max_velocity = 0.0;

    DASampler_MaxVelocity() : DASampler("MaxVelocity") {}

    void Sampling(DynamicAnalysis &da) override;
};

class DASampler_MaxAcceleration : public DASampler
{
public:
    double max_acceleration = 0.0;

    DASampler_MaxAcceleration() : DASampler("MaxAcceleration") {}

    void Sampling(DynamicAnalysis &da) override;
};

class DASampler_MaxDispDirection : public DASampler
{
public:

    double max_displacement = 0.0;
    Vector direction;

    DASampler_MaxDispDirection() : DASampler("MaxDispDirection") {}
    DASampler_MaxDispDirection(Vector direction) : DASampler("MaxDispDirection")
    {
        this->direction = direction;
    }
    
    void Sampling(DynamicAnalysis &da) override;
};

class DASampler_MaxVelocityDirection : public DASampler
{
public:

    double max_velocity = 0.0;
    Vector direction;

    DASampler_MaxVelocityDirection() : DASampler("MaxVelocityDirection"){}
    DASampler_MaxVelocityDirection(Vector direction) : DASampler("MaxVelocityDirection")
    {
        this->direction = direction;
    }

    void Sampling(DynamicAnalysis &da) override;
};

class DASampler_MaxAccelDirection : public DASampler
{
public:
    double max_accel = 0.0;
    Vector direction;

    DASampler_MaxAccelDirection() : DASampler("MaxAccelDirection"){}
    DASampler_MaxAccelDirection(Vector direction) : DASampler("MaxAccelDirection")
    {
        this->direction = direction;
    }

    void Sampling(DynamicAnalysis &da) override;
};

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

/// <summary>
/// 等間隔サンプリングされた離散時系列データ。時刻 t に対する値を線形補間で返す。
/// データ区間 [t0, t0+(N-1)*dt] の外では 0 を返す(荷重の開始前・終了後は無載荷)。
/// </summary>
class TimeSeries
{
public:
    double t0 = 0.0;             // 先頭データ点の時刻
    double dt = 0.0;             // データの時間刻み
    std::vector<double> values;  // データ列

    TimeSeries() = default;
    TimeSeries(double dt, const std::vector<double>& values, double t0 = 0.0)
        : t0(t0), dt(dt), values(values) {}

    bool IsEmpty() const { return values.empty() || dt <= 0.0; }
    /// データ数
    int Count() const { return static_cast<int>(values.size()); }
    /// データの継続時間((N-1)*dt。データ点が1個以下なら0)
    double Duration() const { return IsEmpty() ? 0.0 : (values.size() - 1) * dt; }
    /// 末尾データ点の時刻
    double EndTime() const { return t0 + Duration(); }

    /// 時刻 t における値(線形補間)。データ区間外は 0。
    double Value(double t) const;
};

/// <summary>
/// 時刻歴応答解析における時刻歴荷重の基底クラス。
/// 時刻 t の関数として「全体節点空間」の荷重ベクトルを返す。
/// 縮約(RigidLink変換・固定DOF処理)は DynamicAnalysis 側で一括して行うため、
/// 派生クラスは全節点空間の NodeLoadData(固定DOF成分も含めてよい)を返せばよい。
///
/// 解析の時間刻み・ステップ数は DynamicAnalysis が所有する。荷重側は
/// suggested_timestep() / end_time() で「時間グリッド未指定時のヒント」を返すのみで、
/// 解析時間を拘束しない(ヒント不要な荷重は既定実装のまま 0 を返せばよい)。
///
/// 地震波のように離散データを補間して評価する荷重は has_time_series() が true を返し、
/// time_series() で元データ(データ数・時間刻み・開始時刻・値列)を公開する。波形の
/// プロットなど、解析グリッドではなく荷重固有のグリッドが必要な用途で利用する。
/// 時間グリッドのヒントと代表値は既定でこの時系列から導かれるため、離散データを持つ
/// 派生クラスは has_time_series() / time_series() を実装するだけでよい。
/// </summary>
class DynamicLoad
{
public:
    virtual ~DynamicLoad() = default;

    std::string Name;      // 識別用の名称(任意)
    double Factor = 1.0;   // 荷重倍率

    /// 時刻 t における全体節点荷重ベクトル
    virtual std::vector<NodeLoadData> load_vector(DynamicAnalysis& analysis, double t) = 0;

    /// 離散データ(時系列)を元に評価される荷重かどうか
    virtual bool has_time_series() const { return false; }
    /// 荷重の元データ時系列。has_time_series() が false の場合は空の時系列を返す。
    virtual const TimeSeries& time_series() const;

    /// 時間グリッド自動決定用のヒント: 推奨時間刻み(0 = ヒントなし)
    virtual double suggested_timestep() const { return has_time_series() ? time_series().dt : 0.0; }
    /// 時間グリッド自動決定用のヒント: 荷重の終端時刻(0 = ヒントなし・定常荷重)
    virtual double end_time() const { return has_time_series() ? time_series().EndTime() : 0.0; }

    /// 入力波形の表示等に用いる時刻 t の代表スカラー値
    /// (地動加速度・時刻係数など。荷重倍率は含まない)。既定は時系列の補間値。
    virtual double reference_value(double t) const
    {
        return has_time_series() ? time_series().Value(t) : 0.0;
    }
};

/// <summary>
/// 地震(地動加速度)による慣性外力 -M·ι·a_g(t) を表す時刻歴荷重。
/// 質量は集中質量(対角)のため、各節点の並進成分に -(m/g)·Direction·a_g を
/// 与える(回転成分は0)。地動加速度は離散データを線形補間して評価する。
/// </summary>
class SeismicAccelLoad : public DynamicLoad
{
private:
    Vector direction_;
    TimeSeries accels_;

public:
    SeismicAccelLoad(const Vector& direction, const TimeSeries& accels)
        : direction_(direction), accels_(accels) {}

    const Vector& Direction() const { return direction_; }
    void SetAccels(const TimeSeries& accels) { accels_ = accels; }

    std::vector<NodeLoadData> load_vector(DynamicAnalysis& analysis, double t) override;
    /// 地動加速度の時刻歴
    bool has_time_series() const override { return true; }
    const TimeSeries& time_series() const override { return accels_; }
};

/// <summary>
/// 節点に作用する時刻歴荷重。空間分布 pattern(一定) × 時刻係数 factors(t) で表す。
/// 慣性力以外の一般の時刻歴外力(節点集中荷重の時刻歴)を扱う。
/// 時刻係数は離散データを線形補間して評価する。
/// </summary>
class NodalDynamicLoad : public DynamicLoad
{
private:
    std::vector<NodeLoadData> pattern_; // 空間分布(一定)
    TimeSeries factors_;                // 時刻係数

public:
    NodalDynamicLoad(double dt, const std::vector<NodeLoadData>& pattern,
                     const std::vector<double>& factors, double t0 = 0.0)
        : pattern_(pattern), factors_(dt, factors, t0) {}
    NodalDynamicLoad(const std::vector<NodeLoadData>& pattern, const TimeSeries& factors)
        : pattern_(pattern), factors_(factors) {}

    const std::vector<NodeLoadData>& Pattern() const { return pattern_; }
    void SetFactors(const TimeSeries& factors) { factors_ = factors; }

    std::vector<NodeLoadData> load_vector(DynamicAnalysis& analysis, double t) override;
    /// 時刻係数の時刻歴
    bool has_time_series() const override { return true; }
    const TimeSeries& time_series() const override { return factors_; }
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
    friend class DARecorder_KineticEnergy;
    friend class DARecorder_PotentialEnergy;
    friend class DARecorder_DampingEnergy;
    friend class DARecorder_InputEnergy;
    friend class SeismicAccelLoad;

    // 時刻 t における全荷重の合計を縮約空間へ変換する。
    //   f_reduced: [T^T·f_slave ; f_free]  (RHS用, サイズ master_dof_num + free)
    //   f_fix    : 固定DOF成分            (反力用, サイズ fixed)
    void ReducedLoadVector(double t, Eigen::VectorXd& f_reduced, Eigen::VectorXd& f_fix);

    // 静止状態からの初期加速度 a0 を M·a0 = f0 より求める(質量0のDOFは0)。
    Eigen::VectorXd ComputeInitialAcceleration(const Eigen::VectorXd& f0);

public:
    std::vector<std::shared_ptr<DynamicLoad>> loads;   // 実際に評価する時刻歴荷重(複数登録可)
    double dt = 0.0;                                   // 解析の時間刻み
    int num_steps = 0;                                 // 解析ステップ数
    std::shared_ptr<FEDynamicDampInitializer> damp_initializer = nullptr;

    int current_step = 0;

    std::vector<std::shared_ptr<DASampler>> samplers;
    // 各ステップで呼ばれるレコーダ(複数登録可)。既定でエネルギー系の4種が入っている。
    std::vector<std::shared_ptr<DARecorder>> recorders;
    bool RecordEnabled = true;

    double beta = 0.25; // 平均加速度法
                        // double beta = 1.0/6.0; // 線形加速度法(発散しがち)

    // 汎用コンストラクタ: 任意の時刻歴荷重(DynamicLoad)を与える
    DynamicAnalysis(std::shared_ptr<FEModel> model,
                    std::shared_ptr<DynamicLoad> load,
                    std::shared_ptr<FEDynamicDampInitializer> damp = nullptr);

    // 荷重を持たずに構築し、AddLoad で追加していく
    DynamicAnalysis(std::shared_ptr<FEModel> model,
                    std::shared_ptr<FEDynamicDampInitializer> damp = nullptr);

    /// 時刻歴荷重を追加する(同時に作用する荷重を重ね合わせる)
    void AddLoad(std::shared_ptr<DynamicLoad> load);
    /// 登録済みの時刻歴荷重をすべて削除する
    void ClearLoads() { loads.clear(); }

    /// レコーダを追加する(各ステップの計算後に Record が呼ばれる)
    void AddRecorder(std::shared_ptr<DARecorder> recorder);
    /// 登録済みのレコーダをすべて削除する(既定のエネルギーレコーダも外れる)
    void ClearRecorders() { recorders.clear(); }

    /// 時間刻みとステップ数を直接指定する(解析時刻は step*dt, step = 0..steps)
    void SetTimeGrid(double timestep, int steps);
    /// 時間刻みと継続時間を指定する。ステップ数は継続時間を下回らないよう切り上げる。
    void SetTimeGridByDuration(double timestep, double duration);
    /// 登録済み荷重のヒントから時間グリッドを決定する
    /// (dt = 推奨刻みの最小値、継続時間 = 終端時刻の最大値)。決定できない場合 false。
    bool SetTimeGridFromLoads();

    /// 指定ステップの時刻
    double TimeAt(int step) const { return step * dt; }
    /// 現在ステップの時刻
    double CurrentTime() const { return TimeAt(current_step); }
    /// 解析の終端時刻
    double EndTime() const { return TimeAt(num_steps); }

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
    /// 指定ステップ番号に到達するまで計算を進める
    void ComputeSteps(int steps);
    /// 最終ステップまで計算を進める
    void ComputeAll() { ComputeSteps(num_steps); }
    /// 指定時刻に到達するまで計算を進める
    void ComputeUntil(double t);

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
public:
    virtual bool Initialize(DynamicAnalysis *analysis) = 0;
    virtual bool Initialize(const FEVibrationAnalysis& vibrate_result) = 0;

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
    bool Initialize(const FEVibrationAnalysis& vibrate_result) override;

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
    bool Initialize(const FEVibrationAnalysis& vibrate_result) override;

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

    // α, βを直接指定
    FEDynamicRayleighDampInitializer(double alpha, double beta)
        : alpha(alpha), beta(beta) {}

    // 2つのモードの減衰比を指定(Initialize時の固有値解析でα, βを算出)
    FEDynamicRayleighDampInitializer(double damp_rate1, double damp_rate2, int mode1, int mode2)
        : damp_rate1(damp_rate1), damp_rate2(damp_rate2), mode1(mode1), mode2(mode2) {}

    bool Initialize(DynamicAnalysis *analysis) override;
    bool Initialize(const FEVibrationAnalysis& vibrate_result) override;

    double DampRateAtPeriod(double t) override;
};
#endif
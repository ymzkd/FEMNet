#include "FEDynamic.h"
#include "ReducedSystem.h"

#include <algorithm>
#include <cmath>

void DASampler_MaxDisplacement::Sampling(DynamicAnalysis &da)
{
    bool updated = false;
    std::vector<Displacement> disp = da.GetDisplacements();

    // 任意節点の変位が最大となるステップを記録
    for (int i = 0; i < da.model->Nodes.size(); i++)
    {
        double d_length = disp[i].Translation().norm();
        if (d_length > max_displacement)
        {
            max_displacement = d_length;
            step = da.current_step;
            updated = true;
        }
    }

    if (updated)
    {
        velocity = da.GetVelocities();
        displacement = disp;
        acceleration = da.GetAccelerations();
    }
}

void DAEnergyRecorder::Initialize()
{
    kinetic_energy.clear();
    potential_energy.clear();
    damping_energy.clear();
    input_energy.clear();

    kinetic_energy.push_back(0);
    potential_energy.push_back(0);
    damping_energy.push_back(0);
    input_energy.push_back(0);
}

void DAEnergyRecorder::RecordKineticEnergy(DynamicAnalysis &da)
{
    double energy = 0.5 * da.current_disp.dot(da.matK_aa.selfadjointView<Eigen::Upper>() * da.current_disp);
    kinetic_energy.push_back(energy);
}

void DAEnergyRecorder::RecordPotentialEnergy(DynamicAnalysis &da)
{
    double energy = 0.5 * da.current_vel.dot(da.matM_aa.selfadjointView<Eigen::Upper>() * da.current_vel);
    potential_energy.push_back(energy);
}

void DAEnergyRecorder::RecordDampingEnergy(DynamicAnalysis &da)
{
    double energy = da.current_vel.dot(da.matC_aa.selfadjointView<Eigen::Upper>() * da.current_vel);
    damping_energy.push_back(energy);
}

void DAEnergyRecorder::RecordInputEnergy(DynamicAnalysis &da)
{
    // Record は current_step 更新後に呼ばれる。既存実装に合わせ 1つ前のステップの
    // 外力を参照する。入力エネルギーは現行の符号慣行に合わせ -F_ext·v を積算する
    // (地震では F_ext=-M·ι·a_g なので (M·ι·a_g)·v となり従来と一致)。
    Eigen::VectorXd f_reduced, f_fix;
    da.ReducedLoadVector(da.TimeAt(da.current_step - 1), f_reduced, f_fix);
    input_energy.push_back(-f_reduced.dot(da.current_vel));
}

void DAEnergyRecorder::Record(DynamicAnalysis &da)
{
    RecordKineticEnergy(da);
    RecordPotentialEnergy(da);
    RecordDampingEnergy(da);
    RecordInputEnergy(da);
}

// 時刻 t における値(線形補間)。データ区間外は 0。
double TimeSeries::Value(double t) const
{
    if (IsEmpty())
        return 0.0;

    int n = static_cast<int>(values.size());
    double x = (t - t0) / dt;

    // データ区間外は無載荷。ただし t = k*dt の除算には丸め誤差(インデックス単位で
    // 1e-13 程度)が乗るため、区間端がわずかに外側へはみ出して荷重が消えることがある。
    // 区間の判定にのみ許容差を設ける(補間位置 x そのものは補正しない)。
    const double eps = 1e-12;
    if (x < -eps || x > (n - 1) + eps)
        return 0.0;

    if (x <= 0.0)
        return values[0];
    if (x >= n - 1)
        return values[n - 1];

    int i = static_cast<int>(std::floor(x));
    double s = x - i;
    return values[i] * (1.0 - s) + values[i + 1] * s;
}

// 時系列を持たない荷重が返す空の時系列
const TimeSeries& DynamicLoad::time_series() const
{
    static const TimeSeries empty;
    return empty;
}

DynamicAnalysis::DynamicAnalysis(std::shared_ptr<FEModel> model, std::shared_ptr<DynamicLoad> load, std::shared_ptr<FEDynamicDampInitializer> damp)
    : DynamicAnalysis(model, damp)
{
    AddLoad(load);
}

DynamicAnalysis::DynamicAnalysis(std::shared_ptr<FEModel> model, std::shared_ptr<FEDynamicDampInitializer> damp)
    : FEDeformOperator(model)
{
    if (damp)
    {
        damp_initializer = damp;
    }
    else
    {
        // デフォルトの減衰初期化子を用意
        damp_initializer = std::make_shared<FEDynamicStiffDampInitializer>();
    }
}

void DynamicAnalysis::AddLoad(std::shared_ptr<DynamicLoad> load)
{
    if (load)
        loads.push_back(load);
}

void DynamicAnalysis::SetTimeGrid(double timestep, int steps)
{
    dt = timestep;
    num_steps = steps;
}

void DynamicAnalysis::SetTimeGridByDuration(double timestep, double duration)
{
    dt = timestep;
    if (timestep <= 0.0 || duration <= 0.0)
    {
        num_steps = 0;
        return;
    }
    // 時間刻みは一定に保ち、継続時間を下回らないよう切り上げる
    // (端数がある場合、最終時刻は継続時間をわずかに超える)。
    // 除算誤差でステップが1つ余分に出ないよう微小トレランスを引く。
    num_steps = static_cast<int>(std::ceil(duration / timestep - 1e-9));
    if (num_steps < 1)
        num_steps = 1;
}

bool DynamicAnalysis::SetTimeGridFromLoads()
{
    // dt = 各荷重の推奨刻みの最小値、継続時間 = 終端時刻の最大値
    double min_dt = 0.0;
    double max_end = 0.0;
    for (const auto& ld : loads)
    {
        if (!ld) continue;
        double d = ld->suggested_timestep();
        if (d > 0.0 && (min_dt <= 0.0 || d < min_dt))
            min_dt = d;
        max_end = std::max(max_end, ld->end_time());
    }

    if (min_dt <= 0.0 || max_end <= 0.0)
    {
        // 定常荷重のみ等、荷重から時間グリッドを決められない
        return false;
    }

    // 解析は常に t=0 から始まるため、開始時刻を持つ荷重(TimeSeries.t0 > 0)も
    // 含めて終端時刻までを解析対象とする
    SetTimeGridByDuration(min_dt, max_end);
    return num_steps > 0;
}

// 地震荷重: 各節点の並進成分に -(m/g)·Direction·a_g を与える(回転成分は0)
std::vector<NodeLoadData> SeismicAccelLoad::load_vector(DynamicAnalysis& analysis, double t)
{
    std::vector<NodeLoadData> out;
    if (accels_.IsEmpty())
        return out;

    Vector gacc = direction_ * accels_.Value(t);

    double inv_g = 1.0 / analysis.model->GraityAccel;
    out.reserve(analysis.model->Nodes.size());
    for (size_t i = 0; i < analysis.model->Nodes.size(); i++)
    {
        // -(m/g)·Direction·a_g  (SumMassは重量なのでgで割って真の質量に変換)
        double f = -analysis.model->Nodes[i].MassData.SumMass() * inv_g;
        out.push_back(NodeLoadData(static_cast<int>(i), f * gacc.x, f * gacc.y, f * gacc.z, 0.0, 0.0, 0.0));
    }
    return out;
}

// 節点時刻歴荷重: 空間分布 × 時刻係数
std::vector<NodeLoadData> NodalDynamicLoad::load_vector(DynamicAnalysis& /*analysis*/, double t)
{
    if (factors_.IsEmpty())
        return std::vector<NodeLoadData>();

    double s = factors_.Value(t);

    std::vector<NodeLoadData> out = pattern_;
    for (NodeLoadData& nl : out)
        for (int k = 0; k < 6; k++)
            nl.loads[k] *= s;
    return out;
}

// 時刻 t における全荷重の合計を縮約空間へ変換
void DynamicAnalysis::ReducedLoadVector(double t, Eigen::VectorXd& f_reduced, Eigen::VectorXd& f_fix)
{
    // 全体DOFベクトルへ集約(登録済み荷重をすべて重ね合わせる)
    Eigen::VectorXd f_full = Eigen::VectorXd::Zero(model->Nodes.size() * NODE_DOF);
    for (const auto& ld : loads)
    {
        if (!ld) continue;
        std::vector<NodeLoadData> node_loads = ld->load_vector(*this, t);
        for (const NodeLoadData& nl : node_loads)
        {
            if (nl.id < 0) continue;
            int pos = nl.id * NODE_DOF;
            for (int k = 0; k < 6; k++)
                f_full[pos + k] += ld->Factor * nl.loads[k];
        }
    }

    // slave / free / fix へ分割
    Eigen::VectorXd f_slave(slave_indices.size());
    Eigen::VectorXd f_free(free_indices.size());
    for (size_t i = 0; i < slave_indices.size(); i++)
        f_slave(i) = f_full(slave_indices[i]);
    for (size_t i = 0; i < free_indices.size(); i++)
        f_free(i) = f_full(free_indices[i]);

    // RigidLink変換で縮約(master成分 = T^T·f_slave)
    if (master_dof_num > 0)
    {
        Eigen::VectorXd f_master = linkTransMat.transpose() * f_slave;
        f_reduced.resize(master_dof_num + free_indices.size());
        f_reduced << f_master, f_free;
    }
    else
    {
        f_reduced = f_free;
    }

    // 固定DOF成分(反力計算に使用)
    f_fix.resize(fixed_indices.size());
    for (size_t i = 0; i < fixed_indices.size(); i++)
        f_fix(i) = f_full(fixed_indices[i]);
}

// 静止状態(d=0,v=0)からの初期加速度 M·a0 = f0 を解く。
// 集中質量(対角; 回転DOFは質量0)かつ master-free はブロック対角。
// master ブロック(T^T·M11·T)は剛体回転慣性を含むためPD、質量0の自由DOFは a0=0。
Eigen::VectorXd DynamicAnalysis::ComputeInitialAcceleration(const Eigen::VectorXd& f0)
{
    int reduced_size = master_dof_num + static_cast<int>(free_indices.size());
    Eigen::VectorXd a0 = Eigen::VectorXd::Zero(reduced_size);

    if (master_dof_num > 0)
    {
        Eigen::SparseMatrix<double> Mmm = matM_aa.block(0, 0, master_dof_num, master_dof_num);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Upper> ldlt;
        ldlt.compute(Mmm);
        if (ldlt.info() == Eigen::Success)
            a0.head(master_dof_num) = ldlt.solve(f0.head(master_dof_num));
    }

    for (size_t i = 0; i < free_indices.size(); i++)
    {
        double mi = matM_aa.coeff(master_dof_num + static_cast<int>(i), master_dof_num + static_cast<int>(i));
        a0(master_dof_num + i) = (mi > 0.0) ? f0(master_dof_num + i) / mi : 0.0;
    }

    return a0;
}

bool DynamicAnalysis::Initialize()
{
    // 値を初期化
    current_step = 0;

    // 縮約系の構築（RigidLinkを考慮）
    ReducedSystem rs(*model);
    slave_indices = rs.slave_indices;
    free_indices = rs.free_indices;
    fixed_indices = rs.fixed_indices;
    linkTransMat = rs.linkTransMat;
    master_dof_num = rs.master_dof_num;

    // 時間グリッドは解析が所有する。未指定の場合のみ荷重のヒントから決定する。
    if (dt <= 0.0 || num_steps <= 0)
    {
        if (!SetTimeGridFromLoads())
        {
            std::cerr << "Dynamic analysis: time grid is not set and cannot be determined from the loads."
                << std::endl;
            return false;
        }
    }

    // 縮小空間のサイズ
    int reduced_size = rs.ReducedSize();
    current_disp = Eigen::VectorXd::Zero(reduced_size);
    current_vel = Eigen::VectorXd::Zero(reduced_size);
    current_accel = Eigen::VectorXd::Zero(reduced_size);

    // マトリクスの組み立てと縮約
    rs.Reduce(model->AssembleStiffnessMatrix(), matK_aa, &matK_ab, &matK_bb);
    rs.Reduce(model->AssembleMassMatrix(), matM_aa, &matM_ab, &matM_bb);

    // 減衰マトリクスの組み立て
    bool damp_init = damp_initializer->Initialize(this);
    if (!damp_init)
    {
        std::cerr << "Failed to initialize damping matrix." << std::endl;
        return false;
    }

    // 初期加速度: 静止状態での M·a0 = f(0) を解く(縮約行列を使うため行列組立後に実行)
    Eigen::VectorXd f0_reduced, f0_fix;
    ReducedLoadVector(0.0, f0_reduced, f0_fix);
    current_accel = ComputeInitialAcceleration(f0_reduced);

    // 因数分解しておく
    Eigen::SparseMatrix<double> compute_mat;
    compute_mat = matM_aa + 0.5 * dt * matC_aa + beta * dt * dt * matK_aa;
    solver = createSolver();
    solver->compute(compute_mat);

    // Recorder初期化
    energy_recorder.Initialize();

    return true;
}

void DynamicAnalysis::ComputeStep()
{
    // ステップ数が最大に達した場合は終了
    if (current_step >= num_steps)
    {
        std::cout << "Dynamic analysis completed." << std::endl;
        return;
    }

    // Newmarkは t_{n+1} の釣り合いを解くため、外力も t_{n+1} の値を参照する。
    // 全体節点荷重→縮約空間 f_reduced / f_fix。
    Eigen::VectorXd f_reduced, f_fix;
    ReducedLoadVector(TimeAt(current_step + 1), f_reduced, f_fix);

    // 次ステップの変位、速度、加速度を取得
    Eigen::VectorXd post_accel = f_reduced
        - matC_aa.selfadjointView<Eigen::Upper>() * (current_vel + 0.5 * dt * current_accel)
        - matK_aa.selfadjointView<Eigen::Upper>() * (current_disp + dt * current_vel + (0.5 - beta) * dt * dt * current_accel);
    post_accel = solver->solve(post_accel);
    Eigen::VectorXd post_vel = current_vel + 0.5 * (current_accel + post_accel) * dt;
    Eigen::VectorXd post_disp = current_disp + dt * current_vel + (0.5 - beta) * dt * dt * current_accel + beta * dt * dt * post_accel;

    // Update
    current_step++;
    current_accel = post_accel;
    current_vel = post_vel;
    current_disp = post_disp;

    // 反力: R = K_ab^T·d + M_ab^T·a + C_ab^T·v - F_ext,b (F_ext,b = 外力の固定DOF成分)
    Eigen::VectorXd rf = matK_ab.transpose() * current_disp + matM_ab.transpose() * current_accel +
                         matC_ab.transpose() * current_vel - f_fix;
    Eigen::VectorXd rf_full = Eigen::VectorXd::Zero(model->NodeNum() * 6);
    for (size_t i = 0; i < fixed_indices.size(); i++)
        rf_full(fixed_indices[i]) = rf(i);

    current_react_force.clear();
    for (size_t i = 0; i < model->Nodes.size(); i++)
    {
        if (!model->Nodes[i].Fix.IsAnyFix())
            continue;
        int pos = i * 6;
        current_react_force.push_back(NodeLoad(i, rf_full[pos], rf_full[pos + 1], rf_full[pos + 2], rf_full[pos + 3], rf_full[pos + 4], rf_full[pos + 5]));
    }

    if (RecordEnabled)
    {
        // Sampling
        for (auto sampler : samplers)
        {
            sampler->Sampling(*this);
        }
        // Recording
        energy_recorder.Record(*this);
    }
}

void DynamicAnalysis::ComputeSteps(int steps)
{
    // ステップ数が最大に達した場合は終了
    if (current_step >= steps)
    {
        std::cout << "Dynamic analysis completed." << std::endl;
        return;
    }

    // ステップ数分計算
    for (int i = current_step; i < steps; i++)
        ComputeStep();
}

void DynamicAnalysis::ComputeUntil(double t)
{
    if (dt <= 0.0)
        return;

    // 時刻 t を下回らない最小のステップまで進める
    int target = static_cast<int>(std::ceil(t / dt - 1e-9));
    ComputeSteps(std::min(target, num_steps));
}

bool DynamicAnalysis::SetDisplacements(std::vector<Displacement> disps)
{
    // free DOF の復元（master_dof_num 分のオフセット付き）
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        size_t idx = free_indices[i];
        size_t pos = idx % NODE_DOF;
        size_t node_id = idx / NODE_DOF;

        if (node_id >= disps.size())
            return false;

        current_disp(master_dof_num + i) = disps[node_id].displace[pos];
    }

    // master DOF の復元（RigidLink がある場合: slave変位から逆変換）
    if (master_dof_num > 0)
    {
        Eigen::VectorXd d_slave(slave_indices.size());
        for (size_t i = 0; i < slave_indices.size(); i++)
        {
            size_t idx = slave_indices[i];
            size_t pos = idx % NODE_DOF;
            size_t node_id = idx / NODE_DOF;
            if (node_id >= disps.size())
                return false;
            d_slave(i) = disps[node_id].displace[pos];
        }
        Eigen::SparseMatrix<double> TtT = linkTransMat.transpose() * linkTransMat;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(TtT);
        current_disp.head(master_dof_num) =
            solver.solve(Eigen::VectorXd(linkTransMat.transpose() * d_slave));
    }

    return true;
}

bool DynamicAnalysis::SetVelocities(std::vector<Displacement> vels)
{
    // free DOF の復元（master_dof_num 分のオフセット付き）
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        size_t idx = free_indices[i];
        size_t pos = idx % NODE_DOF;
        size_t node_id = idx / NODE_DOF;

        if (node_id >= vels.size())
            return false;

        current_vel(master_dof_num + i) = vels[node_id].displace[pos];
    }

    // master DOF の復元（RigidLink がある場合: slave速度から逆変換）
    if (master_dof_num > 0)
    {
        Eigen::VectorXd v_slave(slave_indices.size());
        for (size_t i = 0; i < slave_indices.size(); i++)
        {
            size_t idx = slave_indices[i];
            size_t pos = idx % NODE_DOF;
            size_t node_id = idx / NODE_DOF;
            if (node_id >= vels.size())
                return false;
            v_slave(i) = vels[node_id].displace[pos];
        }
        Eigen::SparseMatrix<double> TtT = linkTransMat.transpose() * linkTransMat;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(TtT);
        current_vel.head(master_dof_num) =
            solver.solve(Eigen::VectorXd(linkTransMat.transpose() * v_slave));
    }

    return true;
}

bool DynamicAnalysis::SetAccelerations(std::vector<Displacement> accs)
{
    // free DOF の復元（master_dof_num 分のオフセット付き）
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        size_t idx = free_indices[i];
        size_t pos = idx % NODE_DOF;
        size_t node_id = idx / NODE_DOF;

        if (node_id >= accs.size())
            return false;

        current_accel(master_dof_num + i) = accs[node_id].displace[pos];
    }

    // master DOF の復元（RigidLink がある場合: slave加速度から逆変換）
    if (master_dof_num > 0)
    {
        Eigen::VectorXd a_slave(slave_indices.size());
        for (size_t i = 0; i < slave_indices.size(); i++)
        {
            size_t idx = slave_indices[i];
            size_t pos = idx % NODE_DOF;
            size_t node_id = idx / NODE_DOF;
            if (node_id >= accs.size())
                return false;
            a_slave(i) = accs[node_id].displace[pos];
        }
        Eigen::SparseMatrix<double> TtT = linkTransMat.transpose() * linkTransMat;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(TtT);
        current_accel.head(master_dof_num) =
            solver.solve(Eigen::VectorXd(linkTransMat.transpose() * a_slave));
    }

    return true;
}

std::vector<Displacement> DynamicAnalysis::GetDisplacements()
{
    std::vector<Displacement> disp;
    Eigen::VectorXd d = Eigen::VectorXd::Zero(model->Nodes.size() * 6);

    if (master_dof_num > 0) {
        // RigidLinkがある場合: master → slave に展開
        Eigen::VectorXd d_master = current_disp.head(master_dof_num);
        Eigen::VectorXd d_free = current_disp.tail(free_indices.size());
        Eigen::VectorXd d_slave = linkTransMat * d_master;

        for (size_t i = 0; i < slave_indices.size(); i++)
            d(slave_indices[i]) = d_slave(i);
        for (size_t i = 0; i < free_indices.size(); i++)
            d(free_indices[i]) = d_free(i);
    }
    else {
        // RigidLinkがない場合: 従来通り
        for (size_t i = 0; i < free_indices.size(); i++)
            d(free_indices[i]) = current_disp(i);
    }

    for (size_t i = 0; i < model->Nodes.size(); i++)
    {
        int pos = i * 6;
        disp.push_back(Displacement(d[pos], d[pos + 1], d[pos + 2], d[pos + 3], d[pos + 4], d[pos + 5]));
    }

    return disp;
}

std::vector<Displacement> DynamicAnalysis::GetVelocities()
{
    std::vector<Displacement> vel;
    Eigen::VectorXd d = Eigen::VectorXd::Zero(model->Nodes.size() * 6);

    if (master_dof_num > 0) {
        // RigidLinkがある場合: master → slave に展開
        Eigen::VectorXd d_master = current_vel.head(master_dof_num);
        Eigen::VectorXd d_free = current_vel.tail(free_indices.size());
        Eigen::VectorXd d_slave = linkTransMat * d_master;

        for (size_t i = 0; i < slave_indices.size(); i++)
            d(slave_indices[i]) = d_slave(i);
        for (size_t i = 0; i < free_indices.size(); i++)
            d(free_indices[i]) = d_free(i);
    }
    else {
        // RigidLinkがない場合: 従来通り
        for (size_t i = 0; i < free_indices.size(); i++)
            d(free_indices[i]) = current_vel(i);
    }

    for (size_t i = 0; i < model->Nodes.size(); i++)
    {
        int pos = i * 6;
        vel.push_back(Displacement(d[pos], d[pos + 1], d[pos + 2], d[pos + 3], d[pos + 4], d[pos + 5]));
    }

    return vel;
}

std::vector<Displacement> DynamicAnalysis::GetAccelerations()
{
    std::vector<Displacement> acc;
    Eigen::VectorXd d = Eigen::VectorXd::Zero(model->Nodes.size() * 6);

    if (master_dof_num > 0) {
        // RigidLinkがある場合: master → slave に展開
        Eigen::VectorXd d_master = current_accel.head(master_dof_num);
        Eigen::VectorXd d_free = current_accel.tail(free_indices.size());
        Eigen::VectorXd d_slave = linkTransMat * d_master;

        for (size_t i = 0; i < slave_indices.size(); i++)
            d(slave_indices[i]) = d_slave(i);
        for (size_t i = 0; i < free_indices.size(); i++)
            d(free_indices[i]) = d_free(i);
    }
    else {
        // RigidLinkがない場合: 従来通り
        for (size_t i = 0; i < free_indices.size(); i++)
            d(free_indices[i]) = current_accel(i);
    }

    for (size_t i = 0; i < model->Nodes.size(); i++)
    {
        int pos = i * 6;
        acc.push_back(Displacement(d[pos], d[pos + 1], d[pos + 2], d[pos + 3], d[pos + 4], d[pos + 5]));
    }

    return acc;
}

BeamStressData DynamicAnalysis::GetBeamStress(int eid, double p)
{
    BarElementBase *be = dynamic_cast<BarElementBase *>(model->Elements[eid].get());

    BeamStress b_strs = be->stress(this->GetDisplacements()[be->Nodes[0]->id], this->GetDisplacements()[be->Nodes[1]->id]);
    BeamStressData strs = b_strs.Interpolate(p);
    return strs;
}

PlateStressData DynamicAnalysis::GetPlateStressData(int eid, double xi, double eta)
{
    PlateStressData data;
    std::vector<Displacement> displace = this->GetDisplacements();
    if (model->Elements[eid]->Type() == ElementType::DKT)
    {
        TriPlateElement *el = dynamic_cast<TriPlateElement *>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
                          displace[el->Nodes[2]->id], xi, eta);
    }
    else if (model->Elements[eid]->Type() == ElementType::DKQ)
    {
        QuadPlateElement *el = dynamic_cast<QuadPlateElement *>(model->Elements[eid].get());
        data = el->stress(displace[el->Nodes[0]->id], displace[el->Nodes[1]->id],
                          displace[el->Nodes[2]->id], displace[el->Nodes[3]->id], xi, eta);
    }
    return data;
}

Displacement DynamicAnalysis::GetBeamDisplace(int eid, double p)
{
    BeamElement *elm = model->GetBeamElement(eid);
    std::vector<Displacement> displace = this->GetDisplacements();
    Displacement disp = elm->DisplaceAt(
        displace[elm->Nodes[0]->id], displace[elm->Nodes[1]->id], p);

    return disp;
}

bool FEDynamicStiffDampInitializer::Initialize(DynamicAnalysis *analysis)
{
    // 解析モデルの固有振動数を計算
    FEVibrationAnalysis vib(analysis->model);
    int nconv = vib.Compute(1);
    if (nconv < 0)
    {
        std::cout << "Eigenvalue calculations did not converge." << std::endl;
        return false;
    }

    natural_angle_velocity = vib.EigenValues()[0];

    // 減衰マトリクスの組み立て
    analysis->matC_aa = analysis->matK_aa * (2.0 * damp_rate / natural_angle_velocity);
    analysis->matC_ab = analysis->matK_ab * (2.0 * damp_rate / natural_angle_velocity);
    analysis->matC_bb = analysis->matK_bb * (2.0 * damp_rate / natural_angle_velocity);

    return true;
}

bool FEDynamicStiffDampInitializer::Initialize(const FEVibrationAnalysis& vibrate_result)
{
    const std::vector<double>& eigs = vibrate_result.EigenValues();
    if (eigs.empty())
    {
        std::cerr << "Damping init: no eigenvalues available." << std::endl;
        return false;
    }

    natural_angle_velocity = eigs[0];
    return true;
}

double FEDynamicStiffDampInitializer::DampRateAtPeriod(double t)
{
    // ζ(ω) = ζ0・ω/ω1
    if (t <= 0.0 || natural_angle_velocity <= 0.0)
        return -1.0;
    return damp_rate * (2.0 * PI / t) / natural_angle_velocity;
}

bool FEDynamicMassDampInitializer::Initialize(DynamicAnalysis *analysis)
{
    // 解析モデルの固有振動数を計算
    FEVibrationAnalysis vib(analysis->model);
    int nconv = vib.Compute(1);
    if (nconv < 1)
    {
        std::cout << "Eigenvalue calculations did not converge." << std::endl;
        return false;
    }

    // C = 2ζω1・M (1次モードで減衰比ζとなる質量比例減衰)
    natural_angle_velocity = vib.EigenValues()[0];

    // 減衰マトリクスの組み立て
    double coef = 2.0 * damp_rate * natural_angle_velocity;
    analysis->matC_aa = analysis->matM_aa * coef;
    analysis->matC_ab = analysis->matM_ab * coef;
    analysis->matC_bb = analysis->matM_bb * coef;

    return true;
}

bool FEDynamicMassDampInitializer::Initialize(const FEVibrationAnalysis& vibrate_result)
{
    const std::vector<double>& eigs = vibrate_result.EigenValues();
    if (eigs.empty())
    {
        std::cerr << "Damping init: no eigenvalues available." << std::endl;
        return false;
    }

    natural_angle_velocity = eigs[0];
    return true;
}

double FEDynamicMassDampInitializer::DampRateAtPeriod(double t)
{
    // ζ(ω) = ζ0・ω1/ω
    if (t <= 0.0 || natural_angle_velocity <= 0.0)
        return -1.0;
    return damp_rate * natural_angle_velocity / (2.0 * PI / t);
}

bool FEDynamicRayleighDampInitializer::Initialize(DynamicAnalysis *analysis)
{
    if (mode1 < 1 || mode2 < 1 || mode1 == mode2)
    {
        std::cerr << "Rayleigh damping: invalid mode numbers." << std::endl;
        return false;
    }

    // 対象モードの固有振動数を計算
    int nev = std::max(mode1, mode2);
    FEVibrationAnalysis vib(analysis->model);
    int nconv = vib.Compute(nev);
    if (nconv < nev)
    {
        std::cout << "Eigenvalue calculations did not converge." << std::endl;
        return false;
    }

    std::vector<double> eigen_values = vib.EigenValues();
    natural_angle_velocity1 = eigen_values[mode1 - 1];
    natural_angle_velocity2 = eigen_values[mode2 - 1];

    // 2つのモードで指定減衰比を満たすα, βを算出
    double w1 = natural_angle_velocity1;
    double w2 = natural_angle_velocity2;
    double denom = w2 * w2 - w1 * w1;
    if (std::abs(denom) < 1e-12)
    {
        std::cerr << "Rayleigh damping: natural angle velocities are too close." << std::endl;
        return false;
    }
    alpha = 2.0 * w1 * w2 * (damp_rate1 * w2 - damp_rate2 * w1) / denom;
    beta = 2.0 * (damp_rate2 * w2 - damp_rate1 * w1) / denom;

    // 減衰マトリクスの組み立て C = αM + βK
    analysis->matC_aa = alpha * analysis->matM_aa + beta * analysis->matK_aa;
    analysis->matC_ab = alpha * analysis->matM_ab + beta * analysis->matK_ab;
    analysis->matC_bb = alpha * analysis->matM_bb + beta * analysis->matK_bb;

    return true;
}

bool FEDynamicRayleighDampInitializer::Initialize(const FEVibrationAnalysis& vibrate_result)
{
    if (mode1 < 1 || mode2 < 1 || mode1 == mode2)
    {
        std::cerr << "Rayleigh damping: invalid mode numbers." << std::endl;
        return false;
    }

    const std::vector<double>& eigs = vibrate_result.EigenValues();
    int nev = std::max(mode1, mode2);
    if (static_cast<int>(eigs.size()) < nev)
    {
        std::cerr << "Rayleigh damping: not enough modes in vibrate result." << std::endl;
        return false;
    }

    natural_angle_velocity1 = eigs[mode1 - 1];
    natural_angle_velocity2 = eigs[mode2 - 1];

    // 2つのモードで指定減衰比を満たすα, βを算出
    double w1 = natural_angle_velocity1;
    double w2 = natural_angle_velocity2;
    double denom = w2 * w2 - w1 * w1;
    if (std::abs(denom) < 1e-12)
    {
        std::cerr << "Rayleigh damping: natural angle velocities are too close." << std::endl;
        return false;
    }
    alpha = 2.0 * w1 * w2 * (damp_rate1 * w2 - damp_rate2 * w1) / denom;
    beta = 2.0 * (damp_rate2 * w2 - damp_rate1 * w1) / denom;

    return true;
}

double FEDynamicRayleighDampInitializer::DampRateAtPeriod(double t)
{
    // ζ(ω) = α/(2ω) + βω/2
    if (t <= 0.0)
        return -1.0;
    
    double w = 2.0 * PI / t;
    return alpha / (2.0 * w) + beta * w / 2.0;
}
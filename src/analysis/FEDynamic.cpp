#include "FEDynamic.h"

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
    da.ReducedLoadVector(da.current_step - 1, f_reduced, f_fix);
    input_energy.push_back(-f_reduced.dot(da.current_vel));
}

void DAEnergyRecorder::Record(DynamicAnalysis &da)
{
    RecordKineticEnergy(da);
    RecordPotentialEnergy(da);
    RecordDampingEnergy(da);
    RecordInputEnergy(da);
}

DynamicAnalysis::DynamicAnalysis(std::shared_ptr<FEModel> model, const DynamicAccelLoad& accel_load, std::shared_ptr<FEDynamicDampInitializer> damp)
    : FEDeformOperator(model), accel_load(accel_load)
{
    // 後方互換: 地震入力DTO から地震用の時刻歴荷重を生成する
    load = std::make_shared<SeismicAccelLoad>(accel_load);

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

DynamicAnalysis::DynamicAnalysis(std::shared_ptr<FEModel> model, std::shared_ptr<DynamicLoad> load, std::shared_ptr<FEDynamicDampInitializer> damp)
    : FEDeformOperator(model), load(load)
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

// 地震荷重: Prepare で各節点の SumMass/g をキャッシュ
void SeismicAccelLoad::Prepare(DynamicAnalysis& analysis)
{
    double inv_g = 1.0 / analysis.model->GraityAccel;
    mass_over_g.assign(analysis.model->Nodes.size(), 0.0);
    for (size_t i = 0; i < analysis.model->Nodes.size(); i++)
        mass_over_g[i] = analysis.model->Nodes[i].MassData.SumMass() * inv_g;
}

// 地震荷重: 各節点の並進成分に -(m/g)·Direction·a_g を与える(回転成分は0)
std::vector<NodeLoadData> SeismicAccelLoad::load_vector(DynamicAnalysis& analysis, int step, double /*t*/)
{
    std::vector<NodeLoadData> out;
    if (accel.Accels.empty())
        return out;

    // 範囲外stepはクランプ(末尾ステップの参照や初期条件で使用)
    size_t idx = (step < 0) ? 0
        : std::min(static_cast<size_t>(step), accel.Accels.size() - 1);
    Vector gacc = accel.Direction * accel.Accels[idx];

    out.reserve(analysis.model->Nodes.size());
    for (size_t i = 0; i < analysis.model->Nodes.size(); i++)
    {
        double f = -mass_over_g[i];
        out.push_back(NodeLoadData(static_cast<int>(i), f * gacc.x, f * gacc.y, f * gacc.z, 0.0, 0.0, 0.0));
    }
    return out;
}

// 節点時刻歴荷重: 空間分布 × 時刻係数
std::vector<NodeLoadData> NodalDynamicLoad::load_vector(DynamicAnalysis& /*analysis*/, int step, double /*t*/)
{
    if (factors_.empty())
        return std::vector<NodeLoadData>();

    size_t idx = (step < 0) ? 0
        : std::min(static_cast<size_t>(step), factors_.size() - 1);
    double s = factors_[idx];

    std::vector<NodeLoadData> out = pattern_;
    for (NodeLoadData& nl : out)
        for (int k = 0; k < 6; k++)
            nl.loads[k] *= s;
    return out;
}

// 全体節点荷重ベクトルを縮約空間へ変換
void DynamicAnalysis::ReducedLoadVector(int step, Eigen::VectorXd& f_reduced, Eigen::VectorXd& f_fix)
{
    // 全体DOFベクトルへ集約
    Eigen::VectorXd f_full = Eigen::VectorXd::Zero(model->Nodes.size() * NODE_DOF);
    std::vector<NodeLoadData> node_loads = load->load_vector(*this, step, step * dt);
    for (const NodeLoadData& nl : node_loads)
    {
        if (nl.id < 0) continue;
        int pos = nl.id * NODE_DOF;
        for (int k = 0; k < 6; k++)
            f_full[pos + k] += nl.loads[k];
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

    // インデックスの取得（RigidLinkを考慮）
    slave_indices = model->RigidLinkData->SlaveDOFIndices();
    free_indices = model->FreeIndices(true);  // rigid_link=true
    fixed_indices = model->FixIndices();

    // 変換行列の取得
    linkTransMat = model->RigidLinkData->TransformationMatrix().sparseView(1e-10);
    master_dof_num = linkTransMat.cols();

    // 時間グリッドを荷重から取得(Analysisが所有)
    dt = load->timestep();
    num_steps = load->steps();

    // 縮小空間のサイズ
    int reduced_size = master_dof_num + free_indices.size();
    current_disp = Eigen::VectorXd::Zero(reduced_size);
    current_vel = Eigen::VectorXd::Zero(reduced_size);
    current_accel = Eigen::VectorXd::Zero(reduced_size);

    // マトリクスの組み立て
    if (master_dof_num > 0) {
        // RigidLinkがある場合: 3x3ブロックに分割して縮小
        Eigen::SparseMatrix<double> k11, k12, k13, k22, k23, k33;
        SparseMatrixUtils::splitMatrix3x3(model->AssembleStiffnessMatrix(),
            slave_indices, free_indices, k11, k12, k13, k22, k23, k33);

        Eigen::SparseMatrix<double> m11, m12, m13, m22, m23, m33;
        SparseMatrixUtils::splitMatrix3x3(model->AssembleMassMatrix(),
            slave_indices, free_indices, m11, m12, m13, m22, m23, m33);

        // 剛性行列の縮小
        Eigen::SparseMatrix<double> kaa, kab;
        kaa = (linkTransMat.transpose() * k11.selfadjointView<Eigen::Upper>() * linkTransMat)
              .triangularView<Eigen::Upper>();
        kab = (linkTransMat.transpose() * k12);
        SparseMatrixUtils::mergeMatrixWithResize(kaa, kab, k22, matK_aa);
        matK_ab = SparseMatrixUtils::vstack(linkTransMat.transpose() * k13, k23);
        matK_bb = k33;

        // 質量行列の縮小
        Eigen::SparseMatrix<double> maa, mab;
        maa = (linkTransMat.transpose() * m11.selfadjointView<Eigen::Upper>() * linkTransMat)
              .triangularView<Eigen::Upper>();
        mab = (linkTransMat.transpose() * m12);
        SparseMatrixUtils::mergeMatrixWithResize(maa, mab, m22, matM_aa);
        matM_ab = SparseMatrixUtils::vstack(linkTransMat.transpose() * m13, m23);
        matM_bb = m33;
    }
    else {
        // RigidLinkがない場合: 従来通り2x2分割
        SparseMatrixUtils::splitMatrixWithResize(model->AssembleStiffnessMatrix(),
            fixed_indices, matK_aa, matK_ab, matK_bb);
        SparseMatrixUtils::splitMatrixWithResize(model->AssembleMassMatrix(),
            fixed_indices, matM_aa, matM_ab, matM_bb);
    }

    // 減衰マトリクスの組み立て
    bool damp_init = damp_initializer->Initialize(this);
    if (!damp_init)
    {
        std::cerr << "Failed to initialize damping matrix." << std::endl;
        return false;
    }

    // 荷重の前処理(空間分布のキャッシュ等)
    load->Prepare(*this);

    // 初期加速度: 静止状態での M·a0 = f(0) を解く(縮約行列を使うため行列組立後に実行)
    Eigen::VectorXd f0_reduced, f0_fix;
    ReducedLoadVector(0, f0_reduced, f0_fix);
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

    // Newmarkは t_{n+1} の釣り合いを解くため、外力も t_{n+1} の値を参照する
    // (末尾ステップは load 側でクランプ)。全体節点荷重→縮約空間 f_reduced / f_fix。
    Eigen::VectorXd f_reduced, f_fix;
    ReducedLoadVector(current_step + 1, f_reduced, f_fix);

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
    std::vector<double> eigen_values;
    std::vector<std::vector<Displacement>> mode_vectors;
    int nconv = analysis->model->SolveVibration(1, eigen_values, mode_vectors);
    if (nconv < 0)
    {
        std::cout << "Eigenvalue calculations did not converge." << std::endl;
        return false;
    }

    natural_angle_velocity = eigen_values[0];
    analysis->matC_aa = analysis->matK_aa * (2.0 * damp_rate / natural_angle_velocity);
    analysis->matC_ab = analysis->matK_ab * (2.0 * damp_rate / natural_angle_velocity);
    analysis->matC_bb = analysis->matK_bb * (2.0 * damp_rate / natural_angle_velocity);

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
    std::vector<double> eigen_values;
    std::vector<std::vector<Displacement>> mode_vectors;
    int nconv = analysis->model->SolveVibration(1, eigen_values, mode_vectors);
    if (nconv < 1)
    {
        std::cout << "Eigenvalue calculations did not converge." << std::endl;
        return false;
    }

    // C = 2ζω1・M (1次モードで減衰比ζとなる質量比例減衰)
    natural_angle_velocity = eigen_values[0];
    double coef = 2.0 * damp_rate * natural_angle_velocity;
    analysis->matC_aa = analysis->matM_aa * coef;
    analysis->matC_ab = analysis->matM_ab * coef;
    analysis->matC_bb = analysis->matM_bb * coef;

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
    if (!direct_coefficients)
    {
        if (mode1 < 1 || mode2 < 1 || mode1 == mode2)
        {
            std::cerr << "Rayleigh damping: invalid mode numbers." << std::endl;
            return false;
        }

        // 対象モードの固有振動数を計算
        int nev = std::max(mode1, mode2);
        std::vector<double> eigen_values;
        std::vector<std::vector<Displacement>> mode_vectors;
        int nconv = analysis->model->SolveVibration(nev, eigen_values, mode_vectors);
        if (nconv < nev)
        {
            std::cout << "Eigenvalue calculations did not converge." << std::endl;
            return false;
        }

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
    }

    // C = αM + βK
    analysis->matC_aa = alpha * analysis->matM_aa + beta * analysis->matK_aa;
    analysis->matC_ab = alpha * analysis->matM_ab + beta * analysis->matK_ab;
    analysis->matC_bb = alpha * analysis->matM_bb + beta * analysis->matK_bb;

    return true;
}

double FEDynamicRayleighDampInitializer::DampRateAtPeriod(double t)
{
    // ζ(ω) = α/(2ω) + βω/2
    if (t <= 0.0)
        return -1.0;
    if (!direct_coefficients && natural_angle_velocity1 <= 0.0)
        return -1.0; // モード指定時はInitialize前は算定不能
    double w = 2.0 * PI / t;
    return alpha / (2.0 * w) + beta * w / 2.0;
}
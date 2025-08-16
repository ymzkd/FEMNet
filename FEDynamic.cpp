#include "FEDynamic.h"

void DASampler_MaxDisplacement::Sampling(DynamicAnalysis &da)
{
    bool updated = false;
    std::vector<Displacement> disp = da.GetDisplacements();

    // 任意�?�節点の変位が最大となるス�?�?プを記録
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
    Vector gacc = da.accel_load.Direction * da.accel_load.Accels[da.current_step - 1];
    std::vector<int> free_indices = da.model->FreeIndices();
    Eigen::VectorXd post_accel0 = Eigen::VectorXd::Zero(free_indices.size());
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        int fi = free_indices[i] % NODE_DOF;
        if (fi == 0)
            post_accel0[i] = gacc.x;
        else if (fi == 1)
            post_accel0[i] = gacc.y;
        else if (fi == 2)
            post_accel0[i] = gacc.z;
        else if (fi == 3)
            post_accel0[i] = 0.0;
        else if (fi == 4)
            post_accel0[i] = 0.0;
        else if (fi == 5)
            post_accel0[i] = 0.0;
    }

    Eigen::VectorXd post_accel = da.matM_aa.selfadjointView<Eigen::Upper>() * post_accel0;
    input_energy.push_back(post_accel.dot(da.current_vel));
}

void DAEnergyRecorder::Record(DynamicAnalysis &da)
{
    RecordKineticEnergy(da);
    RecordPotentialEnergy(da);
    RecordDampingEnergy(da);
    RecordInputEnergy(da);
}

DynamicAnalysis::DynamicAnalysis(std::shared_ptr<FEModel> model, DynamicAccelLoad accel_load, FEDynamicDampInitializer *damp)
    : FEDeformOperator(model), accel_load(accel_load)
{
    if (!damp_initializer)
    {
        // デフォルトの減衰初期化子を用意
        static FEDynamicStiffDampInitializer defaultDamp;
        damp_initializer = &defaultDamp;
    }
}

bool DynamicAnalysis::Initialize()
{
    // 値を初期化
    current_step = 0;
    current_disp = Eigen::VectorXd::Zero(model->FreeDOFNum());
    current_vel = Eigen::VectorXd::Zero(model->FreeDOFNum());
    current_accel = Eigen::VectorXd::Zero(model->FreeDOFNum());

    // マトリクスの組み立て

    // StiffnessMatrixの組み立て
    free_indices = model->FreeIndices();
    fixed_indices = model->FixIndices();
    FEModel::splitMatrixWithResize(model->AssembleStiffnessMatrix(), fixed_indices, matK_aa, matK_ab, matK_bb);
    FEModel::splitMatrixWithResize(model->AssembleMassMatrix(), fixed_indices, matM_aa, matM_ab, matM_bb);

    // 減衰マトリクスの組み立て
    bool damp_init = damp_initializer->Initialize(this);
    if (!damp_init)
    {
        std::cerr << "Failed to initialize damping matrix." << std::endl;
        return false;
    }

    // 因数分解しておく
    double dt = accel_load.timestep;
    Eigen::SparseMatrix<double> compute_mat;
    compute_mat = matM_aa + 0.5 * dt * matC_aa + beta * dt * dt * matK_aa;
    solver.compute(compute_mat);

    // Recorder初期化
    energy_recorder.Initialize();

    return true;
}

void DynamicAnalysis::ComputeStep()
{
    // ステップ数が最大に達した場合は終了
    // Note: 加速度をゼロとして続けるという選択肢もある,,,
    if (current_step >= accel_load.Accels.size())
    {
        std::cout << "Dynamic analysis completed." << std::endl;
        return;
    }

    double dt = accel_load.timestep;
    Vector gacc = accel_load.Direction * accel_load.Accels[current_step];

    std::vector<int> free_indices = model->FreeIndices();
    Eigen::VectorXd post_accel0 = Eigen::VectorXd::Zero(free_indices.size());
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        int fi = free_indices[i] % NODE_DOF;
        if (fi == 0)
            post_accel0[i] = gacc.x;
        else if (fi == 1)
            post_accel0[i] = gacc.y;
        else if (fi == 2)
            post_accel0[i] = gacc.z;
    }

    std::vector<int> fixed_indices = model->FixIndices();
    Eigen::VectorXd post_accel0_fixed = Eigen::VectorXd::Zero(fixed_indices.size());
    for (size_t i = 0; i < fixed_indices.size(); i++)
    {
        int fi = fixed_indices[i] % NODE_DOF;
        if (fi == 0)
            post_accel0_fixed[i] = gacc.x;
        else if (fi == 1)
            post_accel0_fixed[i] = gacc.y;
        else if (fi == 2)
            post_accel0_fixed[i] = gacc.z;
    }

    // 次ステップの変位、速度、加速度を取得
    Eigen::VectorXd post_accel = matM_aa.selfadjointView<Eigen::Upper>() * (-post_accel0) - matC_aa.selfadjointView<Eigen::Upper>() * (current_vel + 0.5 * dt * current_accel) - matK_aa.selfadjointView<Eigen::Upper>() * (current_disp + dt * current_vel + (0.5 - beta) * dt * dt * current_accel);
    post_accel = solver.solve(post_accel);
    Eigen::VectorXd post_vel = current_vel + 0.5 * (current_accel + post_accel) * dt;
    Eigen::VectorXd post_disp = current_disp + dt * current_vel + (0.5 - beta) * dt * dt * current_accel + beta * dt * dt * post_accel;

    // Update
    current_step++;
    current_accel = post_accel;
    current_vel = post_vel;
    current_disp = post_disp;

    // Compute reaction force
    Eigen::VectorXd rf = matK_ab.transpose() * current_disp + matM_ab.transpose() * current_accel +
                         matC_ab.transpose() * current_vel + matM_bb * post_accel0_fixed;
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
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        size_t idx = free_indices[i];
        size_t pos = idx % NODE_DOF;
        size_t node_id = idx / NODE_DOF;

        if (node_id >= disps.size())
            return false;

        current_disp(i) = disps[node_id].displace[pos];
    }
    return true;
}

bool DynamicAnalysis::SetVelocities(std::vector<Displacement> vels)
{
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        size_t idx = free_indices[i];
        size_t pos = idx % NODE_DOF;
        size_t node_id = idx / NODE_DOF;

        if (node_id >= vels.size())
            return false;

        current_vel(i) = vels[node_id].displace[pos];
    }
    return true;
}

bool DynamicAnalysis::SetAccelerations(std::vector<Displacement> accs)
{
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        size_t idx = free_indices[i];
        size_t pos = idx % NODE_DOF;
        size_t node_id = idx / NODE_DOF;

        if (node_id >= accs.size())
            return false;

        current_accel(i) = accs[node_id].displace[pos];
    }
    return true;
}

std::vector<Displacement> DynamicAnalysis::GetDisplacements()
{
    std::vector<Displacement> disp;

    // 変形データ整理
    Eigen::VectorXd d = Eigen::VectorXd::Zero(model->Nodes.size() * 6);
    // Eigen::VectorXd d(Nodes.size() * 6) = Eigen::VectorXd::;
    // d.setZero();
    for (size_t i = 0; i < free_indices.size(); i++)
    {
        d(free_indices[i]) = current_disp(i);
        // std::cout << current_disp(i);
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

    // 変形データ整理
    Eigen::VectorXd d = Eigen::VectorXd::Zero(model->Nodes.size() * 6);
    // Eigen::VectorXd d(Nodes.size() * 6) = Eigen::VectorXd::;
    // d.setZero();
    for (size_t i = 0; i < free_indices.size(); i++)
        d(free_indices[i]) = current_vel(i);

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

    // 変形データ整理
    Eigen::VectorXd d = Eigen::VectorXd::Zero(model->Nodes.size() * 6);
    // Eigen::VectorXd d(Nodes.size() * 6) = Eigen::VectorXd::;
    // d.setZero();
    for (size_t i = 0; i < free_indices.size(); i++)
        d(free_indices[i]) = current_accel(i);

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
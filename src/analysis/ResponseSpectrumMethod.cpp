#include "ResponseSpectrumMethod.h"

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseAWA(ResponseValueType vt)
{
    const auto &mode_vectors = VibrateResult.ModeVectors();
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> eigen_values = VibrateResult.EigenValues(); // omega
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    const size_t M = part_facs.size();
    const size_t N = model->NodeNum();

    // Worst Mode（エネルギー寄与が最大のモードを探す）
    size_t worst_mode_idx = 0;
    double max_energy = -1.0;
    for (size_t j = 0; j < M; j++){
        double omega_j = eigen_values[j];
        double part_j = part_facs[j];

        // 参加係数は符号を持つため絶対値で寄与の大きさを比較する
        double energy = fabs(part_j) / omega_j;
        if (energy > max_energy){
            max_energy = energy;
            worst_mode_idx = j;
        }
    }

    std::vector<double> spectrums(M);
    for (size_t i = 0; i < M; i++)
    {
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    // 最悪項
    std::vector<Displacement> responses(N);
    double fac = 0;
    const size_t idx = worst_mode_idx;
    for (size_t i = 0; i < M; i++)
        fac += pow(part_facs[i] * spectrums[i] / eigen_values[i] * eigen_values[idx], 2.0);
    fac = sqrt(fac) * 0.5;
    const auto &uj = mode_vectors[idx];
    for (size_t i = 0; i < N; i++)
    {
        const Displacement &idj = uj[i];
        responses[i] += Displacement(
            fac * idj.Dx(), fac * idj.Dy(), fac * idj.Dz(),
            fac * idj.Rx(), fac * idj.Ry(), fac * idj.Rz());
    }

    // ABS項
    for (size_t j = 0; j < M; j++)
    {
        double fac = part_facs[j] * spectrums[j] * 0.5;
        const auto &uj = mode_vectors[j];
        for (size_t i = 0; i < N; i++)
        {
            const Displacement &idj = uj[i];
            responses[i] += Displacement(
                fac * idj.Dx(), fac * idj.Dy(), fac * idj.Dz(),
                fac * idj.Rx(), fac * idj.Ry(), fac * idj.Rz());
        }
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseCQC(ResponseValueType vt)
{
    const auto& mode_vectors = VibrateResult.ModeVectors();
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    const size_t M = part_facs.size();
    const size_t N = model->NodeNum();

    std::vector<double> spectrums(M);
    for (size_t i = 0; i < M; i++)
    {
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    // correlation(j,k) と (Sa·β)_j (Sa·β)_k は対称なので上三角のみ走査。
    // 対角は1回、非対角は2回ぶん加算する。
    std::vector<Displacement> responses(N);
    const double damp2 = damping_rate * damping_rate;
    for (size_t j = 0; j < M; j++)
    {
        const auto& uj = mode_vectors[j];
        const double sb_j = spectrums[j] * part_facs[j];

        // 対角 k == j : correlation = 1 ((6)式・(7)式とも χ=1 で 1)
        {
            const double fac = sb_j * sb_j;
            for (size_t i = 0; i < N; i++)
            {
                const Displacement& idj = uj[i];
                responses[i] += Displacement(
                    fac * idj.Dx() * idj.Dx(), fac * idj.Dy() * idj.Dy(), fac * idj.Dz() * idj.Dz(),
                    fac * idj.Rx() * idj.Rx(), fac * idj.Ry() * idj.Ry(), fac * idj.Rz() * idj.Rz());
            }
        }

        // 非対角 k > j : 対称性を利用して × 2
        for (size_t k = j + 1; k < M; k++)
        {
            const auto& uk = mode_vectors[k];
            const double rjk = periods[k] / periods[j];
            // (6)式・(7)式で共通の分母項
            const double denom = pow(1.0 - rjk * rjk, 2.0) + 4.0 * damp2 * rjk * pow(1.0 + rjk, 2.0);

            double correlation;
            if (vt == ResponseValueType::Acceleration)
            {
                // 絶対加速度用：論文(7)式
                const double num = 8.0 * damp2 * (1.0 + rjk) *
                                   (1.0 - (1.0 - 4.0 * damp2) * rjk + rjk * rjk) * sqrt(rjk);
                correlation = num / ((1.0 + 4.0 * damp2) * denom);
            }
            else
            {
                // 相対変位・相対速度用：論文(6)式
                correlation = 8.0 * damp2 * (1.0 + rjk) * pow(rjk, 1.5) / denom;
            }
            const double fac = 2.0 * sb_j * spectrums[k] * part_facs[k] * correlation;

            for (size_t i = 0; i < N; i++)
            {
                const Displacement& idj = uj[i];
                const Displacement& idk = uk[i];
                responses[i] += Displacement(
                    fac * idj.Dx() * idk.Dx(), fac * idj.Dy() * idk.Dy(), fac * idj.Dz() * idk.Dz(),
                    fac * idj.Rx() * idk.Rx(), fac * idj.Ry() * idk.Ry(), fac * idj.Rz() * idk.Rz());
            }
        }
    }

    // 平方根を取る
    for (size_t j = 0; j < responses.size(); j++)
    {
        responses[j] = Displacement(
            sqrt(responses[j].Dx()), sqrt(responses[j].Dy()), sqrt(responses[j].Dz()),
            sqrt(responses[j].Rx()), sqrt(responses[j].Ry()), sqrt(responses[j].Rz()));
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseSRSS(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());
    for (size_t i = 0; i < periods.size(); i++)
    {
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    const auto& mode_vectors_srss = VibrateResult.ModeVectors();
    std::vector<Displacement> responses(model->NodeNum());
    for (size_t i = 0; i < part_facs.size(); i++)
    {
        const auto& mode_vector = mode_vectors_srss[i];

        // SRSS Method
        for (size_t j = 0; j < mode_vector.size(); j++)
        {
            Displacement d2 = spectrums[i] * part_facs[i] * mode_vector[j];
            // 変位の二乗和を計算
            d2 = Displacement(d2.Dx() * d2.Dx(), d2.Dy() * d2.Dy(), d2.Dz() * d2.Dz(),
                              d2.Rx() * d2.Rx(), d2.Ry() * d2.Ry(), d2.Rz() * d2.Rz());
            responses[j] += d2;
        }
    }

    // SRSS Methodの結果を平方根で正規化
    for (size_t j = 0; j < responses.size(); j++)
    {
        responses[j] = Displacement(
            sqrt(responses[j].Dx()), sqrt(responses[j].Dy()), sqrt(responses[j].Dz()),
            sqrt(responses[j].Rx()), sqrt(responses[j].Ry()), sqrt(responses[j].Rz()));
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseABS(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());
    for (size_t i = 0; i < periods.size(); i++)
    {
        if (vt == ResponseValueType::Displacement)
            spectrums[i] = SpectrumFunction->Displacement(periods[i]);
        else if (vt == ResponseValueType::Velocity)
            spectrums[i] = SpectrumFunction->Velocity(periods[i]);
        else // (vt == ResponseValueType::Acceleration)
            spectrums[i] = SpectrumFunction->Acceleration(periods[i]);
    }

    const auto& mode_vectors_abs = VibrateResult.ModeVectors();
    std::vector<Displacement> responses(model->NodeNum());
    for (size_t i = 0; i < part_facs.size(); i++)
    {
        const auto& mode_vector = mode_vectors_abs[i];
        // | sd x vector x beta_i |

        // ABS Method
        for (size_t j = 0; j < mode_vector.size(); j++)
        {
            // 変位の絶対値を計算
            Displacement d = spectrums[i] * part_facs[i] * mode_vector[j];
            responses[j] += Displacement(abs(d.Dx()), abs(d.Dy()), abs(d.Dz()),
                                         abs(d.Rx()), abs(d.Ry()), abs(d.Rz()));
        }
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_response(ResponseValueType vt)
{
    std::vector<Displacement> responses(model->NodeNum());
    if (MethodType == ResponseSpectrumMethodType::CQC)
        responses = calculate_responseCQC(vt);
    else if (MethodType == ResponseSpectrumMethodType::SRSS)
        responses = calculate_responseSRSS(vt);
    else if (MethodType == ResponseSpectrumMethodType::AWA)
        responses = calculate_responseAWA(vt);
    else // MethodType == ResponseSpectrumMethodType::ABS
        responses = calculate_responseABS(vt);

    return responses;
}

std::vector<NodeLoad> ResponseSpectrumMethod::calculate_react_forces(const std::vector<Displacement> &disp)
{
    const size_t N = model->NodeNum();
    Eigen::VectorXd u(N * 6);
    for (size_t i = 0; i < N; i++)
    {
        const Displacement &d = disp[i];
        u[i * 6] = d.Dx();
        u[i * 6 + 1] = d.Dy();
        u[i * 6 + 2] = d.Dz();
        u[i * 6 + 3] = d.Rx();
        u[i * 6 + 4] = d.Ry();
        u[i * 6 + 5] = d.Rz();
    }

    // AssembleStiffnessMatrix()は上三角格納
    Eigen::VectorXd r = model->AssembleStiffnessMatrix().selfadjointView<Eigen::Upper>() * u;

    std::vector<NodeLoad> reacts;
    for (size_t i = 0; i < N; i++)
    {
        if (!model->Nodes[i].Fix.IsAnyFix())
            continue;

        // 静的解析(SolveLinearStatic)と同様、固定自由度の成分のみ反力として報告する
        auto fixed = model->Nodes[i].Fix.isdof_fixed();
        double v[6];
        for (int k = 0; k < 6; k++)
            v[k] = fixed[k] ? r[i * 6 + k] : 0.0;
        reacts.push_back(NodeLoad(i, v[0], v[1], v[2], v[3], v[4], v[5]));
    }
    return reacts;
}

ResponseSpectrumMethod::ResponseSpectrumMethod(std::shared_ptr<FEModel> model,
                                               FEVibrateResult vibrate_result, Vector direction, IResponseSpectrum *spectrum_function, ResponseSpectrumMethodType type)
    : FEDeformOperator(model), VibrateResult(vibrate_result), SpectrumFunction(spectrum_function), Direction(direction), MethodType(type)
{

    Compute();
}

void ResponseSpectrumMethod::Compute()
{
    displacements = calculate_response(ResponseValueType::Displacement);
    velocities = calculate_response(ResponseValueType::Velocity);
    accelerations = calculate_response(ResponseValueType::Acceleration);
    react_forces = calculate_react_forces(displacements);
    m_computed = true;
}

std::vector<Displacement> ResponseSpectrumMethod::GetDisplacements()
{
    if (m_computed)
        return displacements;
    else
        return calculate_response(ResponseValueType::Displacement);
}

std::vector<Displacement> ResponseSpectrumMethod::GetVelocities()
{
    if (m_computed)
        return velocities;
    else
        return calculate_response(ResponseValueType::Velocity);
}

std::vector<Displacement> ResponseSpectrumMethod::GetAccelerations()
{
    if (m_computed)
        return accelerations;
    else
        return calculate_response(ResponseValueType::Acceleration);
}

std::vector<NodeLoad> ResponseSpectrumMethod::GetReactForces()
{
    if (m_computed)
        return react_forces;
    else
        return calculate_react_forces(calculate_response(ResponseValueType::Displacement));
}

BeamStressData ResponseSpectrumMethod::GetBeamStress(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("ResponseSpectrumMethod: need to call Compute()");

    BarElementBase *be = dynamic_cast<BarElementBase *>(model->Elements[eid].get());

    BeamStress b_strs = be->stress(this->GetDisplacements()[be->Nodes[0]->id], this->GetDisplacements()[be->Nodes[1]->id]);
    BeamStressData strs = b_strs.Interpolate(p);
    return strs;
}

PlateStressData ResponseSpectrumMethod::GetPlateStressData(int eid, double xi, double eta)
{
    if (!m_computed)
        throw std::runtime_error("ResponseSpectrumMethod: need to call Compute()");

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

Displacement ResponseSpectrumMethod::GetBeamDisplace(int eid, double p)
{
    if (!m_computed)
        throw std::runtime_error("ResponseSpectrumMethod: need to call Compute()");

    BeamElement *elm = model->GetBeamElement(eid);
    std::vector<Displacement> displace = this->GetDisplacements();
    Displacement disp = elm->DisplaceAt(
        displace[elm->Nodes[0]->id], displace[elm->Nodes[1]->id], p);

    return disp;
}

FELinearStaticOp ResponseSpectrumMethod::GetLinearStaticCase()
{
    std::vector<std::shared_ptr<LoadBase>> loads;
    std::vector<Displacement> accels = this->GetAccelerations();
    for (size_t i = 0; i < model->NodeNum(); i++)
    {
        // とりあえず並進だけ
        double x = accels[i].Dx() / model->GraityAccel;
        double y = accels[i].Dy() / model->GraityAccel;
        double z = accels[i].Dz() / model->GraityAccel;
        loads.push_back(std::make_shared<NodeBodyForce>(NodeBodyForce(&model->Nodes[i], x, y, z)));
    }
    return FELinearStaticOp(model, loads);
}

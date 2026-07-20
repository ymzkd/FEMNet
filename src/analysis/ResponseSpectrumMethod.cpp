#include "ResponseSpectrumMethod.h"
#include "FEDynamic.h"

#include <iostream>

size_t ResponseSpectrumMethod::worst_mode_index()
{
    // エネルギー寄与 |β_j|/ω_j が最大のモード(AWAの最悪モード)の添字を返す。
    // 参加係数は符号を持つため絶対値で寄与の大きさを比較する。
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> eigen_values = VibrateResult.EigenValues(); // omega
    size_t wid = 0;
    double max_energy = -1.0;
    for (size_t j = 0; j < part_facs.size(); j++)
    {
        double energy = fabs(part_facs[j]) / eigen_values[j];
        if (energy > max_energy)
        {
            max_energy = energy;
            wid = j;
        }
    }
    return wid;
}

size_t ResponseSpectrumMethod::max_strain_energy_mode_index()
{
    // モードひずみエネルギー sE = ½·β_s²·Sv_s² (Sv=速度応答スペクトル) が最大の添字を返す。
    // argmax は |β_s|·Sv_s の最大と同値。worst_mode_index(|β|/ω) と異なりスペクトル値に依存する。
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    size_t sid = 0;
    double max_e = -1.0;
    for (size_t j = 0; j < part_facs.size(); j++)
    {
        double sv = SpectrumFunction->Velocity(periods[j]);
        double e = fabs(part_facs[j]) * sv; // ∝ sqrt(sE)
        if (e > max_e)
        {
            max_e = e;
            sid = j;
        }
    }
    return sid;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseAWA(ResponseValueType vt)
{
    const auto &mode_vectors = VibrateResult.ModeVectors();
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> eigen_values = VibrateResult.EigenValues(); // omega
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    const size_t mode_num = part_facs.size();
    const size_t N = model->NodeNum();

    // Worst Mode（エネルギー寄与が最大のモード）
    size_t wid = worst_mode_index();

    // Zero Period Response
    Displacement dir(Direction.x, Direction.y, Direction.z);
    std::vector<double> spectrums(mode_num);
    if (vt == ResponseValueType::Displacement)
    {
        for (size_t i = 0; i < mode_num; i++)
            spectrums[i] = SpectrumFunction->displacement_factored(periods[i]);
    }
    else if (vt == ResponseValueType::Velocity)
    {
        for (size_t i = 0; i < mode_num; i++)
            spectrums[i] = SpectrumFunction->velocity_factored(periods[i]);
    }
    else
    { // (vt == ResponseValueType::Acceleration)
        for (size_t i = 0; i < mode_num; i++)
            spectrums[i] = SpectrumFunction->acceleration_factored(periods[i]);
    }

    // 最悪項
    std::vector<Displacement> responses_worst(N);
    double fac = 0;
    // const size_t wid = wid;
    if (vt == ResponseValueType::Acceleration){
        for (size_t i = 0; i < mode_num; i++)
            fac += pow(part_facs[i] * spectrums[i] / eigen_values[i] * eigen_values[wid], 2.0);
        fac = sqrt(fac) * 0.5;
    } else if(vt == ResponseValueType::Velocity){
        for (size_t i = 0; i < mode_num; i++)
            fac += pow(part_facs[i] * spectrums[i], 2.0);
        fac = sqrt(fac) * 0.5;
    } else{ // (vt == ResponseValueType::Displacement)
        for (size_t i = 0; i < mode_num; i++)
            fac += pow(part_facs[i] * spectrums[i] / eigen_values[wid] * eigen_values[i], 2.0);
        fac = sqrt(fac) * 0.5;
    }

    double irfac = 1.0;
    if (EnableRigidResponse)
    {
        double alpha = RigidResponse.RigidResponseFactor(periods[wid]);
        irfac = sqrt(1.0 - alpha * alpha);
    }

    const auto &uj = mode_vectors[wid];
    for (size_t i = 0; i < N; i++)
    {
        const Displacement &idj = uj[i];
        responses_worst[i] += Displacement(
            fac * idj.Dx(), fac * idj.Dy(), fac * idj.Dz(),
            fac * idj.Rx(), fac * idj.Ry(), fac * idj.Rz());
    }

    // 最悪項 W と組み合わせる残差項。いずれも AWA の平均化 ½ を掛ける。
    // AWA_ABS/AWA_CQC の違いは残差に使う関数だけ:
    //   AWA_ABS = calculate_responseABS / AWA_CQC = calculate_responseCQC。
    std::vector<Displacement> responses_abs(N);
    if (MethodType == ResponseSpectrumMethodType::AWA_CQC)
        responses_abs = calculate_responseCQC(vt); // 残差 = CQC 結合場
    else // ResponseSpectrumMethodType::AWA_ABS
        responses_abs = calculate_responseABS(vt); // 残差 = 絶対値和 Σ_j|β_j S_j φ_j|
    for (Displacement &r : responses_abs)
        r = 0.5 * r;                               // AWA の平均化 ½

    // 合成(A/B 共通): 成分ごとローカル符号。worst項 W の符号を基準に残差を必ず外側
    //   (絶対値が増える向き)へ加算する。R_k = W_k + sign(W_k)·|residual_k|
    //   (W≈0 の節点/成分では residual 側の符号)。
    std::vector<Displacement> responses(N);
    auto reinforce = [](double w, double a) -> double {
        double sgn = (w != 0.0) ? ((w > 0.0) ? 1.0 : -1.0) : ((a >= 0.0) ? 1.0 : -1.0);
        return w + sgn * fabs(a);
    };
    for (size_t i = 0; i < N; i++)
    {
        const Displacement &w = responses_worst[i];
        const Displacement &a = responses_abs[i];
        responses[i] = Displacement(
            reinforce(w.Dx(), a.Dx()), reinforce(w.Dy(), a.Dy()), reinforce(w.Dz(), a.Dz()),
            reinforce(w.Rx(), a.Rx()), reinforce(w.Ry(), a.Ry()), reinforce(w.Rz(), a.Rz()));
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseCQC(ResponseValueType vt)
{
    const auto& mode_vectors = VibrateResult.ModeVectors();
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    const size_t mode_num = part_facs.size();
    const size_t N = model->NodeNum();
    
    // Zero Period Response
    std::vector<double> spectrums(mode_num);
    if (vt == ResponseValueType::Displacement)
    {
        for (size_t i = 0; i < mode_num; i++)
            spectrums[i] = SpectrumFunction->displacement_factored(periods[i]);
    }
    else if (vt == ResponseValueType::Velocity)
    {
        for (size_t i = 0; i < mode_num; i++)
            spectrums[i] = SpectrumFunction->velocity_factored(periods[i]);
    }
    else
    { // (vt == ResponseValueType::Acceleration)
        for (size_t i = 0; i < mode_num; i++)
            spectrums[i] = SpectrumFunction->acceleration_factored(periods[i]);
    }

    // correlation(j,k) と (Sa·β)_j (Sa·β)_k は対称なので上三角のみ走査。
    // 対角は1回、非対角は2回ぶん加算する。
    std::vector<Displacement> responses(N);
    const double damp2 = damping_rate * damping_rate;
    for (size_t j = 0; j < mode_num; j++)
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
        for (size_t k = j + 1; k < mode_num; k++)
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
            
            double irfac_j = 1.0;
            double irfac_k = 1.0;
            if (EnableRigidResponse)
            {
                irfac_j = sqrt(1.0 - pow(RigidResponse.RigidResponseFactor(periods[j]), 2.0));
                irfac_k = sqrt(1.0 - pow(RigidResponse.RigidResponseFactor(periods[k]), 2.0));
            }
            const double fac = 2.0 * sb_j * spectrums[k] * part_facs[k] * correlation * irfac_j * irfac_k;

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

    // 平方根を取る(分散 = xᵀCx は C が半正定値なので理論上非負。
    //  加力方向に直交する≈ゼロ成分では対角が相関項でほぼ完全に相殺され、
    //  残差が丸め誤差(~1e-16)で微小負になり sqrt が NaN になり得るため 0 にクランプする)
    auto sqrt_clamp = [](double v) { return sqrt(v > 0.0 ? v : 0.0); };
    for (size_t j = 0; j < responses.size(); j++)
    {
        responses[j] = Displacement(
            sqrt_clamp(responses[j].Dx()), sqrt_clamp(responses[j].Dy()), sqrt_clamp(responses[j].Dz()),
            sqrt_clamp(responses[j].Rx()), sqrt_clamp(responses[j].Ry()), sqrt_clamp(responses[j].Rz()));
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_responseSRSS(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());

    // Zero Period Response
    if (vt == ResponseValueType::Displacement)
    {
        for (size_t i = 0; i < periods.size(); i++)
            spectrums[i] = SpectrumFunction->displacement_factored(periods[i]);
    }
    else if (vt == ResponseValueType::Velocity)
    {
        for (size_t i = 0; i < periods.size(); i++)
            spectrums[i] = SpectrumFunction->velocity_factored(periods[i]);
    }
    else
    { // (vt == ResponseValueType::Acceleration)
        for (size_t i = 0; i < periods.size(); i++)
            spectrums[i] = SpectrumFunction->acceleration_factored(periods[i]);
    }

    const auto& mode_vectors_srss = VibrateResult.ModeVectors();
    std::vector<Displacement> responses(model->NodeNum());
    for (size_t i = 0; i < part_facs.size(); i++)
    {
        const auto& mode_vector = mode_vectors_srss[i];

        // SRSS Method
        for (size_t j = 0; j < mode_vector.size(); j++)
        {
            double irfac = 1.0;
            if (EnableRigidResponse)
            {
                double alpha = RigidResponse.RigidResponseFactor(periods[i]);
                irfac = sqrt(1.0 - alpha * alpha);
            }
            Displacement d2 = irfac * spectrums[i] * part_facs[i] * mode_vector[j];
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
    
    // Zero Period Response
    if (vt == ResponseValueType::Displacement){
        for (size_t i = 0; i < periods.size(); i++)
            spectrums[i] = SpectrumFunction->displacement_factored(periods[i]);
    }
    else if (vt == ResponseValueType::Velocity){
        for (size_t i = 0; i < periods.size(); i++)
            spectrums[i] = SpectrumFunction->velocity_factored(periods[i]);
    }
    else { // (vt == ResponseValueType::Acceleration)
        for (size_t i = 0; i < periods.size(); i++)
            spectrums[i] = SpectrumFunction->acceleration_factored(periods[i]);
    }

    const auto& mode_vectors_abs = VibrateResult.ModeVectors();
    std::vector<Displacement> responses(model->NodeNum());
    for (size_t i = 0; i < part_facs.size(); i++)
    {
        const auto& mode_vector = mode_vectors_abs[i];

        // ABS Method  | sd x vector x beta_i |
        for (size_t j = 0; j < mode_vector.size(); j++)
        {
            // 変位の絶対値を計算
            double irfac = 1.0;
            if (EnableRigidResponse){
                double alpha = RigidResponse.RigidResponseFactor(periods[i]);
                irfac = sqrt(1.0 - alpha * alpha);
            }

            Displacement d = irfac * spectrums[i] * part_facs[i] * mode_vector[j];
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
    else if (MethodType == ResponseSpectrumMethodType::AWA_ABS ||
             MethodType == ResponseSpectrumMethodType::AWA_CQC)
        responses = calculate_responseAWA(vt);
    else // MethodType == ResponseSpectrumMethodType::ABS
        responses = calculate_responseABS(vt);

    std::vector<Displacement> r_responses;
    if (EnableRigidResponse)
    {
        r_responses = calculate_rigid_response(vt);
        for (size_t i = 0; i < responses.size(); i++){
            const Displacement &r = r_responses[i];
            const Displacement &d = responses[i];
            responses[i] = Displacement(
                sqrt(r.Dx() * r.Dx() + d.Dx() * d.Dx()),
                sqrt(r.Dy() * r.Dy() + d.Dy() * d.Dy()),
                sqrt(r.Dz() * r.Dz() + d.Dz() * d.Dz()),
                sqrt(r.Rx() * r.Rx() + d.Rx() * d.Rx()),
                sqrt(r.Ry() * r.Ry() + d.Ry() * d.Ry()),
                sqrt(r.Rz() * r.Rz() + d.Rz() * d.Rz())
            );
        }
    }

    // 符号調整: responses の各成分の絶対値は保持したまま、符号を「符号源」の各成分の符号に合わせる。
    // 符号源は次のいずれか(節点ごとの Displacement 列):
    // - SIGN_NONE               : 符号調整なし(そのまま)
    // - SIGN_SPECIFIED_MODE     : 指定次数(sign_mode_index, 0始まり)のモード形状
    // - SIGN_WORST_MODE         : AWAの最悪モード(|β|/ω 最大)の形状。採用次数を sign_mode_index に記録。
    // - SIGN_STRAIN_ENERGY_MODE : モードひずみエネルギー(β²·Sv² 最大)の形状。採用次数を sign_mode_index に記録。
    if (sign_type != ResponseSignType::SIGN_NONE)
    {
        const std::vector<Displacement> *sign_src = nullptr; // 符号源(節点ごとの Displacement)

        const auto &mode_vectors = VibrateResult.ModeVectors();
        int sid;
        if (sign_type == ResponseSignType::SIGN_WORST_MODE)
        {
            sid = static_cast<int>(worst_mode_index());
            sign_mode_index = sid; // 採用した最悪モード次数を記録(出力)
        }
        else if (sign_type == ResponseSignType::SIGN_STRAIN_ENERGY_MODE)
        {
            sid = static_cast<int>(max_strain_energy_mode_index());
            sign_mode_index = sid; // 採用したひずみエネルギー最大モード次数を記録(出力)
        }
        else // SIGN_SPECIFIED_MODE
        {
            sid = sign_mode_index; // ユーザが指定した次数(入力, 0始まり)
        }
        if (sid >= 0 && static_cast<size_t>(sid) < mode_vectors.size())
            sign_src = &mode_vectors[sid];

        if (sign_src != nullptr && sign_src->size() == responses.size())
        {
            const auto &um = *sign_src;
            auto sgn = [](double v) { return (v < 0.0) ? -1.0 : 1.0; };
            for (size_t i = 0; i < responses.size(); i++)
            {
                const Displacement &s = um[i];
                const Displacement &r = responses[i];
                responses[i] = Displacement(
                    fabs(r.Dx()) * sgn(s.Dx()), fabs(r.Dy()) * sgn(s.Dy()), fabs(r.Dz()) * sgn(s.Dz()),
                    fabs(r.Rx()) * sgn(s.Rx()), fabs(r.Ry()) * sgn(s.Ry()), fabs(r.Rz()) * sgn(s.Rz()));
            }
        }
    }

    return responses;
}

std::vector<Displacement> ResponseSpectrumMethod::calculate_rigid_response(ResponseValueType vt)
{
    std::vector<double> part_facs = VibrateResult.ParticipationFactors(Direction);
    std::vector<double> periods = VibrateResult.NaturalPeriods();
    std::vector<double> spectrums(periods.size());

    // Zero Period Response
    double s0 = 0.0;
    Displacement dir(Direction.x, Direction.y, Direction.z);
    if (vt == ResponseValueType::Displacement)
    {
        s0 = SpectrumFunction->displacement_factored(0.0);
        for (size_t i = 0; i < periods.size(); i++)
            spectrums[i] = SpectrumFunction->displacement_factored(periods[i]);
    }
    else if (vt == ResponseValueType::Velocity)
    {
        s0 = SpectrumFunction->velocity_factored(0.0);
        for (size_t i = 0; i < periods.size(); i++)
            spectrums[i] = SpectrumFunction->velocity_factored(periods[i]);
    }
    else
    { // (vt == ResponseValueType::Acceleration)
        s0 = SpectrumFunction->acceleration_factored(0.0);
        for (size_t i = 0; i < periods.size(); i++)
            spectrums[i] = SpectrumFunction->acceleration_factored(periods[i]);
    }

    const auto &mode_vectors_abs = VibrateResult.ModeVectors();
    const size_t N = model->NodeNum();
    std::vector<Displacement> r_responses(N);
    std::vector<Displacement> modal_influence(N); // Σ_i β_i·φ_i (影響ベクトルのモード展開)

    // 剛成分は同相のためモード間で代数和。欠落質量残差 (dir − Σβφ)·s0 はループ外で1回だけ加算する。
    // (以前は欠落質量項をモードループ内に置いていたため dir·s0 が M 回加算され、
    //  加速度応答が M×ZPA に膨れるバグがあった。)
    for (size_t i = 0; i < part_facs.size(); i++)
    {
        const auto &mode_vector = mode_vectors_abs[i];
        double rfac = 1.0;
        if (EnableRigidResponse)
            rfac = RigidResponse.RigidResponseFactor(periods[i]);

        for (size_t j = 0; j < N; j++)
        {
            // 剛成分: α_i·β_i·S_i·φ_ij (代数和)
            r_responses[j] += rfac * (spectrums[i] * part_facs[i] * mode_vector[j]);
            // 影響ベクトルのモード寄与 Σβφ を蓄積(欠落質量の算出用)
            modal_influence[j] += part_facs[i] * mode_vector[j];
        }
    }

    // 欠落質量(残差)応答: (影響ベクトル dir − Σβφ)·s0 を1回だけ加算
    for (size_t j = 0; j < N; j++)
        r_responses[j] += (dir - modal_influence[j]) * s0;

    return r_responses;
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
    // 減衰考慮が有効かつ初期化子が設定されている場合のみ初期化(既定は無効/null)。
    // 初期化失敗(モード数不足など)時は誤った減衰補正を避けるため減衰を無効化する。
    if (SpectrumFunction->enable_damp_factor && SpectrumFunction->DampInitializer != nullptr)
    {
        if (!SpectrumFunction->DampInitializer->Initialize(VibrateResult))
        {
            std::cerr << "ResponseSpectrumMethod: 減衰初期化に失敗したため減衰補正を無効化します。" << std::endl;
            SpectrumFunction->enable_damp_factor = false;
        }
    }

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

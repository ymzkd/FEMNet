// DynamicAnalysis::Rewind() の検証
//
// Rewind() は組み立て済みの系(縮約情報・質量/剛性/減衰マトリクス・因数分解)を
// そのまま使い、状態量をステップ0へ戻す。Initialize() から
// 「系の再構築」だけを取り除いたものなので、次を検証する。
//   (1) 未初期化の解析では false を返す(何もしない)
//   (2) 巻き戻し後の状態量・反力が Initialize() の場合とビット単位で一致する
//   (3) 巻き戻した後に計算し直すと Initialize() の場合と同じ応答が得られる
//   (4) RecordEnabled == false のとき、Initialize() と同じく記録を保持する
//   (5) RecordEnabled == true のとき、Initialize() と同じく記録を取り直す

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <memory>

#include "Model.h"
#include "Elements/Elements.h"
#include "FEDynamic.h"
#include "DASampler.h"
#include "DARecorder.h"

namespace {

int g_failures = 0;

void Check(bool cond, const std::string &msg)
{
    if (!cond)
    {
        std::cout << "  [FAIL] " << msg << "\n";
        g_failures++;
    }
    else
    {
        std::cout << "  [ ok ] " << msg << "\n";
    }
}

// 片持ち梁(x方向に n 分割)。質量を持たせるため密度を与える。
std::shared_ptr<FEModel> BuildCantileverModel(double length, int n)
{
    auto model = std::make_shared<FEModel>();

    double dl = length / n;
    model->Nodes.push_back(Node(0, 0, 0, 0));
    model->Nodes[0].Fix.FixAll();
    for (int i = 0; i < n; i++)
        model->Nodes.push_back(Node(i + 1, dl * (i + 1), 0, 0));

    Material mat(5000.0, 0.2);
    mat.dense = 5.0 / 1000.0 / 1000.0;
    Section sec(100, 833.33, 833.33, 1406.25);

    model->Materials.push_back(mat);
    model->Sections.push_back(sec);

    for (int i = 0; i < n; i++)
        model->Elements.push_back(std::make_shared<BeamElement>(
            i, &model->Nodes[i], &model->Nodes[i + 1], &model->Sections[0], model->Materials[0]));

    model->ComputeElementNodeMass();
    return model;
}

// 先端を鉛直に揺らす節点時刻歴荷重
std::shared_ptr<NodalDynamicLoad> BuildLoad(int tip_node, double dt, int steps)
{
    std::vector<NodeLoadData> pattern{NodeLoadData(tip_node, 0, 0, -100.0)};
    std::vector<double> factors(steps + 1);
    for (int i = 0; i <= steps; i++)
        factors[i] = std::sin(2.0 * PI * (i * dt) / 0.4);
    return std::make_shared<NodalDynamicLoad>(dt, pattern, factors);
}

// 状態量(変位・速度・加速度・反力)をひとつの列に並べる。
// 巻き戻しの一致はビット単位で成立するはずなので、丸めずそのまま比較する。
std::vector<double> StateVector(DynamicAnalysis &da)
{
    std::vector<double> out;
    auto push_disps = [&out](const std::vector<Displacement> &ds) {
        for (const Displacement &v : ds)
        {
            out.push_back(v.Dx()); out.push_back(v.Dy()); out.push_back(v.Dz());
            out.push_back(v.Rx()); out.push_back(v.Ry()); out.push_back(v.Rz());
        }
    };
    push_disps(da.GetDisplacements());
    push_disps(da.GetVelocities());
    push_disps(da.GetAccelerations());
    for (NodeLoad &nl : da.GetReactForces())
        for (int k = 0; k < 6; k++)
            out.push_back(nl.data.loads[k]);
    return out;
}

bool Identical(const std::vector<double> &a, const std::vector<double> &b)
{
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); i++)
        if (a[i] != b[i])
            return false;
    return true;
}

// 最大相対差(一致しない場合の報告用)
double MaxRelDiff(const std::vector<double> &a, const std::vector<double> &b)
{
    if (a.size() != b.size())
        return 1.0;
    double worst = 0.0;
    for (size_t i = 0; i < a.size(); i++)
    {
        if (a[i] == b[i])
            continue;
        double denom = std::max(std::abs(a[i]), 1e-30);
        worst = std::max(worst, std::abs(a[i] - b[i]) / denom);
    }
    return worst;
}

const double kDt = 0.005;
const int kSteps = 400;
const int kTip = 4;

// 解析を1つ組み立てる(サンプラーとレコーダを1つずつ登録する)
std::shared_ptr<DynamicAnalysis> BuildAnalysis(std::shared_ptr<FEModel> model,
                                               std::shared_ptr<DASampler_MaxResponse> &sampler,
                                               std::shared_ptr<DARecorder_KineticEnergy> &recorder)
{
    auto damp = std::make_shared<FEDynamicStiffDampInitializer>(0.05);
    auto da = std::make_shared<DynamicAnalysis>(model, BuildLoad(kTip, kDt, kSteps), damp);
    da->SetTimeGrid(kDt, kSteps);

    sampler = std::make_shared<DASampler_MaxResponse>(ResponseValueType::Displacement);
    da->samplers.push_back(sampler);

    da->ClearRecorders();
    recorder = std::make_shared<DARecorder_KineticEnergy>();
    da->recorders.push_back(recorder);

    return da;
}

} // namespace

int main()
{
    std::cout << "DynamicAnalysis::Rewind test\n\n";

    auto model = BuildCantileverModel(200.0, 4);

    // --- (1) 未初期化の解析では false を返す ---
    std::cout << "[1] 未初期化の解析\n";
    {
        std::shared_ptr<DASampler_MaxResponse> sampler;
        std::shared_ptr<DARecorder_KineticEnergy> recorder;
        auto da = BuildAnalysis(model, sampler, recorder);
        Check(!da->Rewind(), "Initialize() 前の Rewind() は false を返す");
        Check(da->Initialize(), "Initialize() は成功する");
        Check(da->Rewind(), "Initialize() 後の Rewind() は true を返す");
    }

    // --- (2)(3)(4) RecordEnabled = false: 記録を保持したまま巻き戻す ---
    // DynamicSolverComponent が計算後に行っている使い方をなぞる。
    std::cout << "\n[2] RecordEnabled = false (記録を保持する巻き戻し)\n";
    {
        std::shared_ptr<DASampler_MaxResponse> sampler_i, sampler_r;
        std::shared_ptr<DARecorder_KineticEnergy> recorder_i, recorder_r;

        // 基準: Initialize() で巻き戻す
        auto da_i = BuildAnalysis(model, sampler_i, recorder_i);
        da_i->Initialize();
        da_i->ComputeAll();
        int sampled_step_i = sampler_i->step;
        double sampled_max_i = sampler_i->MaxValue;
        std::vector<double> records_i = recorder_i->Values;
        da_i->RecordEnabled = false;
        da_i->Initialize();
        std::vector<double> state_i = StateVector(*da_i);

        // 比較: Rewind() で巻き戻す
        auto da_r = BuildAnalysis(model, sampler_r, recorder_r);
        da_r->Initialize();
        da_r->ComputeAll();
        da_r->RecordEnabled = false;
        Check(da_r->Rewind(), "Rewind() が成功する");
        std::vector<double> state_r = StateVector(*da_r);

        Check(da_i->current_step == 0 && da_r->current_step == 0, "どちらも current_step = 0 に戻る");
        bool same_state = Identical(state_i, state_r);
        Check(same_state, "巻き戻し後の変位・速度・加速度・反力が完全一致する");
        if (!same_state)
            std::cout << "         最大相対差: " << MaxRelDiff(state_i, state_r) << "\n";

        Check(sampler_r->step == sampled_step_i && sampler_r->MaxValue == sampled_max_i,
              "サンプラーの記録が保持される");
        Check(Identical(recorder_r->Values, records_i),
              "レコーダの記録が保持される(件数・値とも)");
        Check(recorder_r->Values.size() == static_cast<size_t>(kSteps) + 1,
              "レコーダの件数がステップ数+1のまま");

        // (3) 巻き戻した後に計算し直すと同じ応答になる
        da_i->ComputeAll();
        da_r->ComputeAll();
        bool same_recompute = Identical(StateVector(*da_i), StateVector(*da_r));
        Check(same_recompute, "巻き戻し後に再計算した最終状態が完全一致する");
        if (!same_recompute)
            std::cout << "         最大相対差: "
                      << MaxRelDiff(StateVector(*da_i), StateVector(*da_r)) << "\n";
    }

    // --- (5) RecordEnabled = true: 記録を取り直す ---
    std::cout << "\n[3] RecordEnabled = true (記録を取り直す巻き戻し)\n";
    {
        std::shared_ptr<DASampler_MaxResponse> sampler;
        std::shared_ptr<DARecorder_KineticEnergy> recorder;
        auto da = BuildAnalysis(model, sampler, recorder);
        da->Initialize();
        da->ComputeAll();
        size_t count_after_run = recorder->Values.size();

        Check(da->Rewind(), "Rewind() が成功する");
        Check(recorder->Values.size() == 1,
              "レコーダのバッファが初期化される(step 0 の 1 件のみ)");
        Check(count_after_run == static_cast<size_t>(kSteps) + 1,
              "巻き戻し前は全ステップぶん記録されていた");

        da->ComputeAll();
        Check(recorder->Values.size() == static_cast<size_t>(kSteps) + 1,
              "再計算で同じ件数まで記録し直される");
    }

    // --- 参考: 巻き戻しにかかる時間の比較(検証ではなく情報表示) ---
    std::cout << "\n[参考] 巻き戻し所要時間\n";
    {
        std::shared_ptr<DASampler_MaxResponse> sampler;
        std::shared_ptr<DARecorder_KineticEnergy> recorder;
        auto da = BuildAnalysis(model, sampler, recorder);
        da->Initialize();
        da->RecordEnabled = false;

        const int trials = 20;
        auto measure = [&](bool use_rewind) {
            double best = 1e30;
            for (int i = 0; i < trials; i++)
            {
                da->ComputeSteps(50);
                auto t0 = std::chrono::steady_clock::now();
                if (use_rewind) da->Rewind(); else da->Initialize();
                auto t1 = std::chrono::steady_clock::now();
                best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
            }
            return best;
        };
        double ms_init = measure(false);
        double ms_rewind = measure(true);
        std::cout << "  Initialize(): " << ms_init << " ms\n";
        std::cout << "  Rewind()    : " << ms_rewind << " ms\n";
    }

    std::cout << "\n";
    if (g_failures == 0)
    {
        std::cout << "RESULT: PASS (all checks)\n";
        return 0;
    }
    std::cout << "RESULT: FAIL (" << g_failures << " checks failed)\n";
    return 1;
}

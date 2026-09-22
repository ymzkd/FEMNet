// 支点ばね要素(SupportSpringElement)の検証テスト
//   1) 1自由度ばね: u = P/k、反力 R = -P
//   2) ばね支持された片持ち梁: 並列剛性 k + k_beam
//   3) ばね支持の固有値: omega = sqrt(k/m)
//   4) ばね支持と固定支持の混在: 反力の合計が荷重と釣り合う
//   5) 減衰のあるばね支持の動的反力: 外力 + 反力 + 慣性力 = 0
//   6) 応答スペクトル法でのばね支点の反力: R = -k・u
//   7) 固定した自由度にばねを重ねても固定反力が失われない
//   8) 1節点に複数のばね要素: 剛性と反力が足し合わされる
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Model.h"
#include "Elements/Elements.h"
#include "FELinearStaticOp.h"
#include "FEVibrationAnalysis.h"
#include "FEDynamic.h"
#include "ResponseSpectrumMethod.h"

namespace
{
int g_failed = 0;

void Check(bool cond, const std::string &what)
{
    std::cout << (cond ? "  [OK]   " : "  [FAIL] ") << what << std::endl;
    if (!cond)
        g_failed++;
}

// 相対誤差で比較する
void CheckNear(double actual, double expected, double rel_tol, const std::string &what)
{
    double denom = std::abs(expected) > 1e-30 ? std::abs(expected) : 1.0;
    double err = std::abs(actual - expected) / denom;
    std::cout << (err <= rel_tol ? "  [OK]   " : "  [FAIL] ") << what
              << " : actual=" << actual << " expected=" << expected
              << " rel.err=" << err << std::endl;
    if (err > rel_tol)
        g_failed++;
}

// 1自由度だけにばね定数を持つ支点ばね要素を追加する。
// 要素は節点をポインタで参照するため、節点をすべて追加してから呼ぶこと。
void AddSpring(FEModel &m, int node, int dof, double k)
{
    std::array<double, 6> ks{};
    ks[dof] = k;
    m.add_element(SupportSpringElement((int)m.Elements.size(), &m.Nodes[node],
                                       ks[0], ks[1], ks[2], ks[3], ks[4], ks[5]));
}

// 1) 1自由度ばね: 節点1つをばねだけで支え、鉛直荷重を与える
void TestSingleSpring()
{
    std::cout << "[1] 1自由度ばね" << std::endl;
    const double k = 250.0;   // N/mm
    const double P = -1000.0; // N (下向き)

    auto m = std::make_shared<FEModel>();
    m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
    for (int d = 0; d < 6; d++)
        m->Nodes[0].Fix.BoundaryTypes[d] = ConstraintType::Fix;
    m->Nodes[0].Fix.BoundaryTypes[2] = ConstraintType::Free;
    AddSpring(*m, 0, 2, k);

    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.push_back(std::make_shared<NodeLoad>(0, 0.0, 0.0, P));

    FELinearStaticOp op(m, loads);
    op.Compute();

    CheckNear(op.GetDisplacements()[0].Dz(), P / k, 1e-10, "変位 u = P/k");
    std::vector<NodeLoad> react = op.GetReactForces();
    Check(react.size() == 1, "反力の節点数 = 1");
    if (!react.empty())
        CheckNear(react[0].Pz(), -P, 1e-10, "ばね反力 R = -P");
}

// 2) 先端ばね支持の片持ち梁: ばねと梁が並列に効く等価剛性と一致するか
void TestSpringWithBeam()
{
    std::cout << "[2] ばね支持された片持ち梁(並列剛性)" << std::endl;
    const double L = 2000.0, E = 205000.0, I = 1.0e7;
    const double k_beam = 3.0 * E * I / (L * L * L); // 片持ち梁の先端剛性
    const double k = 0.5 * k_beam;                   // ばねは梁の半分の剛性
    const double P = -1000.0;

    auto m = std::make_shared<FEModel>();
    m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
    m->Nodes.push_back(Node(1, L, 0.0, 0.0));
    m->Materials.push_back(Material(E, 0.3));
    m->Sections.push_back(Section(1000.0, I, I, 1.0e7));
    m->Nodes[0].Fix.FixAll();
    // 先端: 鉛直をばね支持、面外方向とねじりは固定して面内の曲げのみ扱う
    for (int d : {0, 1, 3, 5})
        m->Nodes[1].Fix.BoundaryTypes[d] = ConstraintType::Fix;
    m->add_beam_element(0, 0, 1, 0, 0, 0.0);
    AddSpring(*m, 1, 2, k);

    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.push_back(std::make_shared<NodeLoad>(1, 0.0, 0.0, P));

    FELinearStaticOp op(m, loads);
    op.Compute();

    // ばねと梁は並列(同じ節点を両方が支える)なので剛性は和
    CheckNear(op.GetDisplacements()[1].Dz(), P / (k + k_beam), 1e-6, "先端変位 u = P/(k + k_beam)");

    double r_spring = 0.0, r_fix = 0.0;
    std::vector<NodeLoad> react_forces = op.GetReactForces();
    for (NodeLoad &nl : react_forces)
        (nl.id == 1 ? r_spring : r_fix) = nl.Pz();
    CheckNear(r_spring, -k * (P / (k + k_beam)), 1e-6, "ばね反力 R = -k・u");
    CheckNear(r_spring + r_fix, -P, 1e-6, "反力の合計が荷重と釣り合う");
}

// 3) ばね支持の固有値: omega = sqrt(k/m)
// Spectraは求めるモード数より十分に大きな行列を要するため、互いに独立した
// 1質点系(ばね定数 k*(i+1)^2)を5つ並べ、小さい方から3モードを検証する。
void TestSpringVibration()
{
    std::cout << "[3] ばね支持の固有値" << std::endl;
    const double k = 100.0;      // N/mm
    const double weight = 500.0; // N (質量は重量/重力加速度)
    const int n = 5;

    auto m = std::make_shared<FEModel>();
    for (int i = 0; i < n; i++)
    {
        m->Nodes.push_back(Node(i, i * 1000.0, 0.0, 0.0));
        for (int d = 0; d < 6; d++)
            m->Nodes[i].Fix.BoundaryTypes[d] = ConstraintType::Fix;
        m->Nodes[i].Fix.BoundaryTypes[2] = ConstraintType::Free; // 鉛直のみばねで支える
        m->Nodes[i].MassData.Mass = weight;
    }
    for (int i = 0; i < n; i++)
        AddSpring(*m, i, 2, k * (i + 1) * (i + 1));

    FEVibrationAnalysis vib(m);
    int nconv = vib.Compute(3);
    Check(nconv == 3, "3モード求まる");
    if (nconv == 3)
    {
        double mass = weight / m->GraityAccel;
        std::vector<double> eigs = vib.EigenValues();
        std::sort(eigs.begin(), eigs.end());
        for (int i = 0; i < 3; i++)
            CheckNear(eigs[i], std::sqrt(k * (i + 1) * (i + 1) / mass), 1e-8,
                      "omega" + std::to_string(i + 1) + " = sqrt(k/m)");
    }
}

// 4) ばね支持と固定支持の混在: 反力の合計が荷重と釣り合う
void TestMixedSupports()
{
    std::cout << "[4] ばね支持と固定支持の混在" << std::endl;
    const double L = 1000.0, E = 205000.0, I = 1.0e7;
    const double k = 50.0;
    const double P = -3000.0;

    auto m = std::make_shared<FEModel>();
    for (int i = 0; i < 3; i++)
        m->Nodes.push_back(Node(i, i * L, 0.0, 0.0));
    m->Materials.push_back(Material(E, 0.3));
    m->Sections.push_back(Section(1000.0, I, I, 1.0e7));

    m->Nodes[0].Fix.FixAll();                 // 左端: 固定
    for (int d : {0, 1, 3, 5})                // 右端: 鉛直のみばね支持
        m->Nodes[2].Fix.BoundaryTypes[d] = ConstraintType::Fix;

    m->add_beam_element(0, 0, 1, 0, 0, 0.0);
    m->add_beam_element(1, 1, 2, 0, 0, 0.0);
    AddSpring(*m, 2, 2, k);

    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.push_back(std::make_shared<NodeLoad>(1, 0.0, 0.0, P));

    FELinearStaticOp op(m, loads);
    op.Compute();

    double sum_pz = 0.0;
    double r_spring = 0.0;
    std::vector<NodeLoad> react_forces = op.GetReactForces();
    for (NodeLoad &nl : react_forces)
    {
        sum_pz += nl.Pz();
        if (nl.id == 2)
            r_spring = nl.Pz();
    }
    CheckNear(sum_pz, -P, 1e-8, "反力の合計が荷重と釣り合う");
    CheckNear(r_spring, -k * op.GetDisplacements()[2].Dz(), 1e-8, "ばね反力 R = -k・u");
}
// 5) 動的解析: 減衰があるばね支持でも、荷重と反力(弾性+減衰)が釣り合うか
//    剛性比例減衰では C = a·K なので、ばね自由度には c = a·k の減衰が付く。
//    反力が -k·u のみだと速度が非零のステップで釣り合いが崩れる。
void TestSpringDynamicReaction()
{
    std::cout << "[5] 減衰のあるばね支持の動的反力" << std::endl;
    const double k = 200.0;      // N/mm
    const double weight = 500.0; // N
    const int n = 5;
    const double P = -1000.0;

    auto m = std::make_shared<FEModel>();
    for (int i = 0; i < n; i++)
    {
        m->Nodes.push_back(Node(i, i * 1000.0, 0.0, 0.0));
        for (int d = 0; d < 6; d++)
            m->Nodes[i].Fix.BoundaryTypes[d] = ConstraintType::Fix;
        m->Nodes[i].Fix.BoundaryTypes[2] = ConstraintType::Free;
        m->Nodes[i].MassData.Mass = weight;
    }
    for (int i = 0; i < n; i++)
        AddSpring(*m, i, 2, k * (i + 1));

    // 節点0に正弦波の鉛直荷重
    const double dt = 0.002;
    const int steps = 200;
    std::vector<NodeLoadData> pattern{NodeLoadData(0, 0.0, 0.0, P)};
    std::vector<double> factors(steps + 1);
    for (int i = 0; i <= steps; i++)
        factors[i] = std::sin(2.0 * PI * (i * dt) / 0.2);
    auto load = std::make_shared<NodalDynamicLoad>(dt, pattern, factors);

    auto damp = std::make_shared<FEDynamicStiffDampInitializer>(0.05);
    DynamicAnalysis da(m, load, damp);
    da.SetTimeGrid(dt, steps);
    Check(da.Initialize(), "初期化できる");

    // 各ステップで「外力 + 反力 + 慣性力 = 0」(鉛直成分)を確認する
    double max_residual = 0.0, max_scale = 0.0;
    for (int s = 0; s < steps; s++)
    {
        da.ComputeStep();
        double t = da.CurrentTime();
        double f_ext = P * std::sin(2.0 * PI * t / 0.2); // 節点0のみに作用

        double r_sum = 0.0;
        std::vector<NodeLoad> react = da.GetReactForces();
        for (NodeLoad &nl : react)
            r_sum += nl.Pz();

        double inertia = 0.0;
        std::vector<Displacement> acc = da.GetAccelerations();
        for (size_t i = 0; i < m->Nodes.size(); i++)
            inertia += -(m->Nodes[i].MassData.SumMass() / m->GraityAccel) * acc[i].Dz();

        max_residual = std::max(max_residual, std::abs(f_ext + r_sum + inertia));
        max_scale = std::max(max_scale, std::abs(r_sum));
    }
    CheckNear(max_residual / max_scale, 0.0, 1e-8, "各ステップで 外力+反力+慣性力 = 0 (相対)");
}

// 加速度一定の応答スペクトル(Sd = Sa/omega^2)
class ConstantAccelSpectrum : public IResponseSpectrum
{
public:
    double Sa;
    explicit ConstantAccelSpectrum(double sa) : Sa(sa) {}
    double Acceleration(double) override { return Sa; }
    double Velocity(double t) override { return Sa * t / (2.0 * PI); }
    double Displacement(double t) override
    {
        double r = t / (2.0 * PI);
        return Sa * r * r;
    }
};

// 6) 応答スペクトル法: ばね支点の反力 R = -k・u
//    3と同じく互いに独立した1質点系を鉛直方向に加振する。各モードは1節点だけが
//    動くため、節点の応答は |u_i| = Sd(T_i) = Sa/omega_i^2 となる。
void TestSpringResponseSpectrum()
{
    std::cout << "[6] 応答スペクトル法でのばね支点の反力" << std::endl;
    const double k = 100.0;      // N/mm
    const double weight = 500.0; // N
    const double Sa = 3000.0;    // mm/s^2
    const int n = 5;

    auto m = std::make_shared<FEModel>();
    for (int i = 0; i < n; i++)
    {
        m->Nodes.push_back(Node(i, i * 1000.0, 0.0, 0.0));
        for (int d = 0; d < 6; d++)
            m->Nodes[i].Fix.BoundaryTypes[d] = ConstraintType::Fix;
        m->Nodes[i].Fix.BoundaryTypes[2] = ConstraintType::Free;
        m->Nodes[i].MassData.Mass = weight;
    }
    for (int i = 0; i < n; i++)
        AddSpring(*m, i, 2, k * (i + 1) * (i + 1));

    FEVibrationAnalysis vib(m);
    int nconv = vib.Compute(3);
    Check(nconv == 3, "3モード求まる");
    if (nconv != 3)
        return;

    ConstantAccelSpectrum spectrum(Sa);
    ResponseSpectrumMethod rsm(m, vib, Vector(0, 0, 1), &spectrum, ResponseSpectrumMethodType::SRSS);
    std::vector<Displacement> disp = rsm.GetDisplacements();
    std::vector<NodeLoad> react = rsm.GetReactForces();

    const double mass = weight / m->GraityAccel;
    for (int i = 0; i < 3; i++)
    {
        const double ki = k * (i + 1) * (i + 1);
        const std::string tag = "節点" + std::to_string(i);
        CheckNear(std::abs(disp[i].Dz()), Sa * mass / ki, 1e-8, tag + " |u| = Sa/omega^2");

        double r = 0.0;
        for (NodeLoad &nl : react)
            if (nl.id == i)
                r = nl.Pz();
        CheckNear(r, -ki * disp[i].Dz(), 1e-10, tag + " ばね反力 R = -k・u");
    }
}

// 7) 固定した自由度にばねを重ねる: ばねは変位0なので力を持たず、固定反力がそのまま残る
//    (ばね反力を上書きで書き込むと、固定反力が0で消えてしまう)
void TestSpringOnFixedDof()
{
    std::cout << "[7] 固定した自由度に重ねたばね" << std::endl;
    const double L = 2000.0, E = 205000.0, I = 1.0e7;
    const double P = -1000.0;

    auto m = std::make_shared<FEModel>();
    m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
    m->Nodes.push_back(Node(1, L, 0.0, 0.0));
    m->Materials.push_back(Material(E, 0.3));
    m->Sections.push_back(Section(1000.0, I, I, 1.0e7));
    m->Nodes[0].Fix.FixAll(); // 片持ち梁の根元
    m->add_beam_element(0, 0, 1, 0, 0, 0.0);
    AddSpring(*m, 0, 2, 500.0); // 固定済みの鉛直自由度に重ねる

    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.push_back(std::make_shared<NodeLoad>(1, 0.0, 0.0, P));

    FELinearStaticOp op(m, loads);
    op.Compute();

    std::vector<NodeLoad> react = op.GetReactForces();
    Check(react.size() == 1, "反力の節点数 = 1");
    if (!react.empty())
        CheckNear(react[0].Pz(), -P, 1e-10, "根元の鉛直反力 R = -P");
}

// 8) 1節点に複数のばね要素: 剛性は和 (k1 + k2)、反力も両方の和で R = -P
void TestMultipleSpringsOnNode()
{
    std::cout << "[8] 1節点に複数のばね要素" << std::endl;
    const double k1 = 150.0, k2 = 350.0; // N/mm
    const double P = -1000.0;            // N

    auto m = std::make_shared<FEModel>();
    m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
    for (int d = 0; d < 6; d++)
        m->Nodes[0].Fix.BoundaryTypes[d] = ConstraintType::Fix;
    m->Nodes[0].Fix.BoundaryTypes[2] = ConstraintType::Free;
    AddSpring(*m, 0, 2, k1);
    AddSpring(*m, 0, 2, k2);

    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.push_back(std::make_shared<NodeLoad>(0, 0.0, 0.0, P));

    FELinearStaticOp op(m, loads);
    op.Compute();

    CheckNear(op.GetDisplacements()[0].Dz(), P / (k1 + k2), 1e-10, "変位 u = P/(k1 + k2)");
    std::vector<NodeLoad> react = op.GetReactForces();
    Check(react.size() == 1, "反力の節点数 = 1");
    if (!react.empty())
        CheckNear(react[0].Pz(), -P, 1e-10, "反力の合計 R = -P");
}
} // namespace

int main()
{
    std::cout << std::setprecision(15); // 実装変更の前後で出力を比較できる桁数
    std::cout << "Spring support test" << std::endl;
    TestSingleSpring();
    TestSpringWithBeam();
    TestSpringVibration();
    TestMixedSupports();
    TestSpringDynamicReaction();
    TestSpringResponseSpectrum();
    TestSpringOnFixedDof();
    TestMultipleSpringsOnNode();

    if (g_failed == 0)
        std::cout << "\nRESULT: PASS (all checks)" << std::endl;
    else
        std::cout << "\nRESULT: FAIL (" << g_failed << " checks)" << std::endl;
    return g_failed == 0 ? 0 : 1;
}

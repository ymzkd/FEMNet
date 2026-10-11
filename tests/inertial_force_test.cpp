// 慣性力(InertialForce)の梁要素への載荷の検証テスト
//   1) 材端ばね・剛域の無い梁の等価節点荷重が、整合質量 × 加速度と一致する
//      (向き・Beta が様々な部材、加速度の各成分)
//   2) 慣性力と、ρA·a を局所座標の等分布荷重(AxialPolyLoad + BeamPolyLoad)にした部材荷重とが、
//      節点変位・反力・部材応力・部材のたわみで一致する(素の梁、材端ばね・剛域のある梁)
//   3) 両端固定の梁: 軸力・せん断力・モーメント・たわみ・反力の閉形式解(局所 x・y・z の各成分)
//   4) 材端ばね・剛域のある梁: 閉形式解
//   5) 荷重組み合わせ(FromCombination)で係数倍した慣性力の応力が係数倍になる
// 引数に出力ファイルを与えると、実装変更の前後で比較するための節点変位・反力・部材応力を書き出す
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Model.h"
#include "Elements/Elements.h"
#include "LoadComponent.h"
#include "FELinearStaticOp.h"

namespace
{
int g_failed = 0;

const double E = 205000.0;  // N/mm2
const double Rho = 7.85e-9; // t/mm3
const double G = 9806.65;   // mm/s2

using Loads = std::vector<std::shared_ptr<LoadBase>>;

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
    if (!(err <= rel_tol))
        g_failed++;
}

// 値の組どうしを、期待値の組の最大絶対値を基準にした誤差で比較する(有限でない値があれば失敗)
void CheckClose(const std::vector<double> &actual, const std::vector<double> &expected,
                double rel_tol, const std::string &what)
{
    double scale = 0.0, err = 0.0;
    bool finite = true;
    for (size_t i = 0; i < expected.size() && i < actual.size(); i++)
    {
        finite = finite && std::isfinite(actual[i]) && std::isfinite(expected[i]);
        scale = std::max(scale, std::abs(expected[i]));
        err = std::max(err, std::abs(actual[i] - expected[i]));
    }
    double rel = finite ? err / (scale > 0.0 ? scale : 1.0) : NAN;
    bool ok = actual.size() == expected.size() && rel <= rel_tol;
    std::cout << (ok ? "  [OK]   " : "  [FAIL] ") << what << " : rel.err=" << rel << std::endl;
    if (!ok)
        g_failed++;
}

void Append(std::vector<double> &dst, const std::vector<double> &src)
{
    dst.insert(dst.end(), src.begin(), src.end());
}

std::vector<double> Translations(const Displacement &d) { return {d.Dx(), d.Dy(), d.Dz()}; }
std::vector<double> Rotations(const Displacement &d) { return {d.Rx(), d.Ry(), d.Rz()}; }
std::vector<double> Forces(const BeamStressData &s) { return {s.Nx, s.Qy, s.Qz}; }
std::vector<double> Moments(const BeamStressData &s) { return {s.Mx, s.My, s.Mz}; }

// 反力を節点番号順に、力とモーメントに分けて並べる
void ReactionValues(FELinearStaticOp &op, std::vector<double> &f, std::vector<double> &m)
{
    std::vector<NodeLoad> react = op.GetReactForces();
    std::sort(react.begin(), react.end(),
              [](const NodeLoad &a, const NodeLoad &b) { return a.id < b.id; });
    for (NodeLoad &nl : react)
    {
        Append(f, {nl.Px(), nl.Py(), nl.Pz()});
        Append(m, {nl.Mx(), nl.My(), nl.Mz()});
    }
}

// 節点荷重を節点ごとに合計する
std::vector<double> SumByNode(const std::vector<NodeLoadData> &loads, size_t node_num)
{
    std::vector<double> f(node_num * 6, 0.0);
    for (const NodeLoadData &nl : loads)
        for (int d = 0; d < 6; d++)
            f[nl.id * 6 + d] += nl.loads[d];
    return f;
}

std::vector<double> Pick(const std::vector<double> &v, int first)
{
    std::vector<double> out;
    for (size_t i = first; i < v.size(); i += 6)
        out.insert(out.end(), v.begin() + i, v.begin() + i + 3);
    return out;
}

// 加速度の向き(成分ごとと、3成分を同時に与えたもの)
struct AccelCase
{
    std::string name;
    double ax, ay, az;
};
const AccelCase Accels[] = {
    {"ax", 0.3 * G, 0.0, 0.0},
    {"ay", 0.0, -0.2 * G, 0.0},
    {"az", 0.0, 0.0, -G},
    {"axyz", 0.3 * G, -0.2 * G, -G},
};

std::shared_ptr<InertialForce> Inertial(const AccelCase &a)
{
    return std::make_shared<InertialForce>(a.ax, a.ay, a.az);
}

// 材端条件(XY面は z軸まわり曲げ、XZ面は y軸まわり曲げ)
struct EndCond
{
    double lzi = 0, lzj = 0, lyi = 0, lyj = 0;
    double syi = 1, bzi = 1, syj = 1, bzj = 1; // XY面の固定度
    double szi = 1, byi = 1, szj = 1, byj = 1; // XZ面の固定度
};

void SetEnds(ComplexBeamElement &e, const EndCond &c)
{
    e.lzi = c.lzi;
    e.lzj = c.lzj;
    e.lyi = c.lyi;
    e.lyj = c.lyj;
    e.Lambda_syi = c.syi;
    e.Lambda_bzi = c.bzi;
    e.Lambda_syj = c.syj;
    e.Lambda_bzj = c.bzj;
    e.Lambda_szi = c.szi;
    e.Lambda_byi = c.byi;
    e.Lambda_szj = c.szj;
    e.Lambda_byj = c.byj;
}

// ---------------------------------------------------------------------------
// 立体骨組: 柱4本(高さの違う柱頭)、柱頭の梁4本(Y方向は傾斜)、頂部の節点へ集まる登り梁4本、
// ブレース(トラス)1本、三角形の板1枚。柱脚は固定。
// complex = true では梁と登り梁を材端ばね・剛域のある ComplexBeamElement にする
// 節点: 0-3 柱脚 / 4-7 柱頭 / 8 頂部。要素: 0-3 柱 / 4-7 梁 / 8-11 登り梁 / 12 ブレース / 13 板
const int FrameBeamIds[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

std::shared_ptr<FEModel> BuildFrame(bool complex)
{
    auto m = std::make_shared<FEModel>();
    const double pts[9][3] = {
        {0.0, 0.0, 0.0}, {6000.0, 0.0, 0.0}, {0.0, 4000.0, 0.0}, {6000.0, 4000.0, 0.0},
        {0.0, 0.0, 3500.0}, {6000.0, 0.0, 3500.0}, {0.0, 4000.0, 3800.0}, {6000.0, 4000.0, 3800.0},
        {2600.0, 2300.0, 5200.0}};
    for (int i = 0; i < 9; i++)
        m->Nodes.push_back(Node(i, pts[i][0], pts[i][1], pts[i][2]));
    for (int i = 0; i < 4; i++)
        m->Nodes[i].Fix.FixAll();

    m->Materials.push_back(Material(E, 0.3, Rho));
    m->Sections.push_back(Section(1.2e4, 3.0e8, 2.0e8, 5.0e7)); // 柱
    m->Sections.push_back(Section(8000.0, 4.0e7, 2.5e8, 1.0e6)); // 梁
    m->Sections.push_back(Section(2000.0, 1.0e6, 1.0e6, 1.0e5)); // ブレース

    m->add_beam_element(0, 0, 4, 0, 0, 0.0);
    m->add_beam_element(1, 1, 5, 0, 0, 0.3);
    m->add_beam_element(2, 2, 6, 0, 0, 0.0);
    m->add_beam_element(3, 3, 7, 0, 0, 1.2);

    struct Member
    {
        int n0, n1;
        double beta;
        EndCond c;
    };
    std::vector<Member> members(8);
    members[0] = {4, 5, 0.0, {}};
    members[0].c.lzi = 300; members[0].c.lzj = 200; members[0].c.lyi = 250; members[0].c.lyj = 400;
    members[0].c.syi = 0.6; members[0].c.bzi = 0.3; members[0].c.bzj = 0.7;
    members[0].c.byi = 0.5; members[0].c.szj = 0.8; members[0].c.byj = 0.01;
    members[1] = {6, 7, 0.2, {}};
    members[1].c.bzi = members[1].c.bzj = members[1].c.byi = members[1].c.byj = 0.01; // 両端ピン
    members[2] = {4, 6, 0.0, {}};
    members[2].c.lzi = members[2].c.lzj = members[2].c.lyi = members[2].c.lyj = 250.0; // 剛域のみ
    members[3] = {5, 7, -0.4, {}};
    members[3].c.bzi = 0.4; members[3].c.byj = 0.6; // ばねのみ
    members[4] = {4, 8, 0.5, {}};
    members[4].c.lzi = members[4].c.lyi = 400.0; members[4].c.bzj = members[4].c.byj = 0.2;
    members[5] = {5, 8, 0.0, {}};
    members[5].c.szi = 0.7; members[5].c.syj = 0.5; members[5].c.lzj = 150.0;
    members[6] = {6, 8, 1.0, {}}; // ばね・剛域なし
    members[7] = {7, 8, -0.8, {}};
    members[7].c.byi = 0.05; members[7].c.lyj = 300.0;

    for (size_t k = 0; k < members.size(); k++)
    {
        const int id = 4 + (int)k;
        const Member &mb = members[k];
        if (complex)
        {
            ComplexBeamElement e(id, &m->Nodes[mb.n0], &m->Nodes[mb.n1], &m->Sections[1], m->Materials[0], mb.beta);
            SetEnds(e, mb.c);
            m->add_element(e);
        }
        else
            m->add_beam_element(id, mb.n0, mb.n1, 1, 0, mb.beta);
    }
    m->add_truss_element(12, 0, 5, 2, 0);
    m->add_tri_plate_element(13, 4, 5, 7, 150.0, 0);
    return m;
}

// 慣性力と同じ荷重を、梁は ρA·a を局所座標の各軸に射影した等分布の部材荷重で、
// それ以外の要素は要素の等価節点荷重で与える。
// 0 の成分は作らない(AxialPolyLoad は荷重が 0 だと重心位置が 0/0 になる)
Loads EquivalentLoads(FEModel &m, const AccelCase &a)
{
    Loads loads;
    const Eigen::Vector3d acc(a.ax, a.ay, a.az);
    const std::vector<double> whole{0.0, 1.0};
    for (auto &ep : m.Elements)
    {
        if (BeamElement *e = dynamic_cast<BeamElement *>(ep.get()))
        {
            const double rhoA = e->Sec->A * e->Mat.dense;
            auto w = [&](const Vector &v) { return rhoA * (acc.x() * v.x + acc.y() * v.y + acc.z() * v.z); };
            const double wx = w(e->XAxis()), wy = w(e->YAxis()), wz = w(e->ZAxis());
            if (wx != 0.0)
                loads.push_back(std::make_shared<AxialPolyLoad>(std::vector<double>{wx, wx}, whole, e));
            if (wy != 0.0)
                loads.push_back(std::make_shared<BeamPolyLoad>(std::vector<double>{wy, wy}, whole, e, BeamLoadAxis::YAxis));
            if (wz != 0.0)
                loads.push_back(std::make_shared<BeamPolyLoad>(std::vector<double>{wz, wz}, whole, e, BeamLoadAxis::ZAxis));
        }
        else
        {
            for (NodeLoadData &nl : ep->InertialForceToNodeLoadData(acc))
                loads.push_back(std::make_shared<NodeLoad>(nl.id, nl.Px(), nl.Py(), nl.Pz(), nl.Mx(), nl.My(), nl.Mz()));
        }
    }
    return loads;
}

const double StressPoints[] = {0.0, 0.03, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.97, 1.0};

// ---------------------------------------------------------------------------
// 1) 材端ばね・剛域の無い梁の等価節点荷重 = 整合質量 × 加速度
void TestNodeLoadsMatchConsistentMass()
{
    std::cout << "[1] 素の梁の等価節点荷重と、整合質量 × 加速度の一致" << std::endl;
    auto m = BuildFrame(false);
    for (const AccelCase &a : Accels)
    {
        std::vector<NodeLoadData> expected;
        for (auto &ep : m->Elements)
        {
            std::vector<NodeLoadData> el = ep->InertialForceToNodeLoadData(Eigen::Vector3d(a.ax, a.ay, a.az));
            expected.insert(expected.end(), el.begin(), el.end());
        }
        std::vector<double> fa = SumByNode(m->InnertialForceToNodeLoads(InertialForce(a.ax, a.ay, a.az)), m->Nodes.size());
        std::vector<double> fe = SumByNode(expected, m->Nodes.size());
        CheckClose(Pick(fa, 0), Pick(fe, 0), 1e-12, a.name + " 節点の力");
        CheckClose(Pick(fa, 3), Pick(fe, 3), 1e-12, a.name + " 節点のモーメント");
    }
}

// ---------------------------------------------------------------------------
// 2) 慣性力と、等分布の部材荷重にしたものとの一致
void TestEquivalentMemberLoads()
{
    for (bool complex : {false, true})
    {
        std::cout << "[2] 慣性力と等分布の部材荷重の一致(" << (complex ? "材端ばね・剛域あり" : "素の梁") << ")" << std::endl;
        for (const AccelCase &a : Accels)
        {
            auto m1 = BuildFrame(complex);
            auto m2 = BuildFrame(complex);
            FELinearStaticOp o1(m1, Loads{Inertial(a)}), o2(m2, EquivalentLoads(*m2, a));
            o1.Compute();
            o2.Compute();

            std::vector<double> t1, t2, r1, r2;
            for (size_t n = 0; n < m1->Nodes.size(); n++)
            {
                Append(t1, Translations(o1.GetDisplacements()[n]));
                Append(t2, Translations(o2.GetDisplacements()[n]));
                Append(r1, Rotations(o1.GetDisplacements()[n]));
                Append(r2, Rotations(o2.GetDisplacements()[n]));
            }
            CheckClose(t1, t2, 1e-9, a.name + " 節点の変位");
            CheckClose(r1, r2, 1e-9, a.name + " 節点の回転");

            std::vector<double> f1, f2, mm1, mm2;
            ReactionValues(o1, f1, mm1);
            ReactionValues(o2, f2, mm2);
            CheckClose(f1, f2, 1e-9, a.name + " 反力");
            CheckClose(mm1, mm2, 1e-9, a.name + " 反力モーメント");

            std::vector<double> sf1, sf2, sm1, sm2, d1, d2;
            for (int eid : FrameBeamIds)
                for (double p : StressPoints)
                {
                    Append(sf1, Forces(o1.GetBeamStress(eid, p)));
                    Append(sf2, Forces(o2.GetBeamStress(eid, p)));
                    Append(sm1, Moments(o1.GetBeamStress(eid, p)));
                    Append(sm2, Moments(o2.GetBeamStress(eid, p)));
                    Append(d1, Translations(o1.GetBeamDisplace(eid, p)));
                    Append(d2, Translations(o2.GetBeamDisplace(eid, p)));
                }
            CheckClose(sf1, sf2, 1e-9, a.name + " 部材の軸力・せん断力");
            CheckClose(sm1, sm2, 1e-9, a.name + " 部材のモーメント");
            CheckClose(d1, d2, 1e-9, a.name + " 部材のたわみ");
        }
    }
}

// ---------------------------------------------------------------------------
// 全体X軸の向きで Beta = 0 の梁1本(局所座標 = 全体座標)。両端の節点を固定し、
// 荷重と関係しない j端のねじれだけ自由にする
const double L = 6000.0, A = 8000.0, Iy = 4.0e7, Iz = 2.5e8;
const AccelCase Accel3{"axyz", 0.3 * G, -0.2 * G, -G};

std::shared_ptr<FEModel> BuildSingleBeam()
{
    auto m = std::make_shared<FEModel>();
    m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
    m->Nodes.push_back(Node(1, L, 0.0, 0.0));
    m->Materials.push_back(Material(E, 0.3, Rho));
    m->Sections.push_back(Section(A, Iy, Iz, 1.0e6));
    m->Nodes[0].Fix.FixAll();
    m->Nodes[1].Fix.FixAll();
    m->Nodes[1].Fix.BoundaryTypes[3] = ConstraintType::Free;
    return m;
}

// 3) 両端固定の梁の閉形式解
void TestFixedBeamClosedForm()
{
    std::cout << "[3] 両端固定の梁の閉形式解" << std::endl;
    auto m = BuildSingleBeam();
    m->add_beam_element(0, 0, 1, 0, 0, 0.0);

    FELinearStaticOp op(m, Loads{Inertial(Accel3)});
    op.Compute();

    const double wx = Rho * A * Accel3.ax, wy = Rho * A * Accel3.ay, wz = Rho * A * Accel3.az;
    auto S = [&](double p) { return op.GetBeamStress(0, p); };

    CheckNear(S(0.0).Nx, wx * L / 2.0, 1e-10, "i端の軸力");
    CheckNear(S(0.25).Nx, wx * L / 4.0, 1e-10, "1/4点の軸力");
    CheckNear(S(1.0).Nx, -wx * L / 2.0, 1e-10, "j端の軸力");

    CheckNear(S(0.0).Qy, wy * L / 2.0, 1e-10, "XY面 i端せん断力");
    CheckNear(S(1.0).Qy, -wy * L / 2.0, 1e-10, "XY面 j端せん断力");
    CheckNear(S(0.0).Mz, -wy * L * L / 12.0, 1e-10, "XY面 i端モーメント");
    CheckNear(S(0.5).Mz, wy * L * L / 24.0, 1e-10, "XY面 中央モーメント");
    CheckNear(S(1.0).Mz, -wy * L * L / 12.0, 1e-10, "XY面 j端モーメント");

    CheckNear(S(0.0).Qz, wz * L / 2.0, 1e-10, "XZ面 i端せん断力");
    CheckNear(S(1.0).Qz, -wz * L / 2.0, 1e-10, "XZ面 j端せん断力");
    CheckNear(S(0.0).My, -wz * L * L / 12.0, 1e-10, "XZ面 i端モーメント");
    CheckNear(S(0.5).My, wz * L * L / 24.0, 1e-10, "XZ面 中央モーメント");
    CheckNear(S(1.0).My, -wz * L * L / 12.0, 1e-10, "XZ面 j端モーメント");

    Displacement dm = op.GetBeamDisplace(0, 0.5);
    CheckNear(dm.Dy(), wy * std::pow(L, 4) / (384.0 * E * Iz), 1e-9, "XY面 中央のたわみ");
    CheckNear(dm.Dz(), wz * std::pow(L, 4) / (384.0 * E * Iy), 1e-9, "XZ面 中央のたわみ");

    // 節点0の反力(XZ面のモーメントは θy の向きで符号が逆)
    for (NodeLoad &r : op.GetReactForces())
    {
        if (r.id != 0)
            continue;
        CheckNear(r.Px(), -wx * L / 2.0, 1e-10, "節点0の反力 X");
        CheckNear(r.Py(), -wy * L / 2.0, 1e-10, "節点0の反力 Y");
        CheckNear(r.Pz(), -wz * L / 2.0, 1e-10, "節点0の反力 Z");
        CheckNear(r.Mz(), -wy * L * L / 12.0, 1e-10, "節点0の反力モーメント Z");
        CheckNear(r.My(), wz * L * L / 12.0, 1e-10, "節点0の反力モーメント Y");
    }
}

// ---------------------------------------------------------------------------
// 4) 材端ばね・剛域のある梁の閉形式解
void TestComplexBeamClosedForm()
{
    // 4a) 両端ピン相当の回転ばね。面ごとに固定度を変え、面の取り違えがあれば結果が変わるようにする
    //     端部モーメント: wL^2/12 / (1 + (1-λ)/(2λ))、中央はそこから wL^2/8 だけ変わる
    {
        std::cout << "[4a] 材端に回転ばねのある梁の閉形式解" << std::endl;
        const double lam_z = 0.01, lam_y = 0.3;
        auto m = BuildSingleBeam();
        ComplexBeamElement beam(0, &m->Nodes[0], &m->Nodes[1], &m->Sections[0], m->Materials[0], 0.0);
        beam.Lambda_bzi = beam.Lambda_bzj = lam_z;
        beam.Lambda_byi = beam.Lambda_byj = lam_y;
        m->add_element(beam);

        FELinearStaticOp op(m, Loads{Inertial(Accel3)});
        op.Compute();

        const double wx = Rho * A * Accel3.ax, wy = Rho * A * Accel3.ay, wz = Rho * A * Accel3.az;
        const double Mez = -wy * L * L / 12.0 / (1.0 + (1.0 - lam_z) / (2.0 * lam_z));
        const double Mey = -wz * L * L / 12.0 / (1.0 + (1.0 - lam_y) / (2.0 * lam_y));
        auto S = [&](double p) { return op.GetBeamStress(0, p); };

        CheckNear(S(0.0).Nx, wx * L / 2.0, 1e-10, "i端の軸力");
        CheckNear(S(0.0).Mz, Mez, 1e-10, "XY面 i端モーメント");
        CheckNear(S(0.5).Mz, Mez + wy * L * L / 8.0, 1e-10, "XY面 中央モーメント");
        CheckNear(S(1.0).Mz, Mez, 1e-10, "XY面 j端モーメント");
        CheckNear(S(0.0).Qy, wy * L / 2.0, 1e-10, "XY面 i端せん断力");
        CheckNear(S(0.0).My, Mey, 1e-10, "XZ面 i端モーメント");
        CheckNear(S(0.5).My, Mey + wz * L * L / 8.0, 1e-10, "XZ面 中央モーメント");
        CheckNear(S(1.0).My, Mey, 1e-10, "XZ面 j端モーメント");
        CheckNear(S(0.0).Qz, wz * L / 2.0, 1e-10, "XZ面 i端せん断力");
    }

    // 4b) 剛域のある両端固定の梁。内法区間は長さ ln の両端固定梁になり、
    //     剛域上の荷重はそのまま節点へ伝わる。軸方向は剛域を考えず全長で抵抗する
    {
        std::cout << "[4b] 剛域のある両端固定の梁の閉形式解" << std::endl;
        const double lzi = 512.0, lzj = 1024.0, lyi = 768.0, lyj = 256.0;
        auto m = BuildSingleBeam();
        ComplexBeamElement beam(0, &m->Nodes[0], &m->Nodes[1], &m->Sections[0], m->Materials[0], 0.0);
        beam.lzi = lzi;
        beam.lzj = lzj;
        beam.lyi = lyi;
        beam.lyj = lyj;
        m->add_element(beam);

        FELinearStaticOp op(m, Loads{Inertial(Accel3)});
        op.Compute();

        const double wx = Rho * A * Accel3.ax, wy = Rho * A * Accel3.ay, wz = Rho * A * Accel3.az;
        const double lnz = L - lzi - lzj, lny = L - lyi - lyj;
        auto S = [&](double x) { return op.GetBeamStress(0, x / L); };

        CheckNear(S(0.0).Nx, wx * L / 2.0, 1e-10, "i端の軸力");
        CheckNear(S(lzi).Mz, -wy * lnz * lnz / 12.0, 1e-10, "XY面 i側の面のモーメント");
        CheckNear(S(lzi + lnz / 2.0).Mz, wy * lnz * lnz / 24.0, 1e-10, "XY面 内法中央のモーメント");
        CheckNear(S(L - lzj).Mz, -wy * lnz * lnz / 12.0, 1e-10, "XY面 j側の面のモーメント");
        CheckNear(S(lzi).Qy, wy * lnz / 2.0, 1e-10, "XY面 i側の面のせん断力");
        CheckNear(S(lyi).My, -wz * lny * lny / 12.0, 1e-10, "XZ面 i側の面のモーメント");
        CheckNear(S(lyi + lny / 2.0).My, wz * lny * lny / 24.0, 1e-10, "XZ面 内法中央のモーメント");
        CheckNear(S(L - lyj).My, -wz * lny * lny / 12.0, 1e-10, "XZ面 j側の面のモーメント");
        CheckNear(S(lyi).Qz, wz * lny / 2.0, 1e-10, "XZ面 i側の面のせん断力");

        for (NodeLoad &r : op.GetReactForces())
        {
            if (r.id != 0)
                continue;
            CheckNear(r.Px(), -wx * L / 2.0, 1e-10, "節点0の反力 X");
            CheckNear(r.Py(), -wy * (lzi + lnz / 2.0), 1e-10, "節点0の反力 Y");
            CheckNear(r.Pz(), -wz * (lyi + lny / 2.0), 1e-10, "節点0の反力 Z");
        }
    }
}

// ---------------------------------------------------------------------------
// 5) 荷重組み合わせ: FromCombination で係数倍した慣性力の応力が係数倍になる
void TestCombination()
{
    std::cout << "[5] 荷重組み合わせで係数倍した慣性力" << std::endl;
    const double factor = 1.7;
    auto m = BuildFrame(true);
    auto base = std::make_shared<FELinearStaticOp>(m, Loads{Inertial(Accel3)});
    base->Compute();
    auto comb = FELinearStaticOp::FromCombination(m, {LinearStaticDeformFactor(base, factor)});
    comb->Compute();

    std::vector<double> sb, sc, mb, mc;
    for (int eid : FrameBeamIds)
        for (double p : StressPoints)
        {
            for (double v : Forces(base->GetBeamStress(eid, p)))
                sb.push_back(factor * v);
            Append(sc, Forces(comb->GetBeamStress(eid, p)));
            for (double v : Moments(base->GetBeamStress(eid, p)))
                mb.push_back(factor * v);
            Append(mc, Moments(comb->GetBeamStress(eid, p)));
        }
    CheckClose(sc, sb, 1e-12, "部材の軸力・せん断力");
    CheckClose(mc, mb, 1e-12, "部材のモーメント");
}

// ---------------------------------------------------------------------------
// 実装変更の前後で比較するための値を書き出す
void Dump(const std::string &path)
{
    std::ofstream ofs(path);
    ofs << std::setprecision(17);
    for (bool complex : {false, true})
    {
        const std::string tag = complex ? "complex" : "plain";
        for (const AccelCase &a : Accels)
        {
            auto m = BuildFrame(complex);
            const std::string head = tag + " " + a.name + " ";

            std::vector<double> nl = SumByNode(m->InnertialForceToNodeLoads(InertialForce(a.ax, a.ay, a.az)), m->Nodes.size());
            for (size_t n = 0; n < m->Nodes.size(); n++)
            {
                ofs << head << "nodeload " << n;
                for (int d = 0; d < 6; d++)
                    ofs << " " << nl[n * 6 + d];
                ofs << "\n";
            }

            FELinearStaticOp op(m, Loads{Inertial(a)});
            op.Compute();
            std::vector<Displacement> disp = op.GetDisplacements();
            for (size_t n = 0; n < disp.size(); n++)
            {
                ofs << head << "disp " << n;
                for (int d = 0; d < 6; d++)
                    ofs << " " << disp[n].displace[d];
                ofs << "\n";
            }
            std::vector<NodeLoad> react = op.GetReactForces();
            std::sort(react.begin(), react.end(),
                      [](const NodeLoad &x, const NodeLoad &y) { return x.id < y.id; });
            for (NodeLoad &r : react)
                ofs << head << "react " << r.id << " " << r.Px() << " " << r.Py() << " " << r.Pz()
                    << " " << r.Mx() << " " << r.My() << " " << r.Mz() << "\n";
            for (int eid : FrameBeamIds)
                for (double p : StressPoints)
                {
                    BeamStressData s = op.GetBeamStress(eid, p);
                    ofs << head << "stress " << eid << "@" << p << " " << s.Nx << " " << s.Qy << " " << s.Qz
                        << " " << s.Mx << " " << s.My << " " << s.Mz << "\n";
                    Displacement d = op.GetBeamDisplace(eid, p);
                    ofs << head << "bdisp " << eid << "@" << p << " " << d.Dx() << " " << d.Dy() << " " << d.Dz() << "\n";
                }
        }
    }
}
} // namespace

int main(int argc, char **argv)
{
    std::cout << std::setprecision(15);
    std::cout << "InertialForce on beam elements test" << std::endl;
    TestNodeLoadsMatchConsistentMass();
    TestEquivalentMemberLoads();
    TestFixedBeamClosedForm();
    TestComplexBeamClosedForm();
    TestCombination();

    if (argc > 1)
    {
        Dump(argv[1]);
        std::cout << "\nwrote " << argv[1] << std::endl;
    }

    if (g_failed == 0)
        std::cout << "\nRESULT: PASS (all checks)" << std::endl;
    else
        std::cout << "\nRESULT: FAIL (" << g_failed << " checks)" << std::endl;
    return g_failed == 0 ? 0 : 1;
}

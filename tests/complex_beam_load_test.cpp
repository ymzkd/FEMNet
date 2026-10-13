// 材端ばね・剛域のある梁(ComplexBeamElement)にかかる部材荷重(BeamPolyLoad)の検証テスト
//   1) ばね・剛域の無い ComplexBeamElement は BeamElement と同じ結果になる
//   2) 両端固定の節点に両端ピン(λ=0.01)の梁: 閉形式の端部・中央モーメントとたわみ(XY面・XZ面)
//   3) 縮小した節点荷重が、要素の変位変換(Tu・Tb)の転置で整合節点力を変換したものと一致する
//   4) 門形骨組: 梁を荷重の無い区間で切り、ばね・剛域を荷重の無い要素に持たせた参照モデルと一致する
//   5) 剛域上の荷重は、剛体として節点へ移した節点荷重と同じ結果になる
//   6) 面をまたぐ荷重は、面で切り分けた荷重の和と同じ結果になる(段差のある荷重を含む)
//   7) 節点荷重の合計と節点iまわりのモーメントが、荷重の合力と釣り合う
//   8) 荷重0・逆対称の区間と段差を含む軸方向荷重(AxialPolyLoad)が有限で、合力・軸力分布と一致する
#include <algorithm>
#include <cmath>
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

const double E = 205000.0; // N/mm2

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

// 値の組どうしを、期待値の組の最大絶対値を基準にした誤差で比較する
void CheckClose(const std::vector<double> &actual, const std::vector<double> &expected,
                double rel_tol, const std::string &what)
{
    double scale = 0.0, err = 0.0;
    for (size_t i = 0; i < expected.size() && i < actual.size(); i++)
    {
        scale = std::max(scale, std::abs(expected[i]));
        err = std::max(err, std::abs(actual[i] - expected[i]));
    }
    double rel = err / (scale > 0.0 ? scale : 1.0);
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

std::shared_ptr<BeamPolyLoad> Poly(BeamElement *e, BeamLoadAxis axis,
                                   std::vector<double> w, std::vector<double> params)
{
    return std::make_shared<BeamPolyLoad>(w, params, e, axis);
}

// 1次分布の荷重(x0 で w0、x1 で w1)の合力 P と、原点まわりの1次モーメント Px
void Resultant(double w0, double w1, double x0, double x1, double &P, double &Px)
{
    double h = x1 - x0;
    P = (w0 + w1) * h / 2.0;
    Px = h / 6.0 * (w0 * (2.0 * x0 + x1) + w1 * (x0 + 2.0 * x1));
}

// 局所座標の節点力 [Q_i, M_i, Q_j, M_j] を節点荷重から取り出す(θ = dv/dx の向き)。
// 要素が全体X軸の向きで Beta = 0 のときだけ使う(局所座標 = 全体座標)
Eigen::Vector4d LocalNodeForces(LoadBase &load, BendingPlane plane)
{
    std::vector<NodeLoadData> nl = load.NodeLoads();
    if (plane == BendingPlane::XY)
        return Eigen::Vector4d(nl[0].Py(), nl[0].Mz(), nl[1].Py(), nl[1].Mz());
    return Eigen::Vector4d(nl[0].Pz(), -nl[0].My(), nl[1].Pz(), -nl[1].My());
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

// 固定度 λ を、基準剛性 D に対するばね剛性に換算する
double Spring(double lambda, double D) { return lambda / (1.0 - lambda) * D; }

// ---------------------------------------------------------------------------
// 1) ばね・剛域の無い ComplexBeamElement と BeamElement の一致
void TestPlainEquivalence()
{
    std::cout << "[1] ばね・剛域の無い ComplexBeamElement と BeamElement の一致" << std::endl;

    auto build = [](bool complex) {
        auto m = std::make_shared<FEModel>();
        m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
        m->Nodes.push_back(Node(1, 3000.0, 2000.0, 1500.0));
        m->Materials.push_back(Material(E, 0.3));
        m->Sections.push_back(Section(8000.0, 4.0e7, 2.0e8, 1.0e6));
        m->Nodes[0].Fix.FixAll();
        for (int d : {0, 1, 2})
            m->Nodes[1].Fix.BoundaryTypes[d] = ConstraintType::Fix;
        if (complex)
            m->add_element(ComplexBeamElement(0, &m->Nodes[0], &m->Nodes[1], &m->Sections[0], m->Materials[0], 0.4));
        else
            m->add_element(BeamElement(0, &m->Nodes[0], &m->Nodes[1], &m->Sections[0], m->Materials[0], 0.4));
        return m;
    };
    auto make_loads = [](std::shared_ptr<FEModel> m) {
        BeamElement *e = m->GetBeamElement(0);
        return Loads{Poly(e, BeamLoadAxis::YAxis, {-8.0, -12.0, -4.0}, {0.0, 0.4, 1.0}),
                     Poly(e, BeamLoadAxis::ZAxis, {3.0, 3.0}, {0.2, 0.7})};
    };

    auto mb = build(false), mc = build(true);
    Loads lb = make_loads(mb), lc = make_loads(mc);

    bool same = true;
    for (size_t k = 0; k < lb.size(); k++)
    {
        std::vector<NodeLoadData> nb = lb[k]->NodeLoads(), nc = lc[k]->NodeLoads();
        for (size_t n = 0; n < nb.size(); n++)
            for (int d = 0; d < 6; d++)
                same = same && nb[n].loads[d] == nc[n].loads[d];
    }
    Check(same, "節点荷重がビット単位で一致する");

    FELinearStaticOp ob(mb, lb), oc(mc, lc);
    ob.Compute();
    oc.Compute();

    CheckClose(Rotations(oc.GetDisplacements()[1]), Rotations(ob.GetDisplacements()[1]), 1e-9, "節点の回転");

    std::vector<double> fb, fc, mmb, mmc, db, dc;
    for (double p : {0.0, 0.1, 0.25, 0.4, 0.5, 0.7, 0.85, 1.0})
    {
        Append(fb, Forces(ob.GetBeamStress(0, p)));
        Append(fc, Forces(oc.GetBeamStress(0, p)));
        Append(mmb, Moments(ob.GetBeamStress(0, p)));
        Append(mmc, Moments(oc.GetBeamStress(0, p)));
        Append(db, Translations(ob.GetBeamDisplace(0, p)));
        Append(dc, Translations(oc.GetBeamDisplace(0, p)));
    }
    CheckClose(fc, fb, 1e-9, "部材の軸力・せん断力");
    CheckClose(mmc, mmb, 1e-9, "部材のモーメント");
    CheckClose(dc, db, 1e-9, "部材のたわみ");
}

// ---------------------------------------------------------------------------
// 2) 両端固定の節点に両端ピン(λ=0.01)の梁、等分布荷重
//    端部モーメント: wL^2/12 / (1 + (1-λ)/(2λ))、中央はそこから wL^2/8 だけ変わる
void TestPinnedBeamClosedForm()
{
    const double L = 6000.0, w = -10.0, lam = 0.01;
    const double Iy = 4.0e7, Iz = 2.5e8;

    for (BendingPlane plane : {BendingPlane::XY, BendingPlane::XZ})
    {
        const bool xy = plane == BendingPlane::XY;
        std::cout << "[2] 両端固定の節点に両端ピンの梁(" << (xy ? "XY面" : "XZ面") << ")" << std::endl;

        auto m = std::make_shared<FEModel>();
        m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
        m->Nodes.push_back(Node(1, L, 0.0, 0.0));
        m->Materials.push_back(Material(E, 0.3));
        m->Sections.push_back(Section(8000.0, Iy, Iz, 1.0e6));
        m->Nodes[0].Fix.FixAll();
        m->Nodes[1].Fix.FixAll();
        // 自由度が0だと解けないので、荷重と関係しない j端の軸方向だけ自由にする
        m->Nodes[1].Fix.BoundaryTypes[0] = ConstraintType::Free;

        ComplexBeamElement beam(0, &m->Nodes[0], &m->Nodes[1], &m->Sections[0], m->Materials[0], 0.0);
        // 荷重の無い側の面にもばねを入れ、面の取り違えがあれば結果が変わるようにする
        if (xy)
        {
            beam.Lambda_bzi = beam.Lambda_bzj = lam;
            beam.Lambda_byi = beam.Lambda_byj = 0.5;
        }
        else
        {
            beam.Lambda_byi = beam.Lambda_byj = lam;
            beam.Lambda_bzi = beam.Lambda_bzj = 0.5;
        }
        m->add_element(beam);

        BeamLoadAxis axis = xy ? BeamLoadAxis::YAxis : BeamLoadAxis::ZAxis;
        Loads loads{Poly(m->GetBeamElement(0), axis, {w, w}, {0.0, 1.0})};
        FELinearStaticOp op(m, loads);
        op.Compute();

        const double EI = E * (xy ? Iz : Iy);
        const double Me = -w * L * L / 12.0 / (1.0 + (1.0 - lam) / (2.0 * lam));
        auto M = [&](double p) { BeamStressData s = op.GetBeamStress(0, p); return xy ? s.Mz : s.My; };
        auto Q = [&](double p) { BeamStressData s = op.GetBeamStress(0, p); return xy ? s.Qy : s.Qz; };

        CheckNear(M(0.0), Me, 1e-10, "i端モーメント");
        CheckNear(M(1.0), Me, 1e-10, "j端モーメント");
        CheckNear(M(0.5), Me + w * L * L / 8.0, 1e-10, "中央モーメント");
        CheckNear(Q(0.0), w * L / 2.0, 1e-10, "i端せん断力");

        // Mz = -EI v'' から、単純梁のたわみに端部モーメントによる分を足す
        Displacement dm = op.GetBeamDisplace(0, 0.5);
        double v_mid = 5.0 * w * std::pow(L, 4) / (384.0 * EI) + Me * L * L / (8.0 * EI);
        CheckNear(xy ? dm.Dy() : dm.Dz(), v_mid, 1e-9, "中央のたわみ");

        // 節点0の反力: 鉛直は -wL/2、モーメントは端部モーメント(XZ面は θy の向きで符号が逆)
        for (NodeLoad &r : op.GetReactForces())
        {
            if (r.id != 0)
                continue;
            CheckNear(xy ? r.Py() : r.Pz(), -w * L / 2.0, 1e-10, "節点0の反力");
            CheckNear(xy ? r.Mz() : r.My(), xy ? Me : -Me, 1e-10, "節点0の反力モーメント");
        }
    }
}

// ---------------------------------------------------------------------------
// 3) 縮小した節点荷重 = (Tu・Tb)^T f_b
//    Tu・Tb(節点の変位から内側の梁の端の変位への変換)は、要素の DisplaceAt を面の位置で
//    評価して取り出す。f_b は内法区間を両端固定梁とみなした整合節点力。
//    剛域長は位置 p が2進で割り切れる値にして、面の位置ちょうどで評価する。
void TestReductionMatchesDisplacementTransform()
{
    std::cout << "[3] 縮小した節点荷重と変位変換の転置の一致" << std::endl;
    const double L = 8192.0;

    auto m = std::make_shared<FEModel>();
    m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
    m->Nodes.push_back(Node(1, L, 0.0, 0.0));
    m->Materials.push_back(Material(E, 0.3));
    m->Sections.push_back(Section(8000.0, 4.0e7, 2.5e8, 1.0e6));

    EndCond c;
    c.lzi = 512; c.lzj = 1024; c.lyi = 1024; c.lyj = 512;
    c.syi = 0.6; c.bzi = 0.3; c.syj = 0.9; c.bzj = 0.05;
    c.szi = 0.8; c.byi = 0.02; c.szj = 0.4; c.byj = 0.7;
    ComplexBeamElement beam(0, &m->Nodes[0], &m->Nodes[1], &m->Sections[0], m->Materials[0], 0.0);
    SetEnds(beam, c);
    m->add_element(beam);
    BeamElement *e = m->GetBeamElement(0);

    for (BendingPlane plane : {BendingPlane::XY, BendingPlane::XZ})
    {
        const bool xy = plane == BendingPlane::XY;
        const double li = xy ? c.lzi : c.lyi, lj = xy ? c.lzj : c.lyj;
        const double ln = L - li - lj;

        // 単位の節点変位を与えて、面の位置の内側の梁の端の変位を読む(XZ面は ψ = -θy)
        const int dofs[4] = {xy ? 1 : 2, xy ? 5 : 4, xy ? 1 : 2, xy ? 5 : 4};
        const double sgn[4] = {1.0, xy ? 1.0 : -1.0, 1.0, xy ? 1.0 : -1.0};
        Eigen::Matrix4d G;
        for (int k = 0; k < 4; k++)
        {
            Displacement d[2];
            d[k / 2].displace[dofs[k]] = sgn[k];
            Displacement di = e->DisplaceAt(d[0], d[1], li / L);
            Displacement dj = e->DisplaceAt(d[0], d[1], 1.0 - lj / L);
            G(0, k) = di.displace[dofs[0]];
            G(1, k) = sgn[1] * di.displace[dofs[1]];
            G(2, k) = dj.displace[dofs[2]];
            G(3, k) = sgn[3] * dj.displace[dofs[3]];
        }

        // 内法区間の中だけにかかる荷重
        const std::vector<double> w = xy ? std::vector<double>{-6.0, -15.0, -3.0} : std::vector<double>{4.0, 9.0, 2.0};
        const std::vector<double> params = xy ? std::vector<double>{0.1, 0.5, 0.8} : std::vector<double>{0.15, 0.6, 0.85};
        auto load = Poly(e, xy ? BeamLoadAxis::YAxis : BeamLoadAxis::ZAxis, w, params);

        Eigen::Vector4d f_b = Eigen::Vector4d::Zero();
        for (size_t k = 0; k + 1 < w.size(); k++)
        {
            double a = params[k] * L - li, b = params[k + 1] * L - li;
            BeamTrapezoidalLoad t(w[k], w[k + 1], a, b - a, ln - b);
            f_b += Eigen::Vector4d(t.R0(), t.M0(), t.RA(), -t.MA());
        }

        Eigen::Vector4d expected = G.transpose() * f_b;
        Eigen::Vector4d actual = LocalNodeForces(*load, plane);
        std::string tag = xy ? "XY面" : "XZ面";
        CheckClose({actual(0), actual(2)}, {expected(0), expected(2)}, 1e-10, tag + " 節点の力");
        CheckClose({actual(1), actual(3)}, {expected(1), expected(3)}, 1e-10, tag + " 節点のモーメント");
    }
}

// ---------------------------------------------------------------------------
// 門形骨組: 柱脚を固定した柱2本と、柱頭をつなぐ ComplexBeamElement の梁
// 節点: 0,1 柱脚 / 2 梁のi端 / 3 梁のj端。要素: 0,1 柱 / 2 梁
const double FL = 6000.0, FH = 3500.0, FBeta = 0.3;
const double BeamIy = 4.0e7, BeamIz = 2.5e8;

EndCond FrameEnds()
{
    EndCond c;
    c.lzi = 300; c.lzj = 200; c.lyi = 250; c.lyj = 400;
    c.syi = 0.6; c.bzi = 0.3; c.syj = 1.0; c.bzj = 0.7;  // XY面: i端にせん断ばねと回転ばね
    c.szi = 1.0; c.byi = 0.5; c.szj = 0.8; c.byj = 0.01; // XZ面: j端はほぼピン
    return c;
}

void AddFrameBase(FEModel &m)
{
    m.Materials.push_back(Material(E, 0.3));
    m.Sections.push_back(Section(1.2e4, 3.0e8, 2.0e8, 5.0e7));      // 柱
    m.Sections.push_back(Section(8000.0, BeamIy, BeamIz, 1.0e6));    // 梁
}

std::shared_ptr<FEModel> BuildFrame(const EndCond &c)
{
    auto m = std::make_shared<FEModel>();
    m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
    m->Nodes.push_back(Node(1, FL, 0.0, 0.0));
    m->Nodes.push_back(Node(2, 0.0, 0.0, FH));
    m->Nodes.push_back(Node(3, FL, 0.0, FH));
    AddFrameBase(*m);
    m->Nodes[0].Fix.FixAll();
    m->Nodes[1].Fix.FixAll();
    m->add_beam_element(0, 0, 2, 0, 0, 0.0);
    m->add_beam_element(1, 1, 3, 0, 0, 0.0);
    ComplexBeamElement beam(2, &m->Nodes[2], &m->Nodes[3], &m->Sections[1], m->Materials[0], FBeta);
    SetEnds(beam, c);
    m->add_element(beam);
    return m;
}

// 4) の参照モデル: 梁を荷重の無い区間 [0, a] と [L-b, L] で切り、材端ばね・剛域はその区間の
// 荷重の無い ComplexBeamElement に、荷重は中央の BeamElement に持たせる。材端ばねは、
// 元の要素の内法長を基準にした固定度を剛性値に換算して与える。
// 節点: 0,1 柱脚 / 2 i端 / 3 j端 / 4 x=a / 5 x=L-b。要素: 0,1 柱 / 2 i側 / 3 中央 / 4 j側
std::shared_ptr<FEModel> BuildFrameReference(const EndCond &c, double a, double b)
{
    auto m = std::make_shared<FEModel>();
    m->Nodes.push_back(Node(0, 0.0, 0.0, 0.0));
    m->Nodes.push_back(Node(1, FL, 0.0, 0.0));
    m->Nodes.push_back(Node(2, 0.0, 0.0, FH));
    m->Nodes.push_back(Node(3, FL, 0.0, FH));
    m->Nodes.push_back(Node(4, a, 0.0, FH));
    m->Nodes.push_back(Node(5, FL - b, 0.0, FH));
    AddFrameBase(*m);
    m->Nodes[0].Fix.FixAll();
    m->Nodes[1].Fix.FixAll();
    m->add_beam_element(0, 0, 2, 0, 0, 0.0);
    m->add_beam_element(1, 1, 3, 0, 0, 0.0);

    const double lz = FL - c.lzi - c.lzj, ly = FL - c.lyi - c.lyj;
    const double EIz = E * BeamIz, EIy = E * BeamIy;
    const double Dsz = 12.0 * EIz / (lz * lz * lz), Dbz = 4.0 * EIz / lz;
    const double Dsy = 12.0 * EIy / (ly * ly * ly), Dby = 4.0 * EIy / ly;

    // set_k* は剛域を除いた長さで固定度に換算するので、剛域を先に設定する
    ComplexBeamElement si(2, &m->Nodes[2], &m->Nodes[4], &m->Sections[1], m->Materials[0], FBeta);
    si.lzi = c.lzi;
    si.lyi = c.lyi;
    if (c.syi < 1) si.set_ksyi(Spring(c.syi, Dsz));
    if (c.bzi < 1) si.set_kbzi(Spring(c.bzi, Dbz));
    if (c.szi < 1) si.set_kszi(Spring(c.szi, Dsy));
    if (c.byi < 1) si.set_kbyi(Spring(c.byi, Dby));
    m->add_element(si);

    m->add_beam_element(3, 4, 5, 1, 0, FBeta);

    ComplexBeamElement sj(4, &m->Nodes[5], &m->Nodes[3], &m->Sections[1], m->Materials[0], FBeta);
    sj.lzj = c.lzj;
    sj.lyj = c.lyj;
    if (c.syj < 1) sj.set_ksyj(Spring(c.syj, Dsz));
    if (c.bzj < 1) sj.set_kbzj(Spring(c.bzj, Dbz));
    if (c.szj < 1) sj.set_kszj(Spring(c.szj, Dsy));
    if (c.byj < 1) sj.set_kbyj(Spring(c.byj, Dby));
    m->add_element(sj);
    return m;
}

// 4) 門形骨組の参照モデルとの比較
void TestFrameAgainstReference()
{
    std::cout << "[4] 門形骨組: 荷重の無い区間で分割した参照モデルとの比較" << std::endl;
    const double a = 1000.0, b = 1000.0, lm = FL - a - b;
    EndCond c = FrameEnds();

    auto mm = BuildFrame(c);
    auto mr = BuildFrameReference(c, a, b);

    // 荷重は [a, L-b] の中だけ。柱頭には水平力も与えて節点を動かす
    Loads lmm{Poly(mm->GetBeamElement(2), BeamLoadAxis::YAxis, {-6.0, -14.0, -9.0},
                   {1000.0 / FL, 2500.0 / FL, 4800.0 / FL}),
              Poly(mm->GetBeamElement(2), BeamLoadAxis::ZAxis, {4.0, 7.0}, {1500.0 / FL, 4000.0 / FL}),
              std::make_shared<NodeLoad>(2, 20000.0, 5000.0, 0.0)};
    Loads lmr{Poly(mr->GetBeamElement(3), BeamLoadAxis::YAxis, {-6.0, -14.0, -9.0},
                   {0.0, 1500.0 / lm, 3800.0 / lm}),
              Poly(mr->GetBeamElement(3), BeamLoadAxis::ZAxis, {4.0, 7.0}, {500.0 / lm, 3000.0 / lm}),
              std::make_shared<NodeLoad>(2, 20000.0, 5000.0, 0.0)};

    FELinearStaticOp om(mm, lmm), orf(mr, lmr);
    om.Compute();
    orf.Compute();

    std::vector<double> tm, tr, rm, rr;
    for (int n : {2, 3})
    {
        Append(tm, Translations(om.GetDisplacements()[n]));
        Append(tr, Translations(orf.GetDisplacements()[n]));
        Append(rm, Rotations(om.GetDisplacements()[n]));
        Append(rr, Rotations(orf.GetDisplacements()[n]));
    }
    CheckClose(tm, tr, 1e-9, "梁端の節点の変位");
    CheckClose(rm, rr, 1e-9, "梁端の節点の回転");

    std::vector<double> fm, fr, mmm, mmr;
    ReactionValues(om, fm, mmm);
    ReactionValues(orf, fr, mmr);
    CheckClose(fm, fr, 1e-9, "柱脚の反力");
    CheckClose(mmm, mmr, 1e-9, "柱脚の反力モーメント");

    // 参照モデルで、梁の節点iからの距離 x にある要素と要素内の位置
    auto locate = [&](double x) -> std::pair<int, double> {
        if (x <= a) return {2, x / a};
        if (x <= FL - b) return {3, (x - a) / lm};
        return {4, (x - (FL - b)) / b};
    };

    std::vector<double> sfm, sfr, smm, smr;
    for (double x : {0.0, 150.0, 250.0, 280.0, 300.0, 600.0, 1000.0, 2000.0, 3000.0,
                     4500.0, 5000.0, 5500.0, 5600.0, 5800.0, 5900.0, 6000.0})
    {
        auto [eid, p] = locate(x);
        Append(sfm, Forces(om.GetBeamStress(2, x / FL)));
        Append(sfr, Forces(orf.GetBeamStress(eid, p)));
        Append(smm, Moments(om.GetBeamStress(2, x / FL)));
        Append(smr, Moments(orf.GetBeamStress(eid, p)));
    }
    CheckClose(sfm, sfr, 1e-9, "梁の軸力・せん断力");
    CheckClose(smm, smr, 1e-9, "梁のモーメント");

    // たわみ(並進のみ。部材荷重による回転は GetBeamDisplace に含まれない)
    std::vector<double> dm, dr;
    for (double x : {0.0, 150.0, 280.0, 600.0, 1000.0, 1800.0, 2500.0, 3300.0, 4000.0,
                     4700.0, 5000.0, 5500.0, 5700.0, 5900.0, 6000.0})
    {
        auto [eid, p] = locate(x);
        Append(dm, Translations(om.GetBeamDisplace(2, x / FL)));
        Append(dr, Translations(orf.GetBeamDisplace(eid, p)));
    }
    CheckClose(dm, dr, 1e-9, "梁のたわみ");
}

// ---------------------------------------------------------------------------
// 剛域上の1次分布荷重を、剛体として節点へ移した節点荷重
std::shared_ptr<NodeLoad> RigidTransfer(BeamElement *e, BeamLoadAxis axis,
                                        double w0, double w1, double p0, double p1, int end)
{
    const double L = e->length();
    double P, Px;
    Resultant(w0, w1, p0 * L, p1 * L, P, Px);
    const double moment = end == 0 ? Px : Px - L * P; // 節点まわりの1次モーメント
    Vector dir = axis == BeamLoadAxis::YAxis ? e->YAxis() : e->ZAxis();
    // XY面は局所z軸まわり、XZ面は (x 軸) × (z 軸) = -y 軸まわり
    Vector mv = axis == BeamLoadAxis::YAxis ? e->ZAxis() * moment : e->YAxis() * (-moment);
    return std::make_shared<NodeLoad>(e->Nodes[end]->id, dir.x * P, dir.y * P, dir.z * P, mv.x, mv.y, mv.z);
}

// 5) 剛域上の荷重と、剛体として移した節点荷重の比較
void TestRigidZoneLoad()
{
    std::cout << "[5] 剛域上の荷重と、剛体として節点へ移した節点荷重の一致" << std::endl;
    EndCond c = FrameEnds();
    auto m1 = BuildFrame(c);
    auto m2 = BuildFrame(c);
    BeamElement *e1 = m1->GetBeamElement(2);
    BeamElement *e2 = m2->GetBeamElement(2);

    // 剛域: XY面 i 300 / j 200、XZ面 i 250 / j 400
    struct Piece { BeamLoadAxis axis; double w0, w1, x0, x1; int end; };
    const Piece pieces[] = {
        {BeamLoadAxis::YAxis, -20.0, -10.0, 0.0, 200.0, 0},
        {BeamLoadAxis::YAxis, -5.0, -15.0, 5850.0, 6000.0, 1},
        {BeamLoadAxis::ZAxis, 8.0, 8.0, 0.0, 240.0, 0},
        {BeamLoadAxis::ZAxis, 3.0, 9.0, 5700.0, 5900.0, 1},
    };
    Loads l1, l2;
    for (const Piece &pc : pieces)
    {
        l1.push_back(Poly(e1, pc.axis, {pc.w0, pc.w1}, {pc.x0 / FL, pc.x1 / FL}));
        l2.push_back(RigidTransfer(e2, pc.axis, pc.w0, pc.w1, pc.x0 / FL, pc.x1 / FL, pc.end));
    }

    FELinearStaticOp o1(m1, l1), o2(m2, l2);
    o1.Compute();
    o2.Compute();

    std::vector<double> t1, t2, r1, r2;
    for (int n : {2, 3})
    {
        Append(t1, Translations(o1.GetDisplacements()[n]));
        Append(t2, Translations(o2.GetDisplacements()[n]));
        Append(r1, Rotations(o1.GetDisplacements()[n]));
        Append(r2, Rotations(o2.GetDisplacements()[n]));
    }
    CheckClose(t1, t2, 1e-10, "梁端の節点の変位");
    CheckClose(r1, r2, 1e-10, "梁端の節点の回転");

    // 両方の面の内法区間 [300, 5600] では、部材荷重による応力の項が打ち消し合う
    std::vector<double> f1, f2, mm1, mm2, d1, d2;
    for (double x : {300.0, 400.0, 1000.0, 3000.0, 5000.0, 5600.0})
    {
        Append(f1, Forces(o1.GetBeamStress(2, x / FL)));
        Append(f2, Forces(o2.GetBeamStress(2, x / FL)));
        Append(mm1, Moments(o1.GetBeamStress(2, x / FL)));
        Append(mm2, Moments(o2.GetBeamStress(2, x / FL)));
        Append(d1, Translations(o1.GetBeamDisplace(2, x / FL)));
        Append(d2, Translations(o2.GetBeamDisplace(2, x / FL)));
    }
    CheckClose(f1, f2, 1e-10, "内法区間の軸力・せん断力");
    CheckClose(mm1, mm2, 1e-10, "内法区間のモーメント");
    CheckClose(d1, d2, 1e-10, "内法区間のたわみ");
}

// ---------------------------------------------------------------------------
// 多角形荷重を、各区間と切り分け位置 cuts(要素長に対する比)で切った2点の荷重の列にする
Loads SplitPoly(BeamElement *e, BeamLoadAxis axis, const std::vector<double> &w,
                const std::vector<double> &params, const std::vector<double> &cuts)
{
    Loads out;
    for (size_t k = 0; k + 1 < w.size(); k++)
    {
        const double a = params[k], b = params[k + 1];
        if (b <= a)
            continue;
        std::vector<double> pts{a};
        for (double cp : cuts)
            if (cp > a && cp < b)
                pts.push_back(cp);
        pts.push_back(b);
        auto w_at = [&](double p) { return w[k] + (w[k + 1] - w[k]) * (p - a) / (b - a); };
        for (size_t i = 0; i + 1 < pts.size(); i++)
            out.push_back(Poly(e, axis, {w_at(pts[i]), w_at(pts[i + 1])}, {pts[i], pts[i + 1]}));
    }
    return out;
}

// 6) と 7) の荷重: 剛域と内法区間の境目(面)をまたぎ、0.2 の位置に段差がある
const std::vector<double> SplitWy{-4.0, -10.0, -2.0, -12.0, -6.0};
const std::vector<double> SplitPy{0.0, 0.2, 0.2, 0.9, 1.0};
const std::vector<double> SplitWz{5.0, 2.0, 6.0};
const std::vector<double> SplitPz{0.0, 0.5, 1.0};

// 6) 面をまたぐ荷重と、面で切り分けた荷重の和の比較
void TestCrossingLoadSuperposition()
{
    std::cout << "[6] 面をまたぐ荷重と、面で切り分けた荷重の和の一致" << std::endl;
    EndCond c = FrameEnds();
    auto m1 = BuildFrame(c);
    auto m2 = BuildFrame(c);
    BeamElement *e1 = m1->GetBeamElement(2);
    BeamElement *e2 = m2->GetBeamElement(2);

    Loads l1{Poly(e1, BeamLoadAxis::YAxis, SplitWy, SplitPy),
             Poly(e1, BeamLoadAxis::ZAxis, SplitWz, SplitPz)};
    Loads l2 = SplitPoly(e2, BeamLoadAxis::YAxis, SplitWy, SplitPy, {c.lzi / FL, 1.0 - c.lzj / FL});
    Loads l2z = SplitPoly(e2, BeamLoadAxis::ZAxis, SplitWz, SplitPz, {c.lyi / FL, 1.0 - c.lyj / FL});
    l2.insert(l2.end(), l2z.begin(), l2z.end());

    // 節点荷重の合計
    auto total = [](Loads &loads) {
        std::vector<double> f(12, 0.0);
        for (auto &l : loads)
        {
            std::vector<NodeLoadData> nl = l->NodeLoads();
            for (size_t n = 0; n < nl.size(); n++)
                for (int d = 0; d < 6; d++)
                    f[(nl[n].id == 2 ? 0 : 6) + d] += nl[n].loads[d];
        }
        return f;
    };
    CheckClose(total(l2), total(l1), 1e-12, "節点荷重の合計");

    FELinearStaticOp o1(m1, l1), o2(m2, l2);
    o1.Compute();
    o2.Compute();

    std::vector<double> f1, f2, mm1, mm2, d1, d2;
    bool finite = true;
    for (double p : {0.0, 0.03, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 0.9, 0.95, 0.97, 1.0})
    {
        BeamStressData s1 = o1.GetBeamStress(2, p);
        Append(f1, Forces(s1));
        Append(f2, Forces(o2.GetBeamStress(2, p)));
        Append(mm1, Moments(s1));
        Append(mm2, Moments(o2.GetBeamStress(2, p)));
        Append(d1, Translations(o1.GetBeamDisplace(2, p)));
        Append(d2, Translations(o2.GetBeamDisplace(2, p)));
        for (double v : Moments(s1))
            finite = finite && std::isfinite(v);
    }
    Check(finite, "段差の位置でも応力が有限");
    CheckClose(f1, f2, 1e-10, "梁の軸力・せん断力");
    CheckClose(mm1, mm2, 1e-10, "梁のモーメント");
    CheckClose(d1, d2, 1e-10, "梁のたわみ");
}

// ---------------------------------------------------------------------------
// 7) 節点荷重の合計と節点iまわりのモーメントが、荷重の合力と釣り合う
void TestNodeLoadEquilibrium()
{
    std::cout << "[7] 節点荷重と荷重の合力の釣り合い" << std::endl;
    auto m = BuildFrame(FrameEnds());
    BeamElement *e = m->GetBeamElement(2);
    const double L = e->length();
    const Point pi = e->Nodes[0]->Location;

    auto vec = [](const Vector &v) { return Eigen::Vector3d(v.x, v.y, v.z); };
    const Eigen::Vector3d ex = vec(e->XAxis());

    struct Case { BeamLoadAxis axis; std::vector<double> w, params; };
    const Case cases[] = {{BeamLoadAxis::YAxis, SplitWy, SplitPy}, {BeamLoadAxis::ZAxis, SplitWz, SplitPz}};
    for (const Case &cs : cases)
    {
        const bool y = cs.axis == BeamLoadAxis::YAxis;
        const Eigen::Vector3d dir = vec(y ? e->YAxis() : e->ZAxis());

        // 荷重の合力と節点iまわりのモーメント
        double P = 0.0, Px = 0.0;
        for (size_t k = 0; k + 1 < cs.w.size(); k++)
        {
            double dP, dPx;
            Resultant(cs.w[k], cs.w[k + 1], cs.params[k] * L, cs.params[k + 1] * L, dP, dPx);
            P += dP;
            Px += dPx;
        }
        const Eigen::Vector3d f_exp = P * dir;
        const Eigen::Vector3d m_exp = Px * ex.cross(dir);

        Eigen::Vector3d f = Eigen::Vector3d::Zero(), mom = Eigen::Vector3d::Zero();
        for (NodeLoadData &nl : Poly(e, cs.axis, cs.w, cs.params)->NodeLoads())
        {
            const Point pn = m->Nodes[nl.id].Location;
            Eigen::Vector3d fn(nl.Px(), nl.Py(), nl.Pz()), mn(nl.Mx(), nl.My(), nl.Mz());
            f += fn;
            mom += Eigen::Vector3d(pn.x - pi.x, pn.y - pi.y, pn.z - pi.z).cross(fn) + mn;
        }
        std::string tag = y ? "YAxis荷重" : "ZAxis荷重";
        CheckClose({f.x(), f.y(), f.z()}, {f_exp.x(), f_exp.y(), f_exp.z()}, 1e-12, tag + " 力の合計");
        CheckClose({mom.x(), mom.y(), mom.z()}, {m_exp.x(), m_exp.y(), m_exp.z()}, 1e-12, tag + " 節点iまわりのモーメント");
    }
}

// ---------------------------------------------------------------------------
// 8) 荷重0の区間・逆対称の台形(w1 + w2 = 0)・段差を含む軸方向荷重(AxialPolyLoad)が、
//    有限で荷重の合力と釣り合い、軸力が N(x) = N_i - ∫w の分布になる
void TestAxialLoadZeroSum()
{
    std::cout << "[8] 荷重0・逆対称の区間と段差を含む軸方向荷重" << std::endl;
    auto m = BuildFrame(FrameEnds());
    BeamElement *e = m->GetBeamElement(2);
    const double L = e->length();
    const Eigen::Vector3d ex(e->XAxis().x, e->XAxis().y, e->XAxis().z);

    struct Case { std::string tag; std::vector<double> w, params; };
    const Case cases[] = {
        // [0,0.3] 荷重0、0.3 で段差、[0.6,1] は -2→2 の逆対称
        {"荷重0・段差・逆対称", {0.0, 0.0, 4.0, 4.0, -2.0, 2.0}, {0.0, 0.3, 0.3, 0.6, 0.6, 1.0}},
        {"全長で荷重0", {0.0, 0.0}, {0.0, 1.0}},
    };
    for (const Case &cs : cases)
    {
        double P = 0.0, Px = 0.0;
        for (size_t k = 0; k + 1 < cs.w.size(); k++)
        {
            double dP, dPx;
            Resultant(cs.w[k], cs.w[k + 1], cs.params[k] * L, cs.params[k + 1] * L, dP, dPx);
            P += dP;
            Px += dPx;
        }

        AxialPolyLoad load(cs.w, cs.params, e);
        std::vector<NodeLoadData> nl = load.NodeLoads();
        Eigen::Vector3d fi(nl[0].Px(), nl[0].Py(), nl[0].Pz()), fj(nl[1].Px(), nl[1].Py(), nl[1].Pz());
        Check(fi.allFinite() && fj.allFinite(), cs.tag + " 節点荷重が有限");
        Eigen::Vector3d f = fi + fj, f_exp = P * ex;
        CheckClose({f.x(), f.y(), f.z()}, {f_exp.x(), f_exp.y(), f_exp.z()}, 1e-12, cs.tag + " 力の合計");
        // 軸方向の荷重なので、節点iまわりの釣り合いは j端の分担 × L = Px
        CheckClose({fj.dot(ex) * L}, {Px}, 1e-12, cs.tag + " 節点iまわりの1次モーメント");

        // 軸力 N(x) = N_i - (0~x の荷重の合計)
        std::vector<double> n, n_exp;
        for (double p : {0.0, 0.1, 0.3, 0.45, 0.6, 0.8, 1.0})
        {
            const double x = p * L;
            double acc = 0.0;
            for (size_t k = 0; k + 1 < cs.w.size(); k++)
            {
                double x0 = cs.params[k] * L, x1 = cs.params[k + 1] * L;
                if (x1 <= x0 || x <= x0)
                    continue;
                double xe = std::min(x, x1);
                double we = cs.w[k] + (cs.w[k + 1] - cs.w[k]) * (xe - x0) / (x1 - x0);
                double dP, dPx;
                Resultant(cs.w[k], we, x0, xe, dP, dPx);
                acc += dP;
            }
            n.push_back(load.GetBeamStress(p).Nx);
            n_exp.push_back(fi.dot(ex) - acc);
        }
        bool finite = std::all_of(n.begin(), n.end(), [](double v) { return std::isfinite(v); });
        Check(finite, cs.tag + " 段差の位置でも軸力が有限");
        CheckClose(n, n_exp, 1e-12, cs.tag + " 軸力の分布");
    }
}
} // namespace

int main()
{
    std::cout << std::setprecision(15); // 実装変更の前後で出力を比較できる桁数
    std::cout << "ComplexBeamElement member load test" << std::endl;
    TestPlainEquivalence();
    TestPinnedBeamClosedForm();
    TestReductionMatchesDisplacementTransform();
    TestFrameAgainstReference();
    TestRigidZoneLoad();
    TestCrossingLoadSuperposition();
    TestNodeLoadEquilibrium();
    TestAxialLoadZeroSum();

    if (g_failed == 0)
        std::cout << "\nRESULT: PASS (all checks)" << std::endl;
    else
        std::cout << "\nRESULT: FAIL (" << g_failed << " checks)" << std::endl;
    return g_failed == 0 ? 0 : 1;
}

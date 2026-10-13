#include <numeric>

#include "LoadComponent.h"

LinePolyLoad::LinePolyLoad(const std::vector<double>& w, const std::vector<double>& params,
                           Node* node_i, Node* node_j, const Vector& direction)
    : node_i(node_i), node_j(node_j), w(w), params(params), direction(direction)
{}

std::vector<NodeLoadData> LinePolyLoad::NodeLoads()
{
    // 公開フィールド・節点位置の現在の値から計算する。
    const double length = (node_j->Location - node_i->Location).norm();
    const Vector unit_direction = direction * (1.0 / direction.norm());

    double force_i = 0.0, force_j = 0.0;
    for (size_t k = 0; k + 1 < w.size(); ++k)
    {
        const double t0 = params[k], t1 = params[k + 1];
        if (t1 == t0) // 荷重の段差。長さ0なので積分への寄与は0。
            continue;

        // s = length * t、Ni = 1 - t、Nj = t による台形荷重の厳密積分。
        const double coeff = length * (t1 - t0) / 6.0;
        force_i += coeff * (w[k] * (2.0 * (1.0 - t0) + (1.0 - t1)) +
                            w[k + 1] * ((1.0 - t0) + 2.0 * (1.0 - t1)));
        force_j += coeff * (w[k] * (2.0 * t0 + t1) + w[k + 1] * (t0 + 2.0 * t1));
    }

    const Vector load_i = force_i * unit_direction;
    const Vector load_j = force_j * unit_direction;
    return {NodeLoadData(node_i->id, load_i.x, load_i.y, load_i.z),
            NodeLoadData(node_j->id, load_j.x, load_j.y, load_j.z)};
}

std::shared_ptr<LoadBase> LinePolyLoad::scaled(double factor) const
{
    std::vector<double> ws = w;
    for (double& value : ws)
        value *= factor;
    return std::make_shared<LinePolyLoad>(ws, params, node_i, node_j, direction);
}


double BeamPolyLoad::R0() {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.R0();
        });
}

double BeamPolyLoad::RA() {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.RA();
        });
}

double BeamPolyLoad::M0() {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.M0();
        });
}

double BeamPolyLoad::MA() {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.MA();
        });
}

double BeamPolyLoad::shear_force(double x) {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [&x](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.shear_force(x);
        });
}

double BeamPolyLoad::bending_moment(double x) {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [&x](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.bending_moment(x);
        });
}

double BeamPolyLoad::deflection(double x, double EI) {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [&x, &EI](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.deflection(x, EI);
        });
}

// 1次分布の荷重(x0 で w0、x1 で w1)の合力 P と、原点まわりの1次モーメント Px
static void linear_load_resultant(double w0, double w1, double x0, double x1, double &P, double &Px)
{
    double h = x1 - x0;
    P = (w0 + w1) * h / 2.0;
    Px = h / 6.0 * (w0 * (2.0 * x0 + x1) + w1 * (x0 + 2.0 * x1));
}

BendingPlane BeamPolyLoad::plane() const
{
    return axis == BeamLoadAxis::ZAxis ? BendingPlane::XZ : BendingPlane::XY;
}

BeamPolyLoad::ReducedLoad BeamPolyLoad::reduce(const BeamBendingEnds &ends)
{
    ReducedLoad r;
    if (ends.IsPlain())
    {
        r.f_node << R0(), M0(), RA(), -MA();
        r.u_b0.setZero();
        r.inner_traps = traps;
        return r;
    }

    const double L = element->length();
    const double x0 = ends.li;     // 内法区間の始点
    const double x1 = L - ends.lj; // 内法区間の終点

    // 剛域上の荷重は剛体として節点へ伝える
    Eigen::Vector4d f_rz = Eigen::Vector4d::Zero();
    for (size_t k = 0; k + 1 < w.size(); k++)
    {
        const double a = params[k] * L;
        const double b = params[k + 1] * L;
        if (b <= a)
            continue;
        const double wa = w[k], wb = w[k + 1];
        auto w_at = [&](double x) { return wa + (wb - wa) * (x - a) / (b - a); };

        // 内法区間の部分: 切り口の荷重強さは線形補間し、内法長を基準にした台形荷重にする
        double c0 = std::max(a, x0), c1 = std::min(b, x1);
        if (c1 > c0)
            r.inner_traps.push_back(BeamTrapezoidalLoad(w_at(c0), w_at(c1), c0 - x0, c1 - c0, x1 - c1));

        double P, Px;
        // 剛域 i の部分
        c0 = a;
        c1 = std::min(b, x0);
        if (c1 > c0)
        {
            linear_load_resultant(w_at(c0), w_at(c1), c0, c1, P, Px);
            f_rz(0) += P;
            f_rz(1) += Px;
        }
        // 剛域 j の部分(節点 j まわりのモーメントは腕が x - L)
        c0 = std::max(a, x1);
        c1 = b;
        if (c1 > c0)
        {
            linear_load_resultant(w_at(c0), w_at(c1), c0, c1, P, Px);
            f_rz(2) += P;
            f_rz(3) += Px - L * P;
        }
    }

    // 内法区間の両端固定梁としての整合節点力
    Eigen::Vector4d f_b = Eigen::Vector4d::Zero();
    for (BeamTrapezoidalLoad &t : r.inner_traps)
        f_b += Eigen::Vector4d(t.R0(), t.M0(), t.RA(), -t.MA());

    ends.ReduceBeamLoad(f_b, r.f_node, r.u_b0);
    r.f_node += f_rz;
    return r;
}

NodeLoadData BeamPolyLoad::node_load(int end, const Eigen::Vector4d &f_node)
{
    Eigen::Matrix3d tr = trans_matrix3(
        element->Nodes[0]->Location, element->Nodes[1]->Location, element->Beta).transpose();

    const double q = f_node(2 * end);
    const double m = f_node(2 * end + 1);
    Eigen::Vector3d mv, qv;
    if (axis == BeamLoadAxis::ZAxis) {
        // XZ面は θ = dw/dx = -θy として計算しているので、モーメントの符号を反転する
        mv = tr * Eigen::Vector3d(0, -m, 0);
        qv = tr * Eigen::Vector3d(0, 0, q);
    }
    else // YAxis
    {
        mv = tr * Eigen::Vector3d(0, 0, m);
        qv = tr * Eigen::Vector3d(0, q, 0);
    }

    return NodeLoadData(element->Nodes[end]->id, qv[0], qv[1], qv[2], mv[0], mv[1], mv[2]);
}

NodeLoadData BeamPolyLoad::load_i()
{
    return node_load(0, reduce(element->BendingEnds(plane())).f_node);
}

NodeLoadData BeamPolyLoad::load_j()
{
    return node_load(1, reduce(element->BendingEnds(plane())).f_node);
}

std::vector<NodeLoadData> BeamPolyLoad::NodeLoads()
{
    Eigen::Vector4d f = reduce(element->BendingEnds(plane())).f_node;
    return { node_load(0, f), node_load(1, f) };
}

BeamStressData BeamPolyLoad::GetBeamStress(double p)
{
    const double x = p * element->length();
    double m = bending_moment(x);
    double q = shear_force(x);

    // 材端ばね・剛域があるときは、i端の端部力を両端固定の値から縮小した値に置き換える。
    // 節点iの側を切り出した自由体の釣り合いから、分布は i端の力の差分だけ線形に変わる
    BeamBendingEnds ends = element->BendingEnds(plane());
    if (!ends.IsPlain())
    {
        Eigen::Vector4d f = reduce(ends).f_node;
        const double dq = f(0) - R0();
        const double dm = f(1) - M0();
        m += dq * x - dm;
        q += dq;
    }

    BeamStressData bsd;
    if (axis == BeamLoadAxis::ZAxis) {
        bsd.My = m;
        bsd.Qz = q;
    }
    else {
        bsd.Mz = m;
        bsd.Qy = q;
    }
    return bsd;
}

Displacement BeamPolyLoad::GetDisplacement(double p)
{
    const double x = p * element->length();
    BeamBendingEnds ends = element->BendingEnds(plane());

    double d = 0;
    if (ends.IsPlain())
        d = deflection(x, ends.EI);
    else
    {
        // 節点を固定したときのたわみ。剛域は節点と一体なので動かない。
        // 内法区間は、材端ばねの変形による分と、内法長の両端固定梁のたわみの和
        const double xi = x - ends.li;
        const double ln = ends.l_inner;
        if (xi >= 0 && xi <= ln)
        {
            ReducedLoad r = reduce(ends);
            const double s = xi / ln, s2 = s * s, s3 = s2 * s;
            const Eigen::Vector4d N(1 - 3 * s2 + 2 * s3, ln * (s - 2 * s2 + s3), 3 * s2 - 2 * s3, ln * (-s2 + s3));
            d = N.dot(r.u_b0);
            for (BeamTrapezoidalLoad &t : r.inner_traps)
                d += t.deflection(xi, ends.EI);
        }
    }

    Eigen::Matrix3d tr = trans_matrix3(
        element->Nodes[0]->Location,
        element->Nodes[1]->Location, element->Beta);

    int vid = (axis == BeamLoadAxis::ZAxis) ? 2 : 1;
    Eigen::Vector3d dt = tr.row(vid) * d;
    return Displacement(dt[0], dt[1], dt[2], 0, 0, 0);
}

double AxialPolyLoad::axial_force(double x)
{
    return std::accumulate(traps.begin(), traps.end(), 0.0, [&x](double sum, AxialTrapezoidalLoad trap) {
        return sum + trap.axial_force(x);
        });
}
double AxialPolyLoad::N0()
{
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, AxialTrapezoidalLoad trap) {
        return sum + trap.n1();
        });
}
double AxialPolyLoad::N1()
{
    return -std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, AxialTrapezoidalLoad trap) {
        return sum + trap.n2();
        });
}

NodeLoadData AxialPolyLoad::load_i()
{
    Vector dir = element->Nodes[1]->Location - element->Nodes[0]->Location;
    dir = dir * (N0() / dir.norm());
    return NodeLoadData(element->Nodes[0]->id, dir.x, dir.y, dir.z);
}

NodeLoadData AxialPolyLoad::load_j()
{
    Vector dir = element->Nodes[1]->Location - element->Nodes[0]->Location;
    dir = dir * (N1() / dir.norm());
    return NodeLoadData(element->Nodes[1]->id, dir.x, dir.y, dir.z);
}

BeamStressData AxialPolyLoad::GetBeamStress(double p)
{
    BeamStressData bsd;
    bsd.Nx = axial_force(p * element->length());
    return bsd;
}

Displacement AxialPolyLoad::GetDisplacement(double p)
{
    return Displacement();
}

std::vector<std::shared_ptr<BeamLoadBase>> InertialForce::ToBeamLoads(BeamElement* element) const
{
    // 加速度を局所座標の成分にして、単位長さあたりの質量 ρA を掛ける
    Eigen::Matrix3d tr = trans_matrix3(
        element->Nodes[0]->Location, element->Nodes[1]->Location, element->Beta);
    Eigen::Vector3d w = element->Sec->A * element->Mat.dense * (tr * Eigen::Vector3d(accels.x, accels.y, accels.z));

    // 0 の成分は作らない
    const std::vector<double> whole{0.0, 1.0};
    std::vector<std::shared_ptr<BeamLoadBase>> loads;
    if (w.x() != 0.0)
        loads.push_back(std::make_shared<AxialPolyLoad>(std::vector<double>{w.x(), w.x()}, whole, element));
    if (w.y() != 0.0)
        loads.push_back(std::make_shared<BeamPolyLoad>(std::vector<double>{w.y(), w.y()}, whole, element, BeamLoadAxis::YAxis));
    if (w.z() != 0.0)
        loads.push_back(std::make_shared<BeamPolyLoad>(std::vector<double>{w.z(), w.z()}, whole, element, BeamLoadAxis::ZAxis));
    return loads;
}

NodeLoadData::NodeLoadData(int _id, double px, double py, double pz)
{
    id = _id;
    loads[0] = px;
    loads[1] = py;
    loads[2] = pz;
    loads[3] = 0;
    loads[4] = 0;
    loads[5] = 0;
}

NodeLoadData::NodeLoadData(int _id, double px, double py, double pz, double mx, double my, double mz)
{
    id = _id;
    loads[0] = px;
    loads[1] = py;
    loads[2] = pz;
    loads[3] = mx;
    loads[4] = my;
    loads[5] = mz;
}

NodeLoad::NodeLoad(int _id, double px, double py, double pz)
{
    id = _id;
    data = NodeLoadData(_id, px, py, pz);
}

NodeLoad::NodeLoad(int _id, double px, double py, double pz, double mx, double my, double mz)
{
    id = _id;
    data = NodeLoadData(_id, px, py, pz, mx, my, mz);
}

std::vector<NodeLoadData> PlateLoad::NodeLoads()
{
	return element->AreaForceToNodeLoadData(load_vecs);
    //return std::vector<NodeLoadData>();
}

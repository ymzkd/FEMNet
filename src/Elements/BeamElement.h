#ifndef _BEAM_ELEMENT_
#define _BEAM_ELEMENT_

#include "BarElement.h"

// 梁の曲げの面
//   XY: 局所y方向のたわみ v と z軸まわりの曲げ Mz
//   XZ: 局所z方向のたわみ w と y軸まわりの曲げ My
enum class BendingPlane
{
    XY,
    XZ
};

// 曲げの面ごとの材端条件
// 要素を 節点i →(剛域 li)→ 材端ばね → 内側の梁(長さ l_inner) → 材端ばね →(剛域 lj)→ 節点j
// の直列とみなす。固定度 lambda は内側の梁の剛性の対角成分(せん断 12EI/l^3、回転 4EI/l)を
// 基準にした値で、1 が剛、0 が完全な解放。
struct BeamBendingEnds
{
    double li = 0, lj = 0;
    double lambda_si = 1, lambda_bi = 1, lambda_sj = 1, lambda_bj = 1;
    double EI = 0;
    double l_inner = 0;

    // 剛域も材端ばねも無い(通常の梁と同じ)
    bool IsPlain() const
    {
        return li == 0 && lj == 0 &&
               lambda_si == 1 && lambda_bi == 1 && lambda_sj == 1 && lambda_bj == 1;
    }

    // 内側の梁にかかる荷重の整合節点力 f_b を材端ばねで静的縮小し、剛域を通して節点へ移す。
    // ベクトルの並びは [v_i, θ_i, v_j, θ_j](θ = dv/dx)。XZ面は θ を -θy に読み替えて使う。
    //   f_node: 節点を固定したときに節点位置へ伝わる等価節点力(剛域上の荷重は含まない)
    //   u_b0  : 節点を固定したときの内側の梁の端の変位(材端ばねの変形)
    void ReduceBeamLoad(const Eigen::Vector4d &f_b, Eigen::Vector4d &f_node, Eigen::Vector4d &u_b0) const;
};

class BeamElement : public BarElementBase
{
protected:
    static constexpr int total_dof = 12;
    double element_length();
    Eigen::MatrixXd trans_matrix();
    virtual Eigen::MatrixXd stiffness_matrix_local();

    Eigen::MatrixXd geometric_local_stiffness_matrix(const std::vector<Displacement> &disp) override;

private:
    // static constexpr ElementType type = ElementType::Beam;

public:
    double Beta;
    // Node* Nodes[2];
    // Section* Sec;
    Vector XAxis()
    {
        Eigen::Matrix3d tr0 = trans_matrix3(Nodes[0]->Location, Nodes[1]->Location, Beta);
        Eigen::Vector3d vec = tr0.row(0);
        return Vector(vec.x(), vec.y(), vec.z());
    }

    Vector YAxis()
    {
        Eigen::Matrix3d tr0 = trans_matrix3(Nodes[0]->Location, Nodes[1]->Location, Beta);
        Eigen::Vector3d vec = tr0.row(1);
        return Vector(vec.x(), vec.y(), vec.z());
    }

    Vector ZAxis()
    {
        Eigen::Matrix3d tr0 = trans_matrix3(Nodes[0]->Location, Nodes[1]->Location, Beta);
        Eigen::Vector3d vec = tr0.row(2);
        return Vector(vec.x(), vec.y(), vec.z());
    }

    BeamElement() {}
    // BeamElement() : ElementBase(ElementType::Beam) {}
    BeamElement(Node *n0, Node *n1, Section *sec, Material mat, double beta = 0);
    BeamElement(int _id, Node *n0, Node *n1, Section *sec, Material mat, double beta = 0)
        : BeamElement(n0, n1, sec, mat, beta)
    {
        id = _id;
    };

    virtual Eigen::MatrixXd StiffnessMatrix();
    bool hasRotate() { return true; }
    int NodeNum() override { return node_num; }
    // ElementType Type() { return type; }
    ElementType Type() { return ElementType::Beam; }
    int TotalDof() { return total_dof; }

    // Triplet方式での行列組立
    void GetStiffnessTriplets(std::vector<Eigen::Triplet<double>>& triplets) override;
    void GetGeometricStiffnessTriplets(const std::vector<Displacement>& disp,
                                        std::vector<Eigen::Triplet<double>>& triplets) override;

    std::vector<NodeLoadData> InertialForceToNodeLoadData(Eigen::Vector3d accel_vec) override;
    Eigen::MatrixXd NodeConsistentMass() override;

    // double length();
    BeamStress stress(Displacement d0, Displacement d1);
    // p は節点iから測った位置を全長で割った値(0~1)
    virtual Displacement DisplaceAt(Displacement d0, Displacement d1, double p);

    // 曲げの面ごとの材端条件。通常の梁は剛域なし・全固定
    virtual BeamBendingEnds BendingEnds(BendingPlane plane);
};

#endif
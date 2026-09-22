#ifndef _SUPPORT_SPRING_ELEMENT_
#define _SUPPORT_SPRING_ELEMENT_

#include "ElementBase.h"

// 支点ばね要素。節点と地面を結ぶ1節点のばねで、全体座標系の自由度ごとに
// 独立したばね定数を持つ。
//
// 地面側の端は変位0で消去済みのため、剛性行列は節点の6自由度の対角に
// ばね定数を並べた 6x6 となり、全体自由度は増えない。ばねの自由度は拘束されず
// 解く側に入るので、固定自由度の反力には現れない。反力は AddReaction() で補う。
//
// 質量・慣性力・幾何剛性は持たない。値の妥当性(負のばね定数など)は検査しない。
class SupportSpringElement : public ElementBase
{
private:
    static constexpr int total_dof = 6;
    static constexpr int node_num = 1;

    // 幾何剛性は持たない(ゼロ行列)
    Eigen::MatrixXd geometric_local_stiffness_matrix(const std::vector<Displacement> &disp) override;

public:
    Node *Nodes[node_num];
    std::array<double, 6> K{}; // 全体座標系・自由度ごとのばね定数(Dx, Dy, Dz, Rx, Ry, Rz)

    SupportSpringElement() : Nodes{nullptr} {}
    SupportSpringElement(Node *n, double kx, double ky, double kz,
                         double krx, double kry, double krz);
    SupportSpringElement(int _id, Node *n, double kx, double ky, double kz,
                         double krx, double kry, double krz)
        : SupportSpringElement(n, kx, ky, kz, krx, kry, krz)
    {
        id = _id;
    };

    int NodeNum() override { return node_num; }
    std::vector<Node *> NodesList() override { return std::vector<Node *>{Nodes[0]}; }
    ElementType Type() override { return ElementType::SupportSpring; }
    int TotalDof() override { return total_dof; }

    // ばね定数を対角に並べた剛性行列(6x6)
    Eigen::MatrixXd StiffnessMatrix() override;

    [[deprecated("Use GetStiffnessTriplets() instead. This coeffRef-based method is less efficient and will be removed in a future version.")]]
    void AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat) override;
    [[deprecated("Use GetGeometricStiffnessTriplets() instead. This coeffRef-based method is less efficient and will be removed in a future version.")]]
    void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp) override {}
    void AssembleMassMatrix(Eigen::SparseMatrix<double> &mat) override {}

    // Triplet方式での行列組立
    void GetStiffnessTriplets(std::vector<Eigen::Triplet<double>> &triplets) override;
    void GetGeometricStiffnessTriplets(const std::vector<Displacement> &disp,
                                       std::vector<Eigen::Triplet<double>> &triplets) override {}

    // 質量・慣性力は持たない
    Eigen::VectorXd NodeLumpedMass() override { return Eigen::VectorXd::Zero(node_num); }
    Eigen::MatrixXd NodeConsistentMass() override { return Eigen::MatrixXd::Zero(total_dof, total_dof); }
    std::vector<NodeLoadData> InertialForceToNodeLoadData(Eigen::Vector3d accel_vec) override { return {}; }

    // ばね定数が非零の自由度(反力を報告する自由度)
    std::array<bool, 6> ActiveDOFs() const;

    // ばねが構造へ及ぼす力 R = -(k・u + c・v) を全体反力ベクトル r_full に加算する。
    //   v_full    : 速度ベクトル(減衰力を含める場合。不要なら nullptr)
    //   damp_coef : 減衰マトリクスの剛性比例成分の係数 a (C = a・K + ...)。c = a・k となる
    void AddReaction(const Eigen::VectorXd &u_full, Eigen::VectorXd &r_full,
                     const Eigen::VectorXd *v_full = nullptr, double damp_coef = 0.0) const;
};

#endif

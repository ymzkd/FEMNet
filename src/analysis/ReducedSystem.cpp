#include "ReducedSystem.h"

#include "SparseMatrixUtils.h"

namespace
{
// 回転自由度に剛性が「付いていない」とみなすしきい値(回転自由度の最大対角成分に対する比)。
// 剛性行列の対角成分がこの値以下の回転自由度は、剛性のない死んだ自由度として
// 自動的に拘束する(トラスのみが接続する節点の回転、板要素の面内回転など)。
//   - 絶対値ではなく比で判定する(単位系に依存しないため)
//   - 比較の基準は「回転自由度だけの最大対角成分」とする。並進(N/mm)と回転
//     (N・mm/rad)では単位が異なり、全体の最大値を基準にすると柔らかい回転自由度を
//     取りこぼす/拾いすぎるため
//   - 丸め誤差程度の剛性しか付かない自由度まで拾うための余裕。厳密にゼロのみを
//     対象にしたい場合は 0.0 にする
//   - ごく小さい実剛性を持つ自由度まで拘束してしまう場合は小さくする
constexpr double kDeadRotationDiagRelTol = 1e-12;
} // namespace

ReducedSystem::ReducedSystem(FEModel &model, const Eigen::SparseMatrix<double> &stiffness)
    : slave_indices(model.RigidLinkData->SlaveDOFIndices()),
      linkTransMat(model.RigidLinkData->TransformationMatrix().sparseView(1e-10))
{
    master_dof_num = (int)linkTransMat.cols();
    Classify(model, stiffness);
}

// fixed と free は昇順で構築する(splitMatrix3x3 の残りグループ、および
// splitMatrixWithResize の自由側が昇順であることと整合させるため)。
void ReducedSystem::Classify(FEModel &model, const Eigen::SparseMatrix<double> &stiffness)
{
    const int dof_num = model.DOFNum();
    std::vector<bool> is_slave(dof_num, false);
    for (const int idx : slave_indices)
        if (idx >= 0 && idx < dof_num)
            is_slave[idx] = true;

    // 剛性が付かない回転自由度のしきい値(回転自由度の最大対角成分に対する比)
    const Eigen::VectorXd kdiag = stiffness.diagonal();
    double rot_diag_max = 0.0;
    for (Eigen::Index i = 0; i < kdiag.size(); i++)
        if ((i % NODE_DOF) >= 3 && kdiag(i) > rot_diag_max)
            rot_diag_max = kdiag(i);
    const double dead_tol = kDeadRotationDiagRelTol * rot_diag_max;

    for (int i = 0; i < dof_num; i++)
    {
        if (is_slave[i])
            continue;

        const Support &sup = model.Nodes[i / NODE_DOF].Fix;
        bool fixed = (sup.BoundaryTypes[i % NODE_DOF] == ConstraintType::Fix);
        if (!fixed && (i % NODE_DOF) >= 3)
        {
            // 回転自由度で剛性が付かないものは自動拘束する。支点ばね要素の剛性も
            // 対角に含まれるので、回転ばねで支えた自由度は拘束されない。
            // (並進自由度は本当に不安定なモデルなので各解析側で例外にする)
            if (i < kdiag.size() && kdiag(i) <= dead_tol)
                fixed = true;
        }

        if (fixed)
            fixed_indices.push_back(i);
        else
            free_indices.push_back(i);
    }
}

std::string DescribeReducedDOF(const ReducedSystem &rs, int reduced_index, FEModel &model)
{
    if (reduced_index < 0)
        return "invalid DOF";
    if (reduced_index < rs.master_dof_num)
        return "rigid link master DOF " + std::to_string(reduced_index);

    int i = reduced_index - rs.master_dof_num;
    if (i >= (int)rs.free_indices.size())
        return "reduced DOF " + std::to_string(reduced_index);

    static const char *kDofName[NODE_DOF] = {"Ux", "Uy", "Uz", "Rx", "Ry", "Rz"};
    int global = rs.free_indices[i];
    return "node " + std::to_string(global / NODE_DOF) + ", " + kDofName[global % NODE_DOF];
}

void ReducedSystem::Reduce(const Eigen::SparseMatrix<double> &full,
                           Eigen::SparseMatrix<double> &Aaa,
                           Eigen::SparseMatrix<double> *Aab,
                           Eigen::SparseMatrix<double> *Abb) const
{
    if (master_dof_num > 0)
    {
        // RigidLinkがある場合: 3x3ブロックに分割して縮小
        Eigen::SparseMatrix<double> a11, a12, a13, a22, a23, a33;
        SparseMatrixUtils::splitMatrix3x3(full, slave_indices, free_indices,
                                          a11, a12, a13, a22, a23, a33);

        Eigen::SparseMatrix<double> aaa, aab;
        aaa = (linkTransMat.transpose() * a11.selfadjointView<Eigen::Upper>() * linkTransMat)
                  .triangularView<Eigen::Upper>();
        aab = (linkTransMat.transpose() * a12);
        SparseMatrixUtils::mergeMatrixWithResize(aaa, aab, a22, Aaa);

        if (Aab)
            *Aab = SparseMatrixUtils::vstack(linkTransMat.transpose() * a13, a23);
        if (Abb)
            *Abb = a33;
    }
    else
    {
        // RigidLinkがない場合: 従来通り2x2分割
        if (Aab || Abb)
        {
            Eigen::SparseMatrix<double> ab, bb;
            SparseMatrixUtils::splitMatrixWithResize(full, fixed_indices, Aaa, ab, bb);
            if (Aab)
                *Aab = ab;
            if (Abb)
                *Abb = bb;
        }
        else
        {
            SparseMatrixUtils::splitMatrixWithResize(full, fixed_indices, Aaa);
        }
    }
}

void ReducedSystem::ReduceVector(const Eigen::VectorXd &full,
                                 Eigen::VectorXd &reduced, Eigen::VectorXd &fix) const
{
    Eigen::VectorXd f_slave(slave_indices.size());
    Eigen::VectorXd f_free(free_indices.size());
    for (size_t i = 0; i < slave_indices.size(); i++)
        f_slave(i) = full(slave_indices[i]);
    for (size_t i = 0; i < free_indices.size(); i++)
        f_free(i) = full(free_indices[i]);

    Eigen::VectorXd f_master = linkTransMat.transpose() * f_slave;
    reduced.resize(f_master.size() + f_free.size());
    reduced << f_master, f_free;

    fix.resize(fixed_indices.size());
    for (size_t i = 0; i < fixed_indices.size(); i++)
        fix(i) = full(fixed_indices[i]);
}

Eigen::VectorXd ReducedSystem::ExpandVector(const Eigen::VectorXd &reduced, int full_size) const
{
    Eigen::VectorXd d = Eigen::VectorXd::Zero(full_size);
    Eigen::VectorXd d_slave = linkTransMat * reduced.head(master_dof_num);
    Eigen::VectorXd d_free = reduced.tail(free_indices.size());

    for (size_t i = 0; i < slave_indices.size(); i++)
        d(slave_indices[i]) = d_slave(i);
    for (size_t i = 0; i < free_indices.size(); i++)
        d(free_indices[i]) = d_free(i);
    return d;
}

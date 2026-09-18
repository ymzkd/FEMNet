#include "ReducedSystem.h"

#include "SparseMatrixUtils.h"

namespace
{
// 回転自由度に剛性が「付いていない」とみなすしきい値(最大対角成分に対する比)。
// 剛性行列の対角成分がこの値以下の回転自由度は、剛性のない死んだ自由度として
// 自動的に拘束する(トラスのみが接続する節点の回転、板要素の面内回転など)。
//   - 絶対値ではなく最大対角成分との比で判定する(単位系に依存しないため)
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

    // 剛性が付かない回転自由度のしきい値(全体最大の対角成分に対する比)
    const Eigen::VectorXd kdiag = stiffness.diagonal();
    const double diag_max = (kdiag.size() > 0) ? kdiag.maxCoeff() : 0.0;
    const double dead_tol = kDeadRotationDiagRelTol * ((diag_max > 0.0) ? diag_max : 0.0);

    for (int i = 0; i < dof_num; i++)
    {
        if (is_slave[i])
            continue;

        const Support &sup = model.Nodes[i / NODE_DOF].Fix;
        const ConstraintType type = sup.BoundaryTypes[i % NODE_DOF];

        bool fixed = (type == ConstraintType::Fix);
        if (!fixed && type == ConstraintType::Free && (i % NODE_DOF) >= 3)
        {
            // 回転自由度で剛性が付かないものは自動拘束する。
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

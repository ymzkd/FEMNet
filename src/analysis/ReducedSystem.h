#ifndef _REDUCED_SYSTEM_
#define _REDUCED_SYSTEM_

#include <vector>
#include <Eigen/Sparse>

#include "Model.h"
#include "SparseMatrixUtils.h"

// 回転自由度に剛性が「付いていない」とみなすしきい値(最大対角成分に対する比)。
// 剛性行列の対角成分がこの値以下の回転自由度は、剛性のない死んだ自由度として
// 自動的に拘束する(トラスのみが接続する節点の回転、板要素の面内回転など)。
//   - 絶対値ではなく最大対角成分との比で判定する(単位系に依存しないため)
//   - 丸め誤差程度の剛性しか付かない自由度まで拾うための余裕。厳密にゼロのみを
//     対象にしたい場合は 0.0 にする
//   - ごく小さい実剛性を持つ自由度まで拘束してしまう場合は小さくする
static constexpr double kDeadRotationDiagRelTol = 1e-12;

// RigidLink縮約系のヘルパー(SWIG非公開・内部実装専用)。
// 全体自由度を slave(剛体リンク従属) / free / fixed に分割し、
// 全体対称行列(上三角格納)や荷重ベクトルを縮約空間(master+free)へ縮小する。
// 各解析Operatorで重複していた縮約処理を集約したもの。
//
// 自由度の分類(優先順位): slave > fixed > free
//   slave: 剛体リンクの従属自由度(マスタに従うため、拘束指定があっても従属が優先)
//   fixed: ユーザーが Fix を指定した自由度 + 剛性が付かない回転自由度(自動拘束)
//   free : 残り(ばね支持の自由度もここに入り、剛性行列側でばね剛性を受け持つ)
class ReducedSystem
{
public:
    std::vector<int> slave_indices;
    std::vector<int> free_indices;
    std::vector<int> fixed_indices;
    Eigen::SparseMatrix<double> linkTransMat; // 剛体リンク変換行列 T
    int master_dof_num = 0;

    // stiffness: 組立済みの全体剛性マトリクス(上三角格納)。
    // 剛性が付かない回転自由度の自動拘束判定に用いる。
    ReducedSystem(FEModel &model, const Eigen::SparseMatrix<double> &stiffness)
        : slave_indices(model.RigidLinkData->SlaveDOFIndices()),
          linkTransMat(model.RigidLinkData->TransformationMatrix().sparseView(1e-10))
    {
        master_dof_num = (int)linkTransMat.cols();
        Classify(model, stiffness);
    }

private:
    // 全体自由度を slave / fixed / free に分類する。
    // fixed と free は昇順で構築する(splitMatrix3x3 の残りグループ、および
    // splitMatrixWithResize の自由側が昇順であることと整合させるため)。
    void Classify(FEModel &model, const Eigen::SparseMatrix<double> &stiffness)
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

public:

    // 縮約空間(master+free)の自由度数
    int ReducedSize() const { return master_dof_num + (int)free_indices.size(); }

    // 全体対称行列(上三角格納)を縮約する。
    //   Aaa: 縮約空間(master+free)ブロック
    //   Aab: 縮約空間×固定ブロック(反力計算用, 不要ならnullptr)
    //   Abb: 固定×固定ブロック(不要ならnullptr)
    void Reduce(const Eigen::SparseMatrix<double> &full,
                Eigen::SparseMatrix<double> &Aaa,
                Eigen::SparseMatrix<double> *Aab = nullptr,
                Eigen::SparseMatrix<double> *Abb = nullptr) const
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

    // 全体荷重ベクトルを縮約する。
    //   reduced: [T^T・f_slave ; f_free] (右辺用, サイズ ReducedSize())
    //   fix    : 固定DOF成分(反力計算用)
    void ReduceVector(const Eigen::VectorXd &full,
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

    // 縮約空間の解ベクトルを全体空間へ展開する(slave = T・master, 固定DOFは0)
    Eigen::VectorXd ExpandVector(const Eigen::VectorXd &reduced, int full_size) const
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
};

#endif

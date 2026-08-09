#ifndef _REDUCED_SYSTEM_
#define _REDUCED_SYSTEM_

#include <vector>
#include <Eigen/Sparse>

#include "Model.h"
#include "SparseMatrixUtils.h"

// RigidLink縮約系のヘルパー(SWIG非公開・内部実装専用)。
// 全体自由度を slave(剛体リンク従属) / free / fixed に分割し、
// 全体対称行列(上三角格納)や荷重ベクトルを縮約空間(master+free)へ縮小する。
// 各解析Operatorで重複していた縮約処理を集約したもの。
class ReducedSystem
{
public:
    std::vector<int> slave_indices;
    std::vector<int> free_indices;
    std::vector<int> fixed_indices;
    Eigen::SparseMatrix<double> linkTransMat; // 剛体リンク変換行列 T
    int master_dof_num = 0;

    explicit ReducedSystem(FEModel &model)
        : slave_indices(model.RigidLinkData->SlaveDOFIndices()),
          free_indices(model.FreeIndices(true)),
          fixed_indices(model.FixIndices()),
          linkTransMat(model.RigidLinkData->TransformationMatrix().sparseView(1e-10))
    {
        master_dof_num = (int)linkTransMat.cols();
    }

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

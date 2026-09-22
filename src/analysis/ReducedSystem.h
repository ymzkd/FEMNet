#ifndef _REDUCED_SYSTEM_
#define _REDUCED_SYSTEM_

#include <string>
#include <vector>
#include <Eigen/Sparse>

#include "Model.h"

// RigidLink縮約系のヘルパー(SWIG非公開・内部実装専用)。
// 全体自由度を slave(剛体リンク従属) / free / fixed に分割し、
// 全体対称行列(上三角格納)や荷重ベクトルを縮約空間(master+free)へ縮小する。
// 各解析Operatorで重複していた縮約処理を集約したもの。
//
// 自由度の分類(優先順位): slave > fixed > free
//   slave: 剛体リンクの従属自由度(マスタに従うため、拘束指定があっても従属が優先)
//   fixed: ユーザーが Fix を指定した自由度 + 剛性が付かない回転自由度(自動拘束)
//   free : 残り(支点ばね要素が付く自由度もここに入り、剛性行列側でばね剛性を受け持つ)
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
    ReducedSystem(FEModel &model, const Eigen::SparseMatrix<double> &stiffness);

    // 縮約空間(master+free)の自由度数
    int ReducedSize() const { return master_dof_num + (int)free_indices.size(); }

    // 全体対称行列(上三角格納)を縮約する。
    //   Aaa: 縮約空間(master+free)ブロック
    //   Aab: 縮約空間×固定ブロック(反力計算用, 不要ならnullptr)
    //   Abb: 固定×固定ブロック(不要ならnullptr)
    void Reduce(const Eigen::SparseMatrix<double> &full,
                Eigen::SparseMatrix<double> &Aaa,
                Eigen::SparseMatrix<double> *Aab = nullptr,
                Eigen::SparseMatrix<double> *Abb = nullptr) const;

    // 全体荷重ベクトルを縮約する。
    //   reduced: [T^T・f_slave ; f_free] (右辺用, サイズ ReducedSize())
    //   fix    : 固定DOF成分(反力計算用)
    void ReduceVector(const Eigen::VectorXd &full,
                      Eigen::VectorXd &reduced, Eigen::VectorXd &fix) const;

    // 縮約空間の解ベクトルを全体空間へ展開する(slave = T・master, 固定DOFは0)
    Eigen::VectorXd ExpandVector(const Eigen::VectorXd &reduced, int full_size) const;

private:
    // 全体自由度を slave / fixed / free に分類する。
    void Classify(FEModel &model, const Eigen::SparseMatrix<double> &stiffness);
};

// 縮約空間のインデックスが、どの節点のどの自由度に対応するかを表す文字列を返す
// (剛性ゼロ等のエラーメッセージ用。例: "node 12, Rz")。
std::string DescribeReducedDOF(const ReducedSystem &rs, int reduced_index, FEModel &model);

#endif

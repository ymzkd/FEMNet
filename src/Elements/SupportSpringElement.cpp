#include "SupportSpringElement.h"

SupportSpringElement::SupportSpringElement(Node *n, double kx, double ky, double kz,
                                           double krx, double kry, double krz)
{
    Nodes[0] = n;
    K = {kx, ky, kz, krx, kry, krz};
}

Eigen::MatrixXd SupportSpringElement::geometric_local_stiffness_matrix(
    const std::vector<Displacement> &disp)
{
    return Eigen::MatrixXd::Zero(total_dof, total_dof);
}

// ばね定数を対角に並べた剛性行列(6x6)
Eigen::MatrixXd SupportSpringElement::StiffnessMatrix()
{
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(total_dof, total_dof);
    for (int i = 0; i < total_dof; i++)
        mat(i, i) = K[i];
    return mat;
}

// 非推奨: GetStiffnessTriplets()を使用してください（coeffRef方式は非効率）
void SupportSpringElement::AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat)
{
    Eigen::MatrixXd k = StiffnessMatrix();
    for (int i = 0; i < total_dof; i++)
    {
        int idx = Nodes[0]->id * 6 + i;
        mat.coeffRef(idx, idx) += k(i, i);
    }
}

void SupportSpringElement::GetStiffnessTriplets(std::vector<Eigen::Triplet<double>> &triplets)
{
    Eigen::MatrixXd k = StiffnessMatrix();

    // 1節点6自由度をそのまま全体自由度へ対応付ける(上三角のみ追加)
    for (int i = 0; i < total_dof; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double value = k(i, j);
            if (value != 0.0)
            {
                int row = Nodes[0]->id * 6 + j;
                int col = Nodes[0]->id * 6 + i;
                triplets.emplace_back(row, col, value);
            }
        }
    }
}

std::array<bool, 6> SupportSpringElement::ActiveDOFs() const
{
    std::array<bool, 6> active{};
    for (int i = 0; i < total_dof; i++)
        active[i] = (K[i] != 0.0);
    return active;
}

void SupportSpringElement::AddReaction(const Eigen::VectorXd &u_full, Eigen::VectorXd &r_full,
                                       const Eigen::VectorXd *v_full, double damp_coef) const
{
    // つり合い (M・a + C・v + K・u)_i = f_i より、ばねが構造へ及ぼす力は
    //   R = -(k・u + c・v),  c = damp_coef・k (剛性比例成分による付随減衰)
    // 静的解析など減衰を考慮しない場合は v_full = nullptr で弾性分のみとなる。
    const bool with_damping = (v_full != nullptr && damp_coef != 0.0);
    for (int i = 0; i < total_dof; i++)
    {
        if (K[i] == 0.0)
            continue;
        int idx = Nodes[0]->id * 6 + i;
        if (idx >= u_full.size() || idx >= r_full.size())
            continue;
        double r = -K[i] * u_full(idx);
        if (with_damping && idx < v_full->size())
            r -= damp_coef * K[i] * (*v_full)(idx);
        r_full(idx) += r;
    }
}

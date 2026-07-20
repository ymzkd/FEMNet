#include "TensionTrussElement.h"

Eigen::MatrixXd TensionTrussElement::tangentstiffness_matrix_local()
{
    return TrussElement::stiffness_matrix_local() * (IsActive ? 1.0 : ReductionFactor);
}

Eigen::MatrixXd TensionTrussElement::TangentStiffnessMatrix()
{
    Eigen::Matrix2d local_stiffMat = tangentstiffness_matrix_local();
    Eigen::MatrixXd transMat = trans_matrix();
    return transMat.transpose() * local_stiffMat * transMat;
}

void TensionTrussElement::GetTangentStiffnessTriplets(std::vector<Eigen::Triplet<double>> &triplets)
{
    Eigen::MatrixXd K = TangentStiffnessMatrix();
    int total_dof = TotalDof(); // 6 for TrussElement

    // Map to global DOFs - TrussElement uses only 3 DOF/node (X,Y,Z)
    int indices[6];
    for (size_t i = 0; i < 3; i++)
        indices[i] = Nodes[0]->id * 6 + i;
    for (size_t i = 0; i < 3; i++)
        indices[i + 3] = Nodes[1]->id * 6 + i;

    // Add upper triangle to triplets
    for (int i = 0; i < total_dof; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double value = K(i, j);
            if (std::abs(value) > 1e-20)
            {
                int row = indices[i];
                int col = indices[j];
                if (row > col)
                    std::swap(row, col);
                triplets.emplace_back(row, col, value);
            }
        }
    }
}

BeamStress TensionTrussElement::tangent_stress(Displacement d0, Displacement d1)
{
    Eigen::VectorXd disp(6);
    disp << d0.Dx(), d0.Dy(), d0.Dz(), d1.Dx(), d1.Dy(), d1.Dz();

    Eigen::Vector2d force = tangentstiffness_matrix_local() * trans_matrix() * disp;

    BeamStressData str0(-force(0), 0, 0, 0, 0, 0);
    BeamStressData str1(force(1), 0, 0, 0, 0, 0);
    return BeamStress(str0, str1);
}

bool TensionTrussElement::update(const std::vector<Displacement> &disp)
{
    bool any_change = false;
    BeamStress stress = this->tangent_stress(disp[Nodes[0]->id], disp[Nodes[1]->id]);
    if (stress.S0.Nx < this->CutoffTension && this->IsActive) // 圧縮状態に移行
    {
        this->IsActive = false; // 要素を無効化
        any_change = true;
    }
    else if (stress.S0.Nx > this->CutoffTension && !this->IsActive) // 引張状態に移行
    {
        this->IsActive = true; // 要素を有効化
        any_change = true;
    }
    return any_change;
}

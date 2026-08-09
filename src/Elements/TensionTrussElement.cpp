#include "TensionTrussElement.h"

Eigen::MatrixXd TensionTrussElement::TangentStiffnessMatrix(bool active)
{
    Eigen::Matrix2d local_stiffMat =
        TrussElement::stiffness_matrix_local() * (active ? 1.0 : ReductionFactor);
    Eigen::MatrixXd transMat = trans_matrix();
    return transMat.transpose() * local_stiffMat * transMat;
}

void TensionTrussElement::GetTangentStiffnessTriplets(
    std::vector<Eigen::Triplet<double>> &triplets, bool active)
{
    Eigen::MatrixXd K = TangentStiffnessMatrix(active);
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

BeamStress TensionTrussElement::tangent_stress(Displacement d0, Displacement d1, bool active)
{
    Eigen::VectorXd disp(6);
    disp << d0.Dx(), d0.Dy(), d0.Dz(), d1.Dx(), d1.Dy(), d1.Dz();

    Eigen::MatrixXd k_local =
        TrussElement::stiffness_matrix_local() * (active ? 1.0 : ReductionFactor);
    Eigen::Vector2d force = k_local * trans_matrix() * disp;

    BeamStressData str0(-force(0), 0, 0, 0, 0, 0);
    BeamStressData str1(force(1), 0, 0, 0, 0, 0);
    return BeamStress(str0, str1);
}

bool TensionTrussElement::NextState(const std::vector<Displacement> &disp, bool current)
{
    BeamStress stress = tangent_stress(disp[Nodes[0]->id], disp[Nodes[1]->id], current);
    if (stress.S0.Nx < CutoffTension && current)
        return false; // 圧縮状態に移行 → 無効化
    if (stress.S0.Nx > CutoffTension && !current)
        return true; // 引張状態に移行 → 有効化
    return current;
}

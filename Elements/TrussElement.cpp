#include "TrussElement.h"

Eigen::Matrix<double, 2, 6> TrussElement::trans_matrix()
{
    Vector vij = Nodes[1]->Location - Nodes[0]->Location;
    double l = vij.norm();
    // double EAl = Mat->Young * Sec->A / l;
    double alpha = vij * Vector::XAxis() / l;
    double beta = vij * Vector::YAxis() / l;
    double gamma = vij * Vector::ZAxis() / l;

    Eigen::MatrixXd transMat(2, 6);
    transMat << alpha, beta, gamma, 0, 0, 0,
        0, 0, 0, alpha, beta, gamma;
    return transMat;
}

Eigen::MatrixXd TrussElement::stiffness_matrix_local()
{
    double l = (Nodes[1]->Location - Nodes[0]->Location).norm();
    double EAl = Mat.Young * Sec->A / l;
    // double alpha = vij * Vector::XAxis() / l;
    // double beta = vij * Vector::YAxis() / l;
    // double gamma = vij * Vector::ZAxis() / l;

    Eigen::Matrix2d local_stiffMat;
    local_stiffMat << EAl, -EAl,
        -EAl, EAl;
    return local_stiffMat;
}

Eigen::MatrixXd TrussElement::geometric_local_stiffness_matrix(
    const std::vector<Displacement> &disp)
{
    Eigen::MatrixXd Kg_mat = Eigen::MatrixXd::Zero(total_dof, total_dof);
    double length = (Nodes[1]->Location - Nodes[0]->Location).norm();
    double Nx = stress(disp[0], disp[1]).S0.Nx;

    Kg_mat << 1, 0, 0, -1, 0, 0,
        0, 1, 0, 0, -1, 0,
        0, 0, 1, 0, 0, -1,
        -1, 0, 0, 1, 0, 0,
        0, -1, 0, 0, 1, 0,
        0, 0, -1, 0, 0, 1;

    Kg_mat = Kg_mat * Nx / length;
    // Eigen::MatrixXd transMat = trans_matrix();
    // return transMat.transpose() * Kg_mat * transMat;
    return Kg_mat;
}

TrussElement::TrussElement(Node *n0, Node *n1, Section *sec, Material mat)
{
    Nodes[0] = n0;
    Nodes[1] = n1;

    Sec = sec;
    Mat = mat;
}

// トラス要素の剛性行列(6x6)
Eigen::MatrixXd TrussElement::StiffnessMatrix()
{
    // Eigen::MatrixXd matrix(total_dof, total_dof);
    // Vector vij = ;
    // double l = (Nodes[1]->Location - Nodes[0]->Location).norm();
    // double EAl = Mat->Young * Sec->A / l;
    // double alpha = vij * Vector::XAxis() / l;
    // double beta = vij * Vector::YAxis() / l;
    // double gamma = vij * Vector::ZAxis() / l;

    Eigen::Matrix2d local_stiffMat = stiffness_matrix_local();
    /*local_stiffMat << EAl, -EAl,
                    - EAl, EAl;*/

    Eigen::MatrixXd transMat = trans_matrix();

    return transMat.transpose() * local_stiffMat * transMat;
}

void TrussElement::AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K)
{
    int indices[6];
    for (size_t i = 0; i < 3; i++)
        indices[i] = Nodes[0]->id * 6 + i;
    for (size_t i = 0; i < 3; i++)
        indices[i + 3] = Nodes[1]->id * 6 + i;

    // Eigen::MatrixXd smat = StiffnessMatrix();
    for (size_t i = 0; i < 6; i++)
    {

        for (size_t j = 0; j < i + 1; j++)
        {
            int ci, rj;
            if (indices[j] <= indices[i])
            {
                ci = indices[i];
                rj = indices[j];
            }
            else
            {
                ci = indices[j];
                rj = indices[i];
            }
            mat.coeffRef(rj, ci) += K(i, j);
        }
    }
}

// 非推奨: GetStiffnessTriplets()を使用してください（coeffRef方式は非効率）
void TrussElement::AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat)
{
    AssembleMatrix(mat, StiffnessMatrix());
    // int indices[6];
    // for (size_t i = 0; i < 3; i++)
    // 	indices[i] = Nodes[0]->id * 6 + i;
    // for (size_t i = 0; i < 3; i++)
    // 	indices[i + 3] = Nodes[1]->id * 6 + i;

    // Eigen::MatrixXd smat = StiffnessMatrix();
    // for (size_t i = 0; i < 6; i++)
    // {

    // 	for (size_t j = 0; j < i + 1; j++)
    // 	{
    // 		int ci, rj;
    // 		if (indices[j] <= indices[i]) {
    // 			ci = indices[i];
    // 			rj = indices[j];
    // 		}
    // 		else {
    // 			ci = indices[j];
    // 			rj = indices[i];
    // 		}
    // 		mat.coeffRef(rj, ci) += smat(i, j);
    // 	}
    // }
}

// 非推奨: GetGeometricStiffnessTriplets()を使用してください（coeffRef方式は非効率）
void TrussElement::AssembleGeometricStiffMatrix(
    Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp)
{
    AssembleMatrix(mat, geometric_local_stiffness_matrix(disp));
}

void TrussElement::GetStiffnessTriplets(std::vector<Eigen::Triplet<double>>& triplets)
{
    Eigen::MatrixXd K = StiffnessMatrix();
    int total_dof = TotalDof();  // 6 for TrussElement

    // Map to global DOFs - TrussElement uses only 3 DOF/node (X,Y,Z)
    int indices[6];
    for (size_t i = 0; i < 3; i++)
        indices[i] = Nodes[0]->id * 6 + i;
    for (size_t i = 0; i < 3; i++)
        indices[i + 3] = Nodes[1]->id * 6 + i;

    // Add lower triangle to triplets
    for (int i = 0; i < total_dof; i++) {
        for (int j = 0; j <= i; j++) {
            double value = K(i, j);
            if (std::abs(value) > 1e-20) {
                int row = indices[i];
                int col = indices[j];
                if (col > row) std::swap(row, col);
                triplets.emplace_back(row, col, value);
            }
        }
    }
}

void TrussElement::GetGeometricStiffnessTriplets(
    const std::vector<Displacement>& disp,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    Eigen::MatrixXd Kg = geometric_local_stiffness_matrix(disp);
    int total_dof = TotalDof();  // 6 for TrussElement

    // Map to global DOFs - TrussElement uses only 3 DOF/node (X,Y,Z)
    int indices[6];
    for (size_t i = 0; i < 3; i++)
        indices[i] = Nodes[0]->id * 6 + i;
    for (size_t i = 0; i < 3; i++)
        indices[i + 3] = Nodes[1]->id * 6 + i;

    // Add lower triangle to triplets
    for (int i = 0; i < total_dof; i++) {
        for (int j = 0; j <= i; j++) {
            double value = Kg(i, j);
            if (std::abs(value) > 1e-20) {
                int row = indices[i];
                int col = indices[j];
                if (col > row) std::swap(row, col);
                triplets.emplace_back(row, col, value);
            }
        }
    }
}

Eigen::MatrixXd TrussElement::NodeConsistentMass()
{
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(total_dof, total_dof);
    m << 2, 0, 0, 1, 0, 0,
        0, 2, 0, 0, 1, 0,
        0, 0, 2, 0, 0, 1,
        1, 0, 0, 2, 0, 0,
        0, 1, 0, 0, 2, 0,
        0, 0, 1, 0, 0, 2;

    m *= (Sec->A * Mat.dense / 6.0);
    return m;
}

std::vector<NodeLoadData> TrussElement::InertialForceToNodeLoadData(Eigen::Vector3d accel_vec)
{
    Eigen::Matrix<double, total_dof, total_dof> m = NodeConsistentMass();
    Eigen::Vector<double, total_dof> f;
    f << accel_vec.x(), accel_vec.y(), accel_vec.z(),
        accel_vec.x(), accel_vec.y(), accel_vec.z();

    f = m * f;
    std::vector<NodeLoadData> loads;
    loads.push_back(NodeLoadData(Nodes[0]->id, f(0), f(1), f(2)));
    loads.push_back(NodeLoadData(Nodes[1]->id, f(3), f(4), f(5)));
    return loads;
}

BeamStress TrussElement::stress(Displacement d0, Displacement d1)
{
    Eigen::VectorXd disp(6);
    disp << d0.Dx(), d0.Dy(), d0.Dz(), d1.Dx(), d1.Dy(), d1.Dz();

    Eigen::Vector2d force = stiffness_matrix_local() * trans_matrix() * disp;

    BeamStressData str0(-force(0), 0, 0, 0, 0, 0);
    BeamStressData str1(force(1), 0, 0, 0, 0, 0);
    return BeamStress(str0, str1);
}

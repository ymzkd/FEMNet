#include "BarElement.h"

std::ostream &operator<<(std::ostream &os, const BeamStressData &bsd)
{
    /*std::string force = std::format("Nx: {}, Qy: {}, Qz: {}, Mx: {}, My: {}, Mz: {}",
        bsd.Nx, bsd.Qy, bsd.Qz, bsd.Mx, bsd.My, bsd.Mz);
    os << force << std::endl;*/

    os << "Nx: " << bsd.Nx << ", Qy: " << bsd.Qy << ", Qz: " << bsd.Qz
       << ", Mx: " << bsd.Mx << ", My: " << bsd.My << ", Mz: " << bsd.Mz << std::endl;
    return os;
}

std::ostream &operator<<(std::ostream &os, const BeamStress &bsd)
{
    os << "Beam Stress" << std::endl;
    os << "Start  " << bsd.S0;
    os << "End    " << bsd.S1;
    return os;
}

double BarElementBase::length()
{
    return Nodes[0]->Location.distance_to(Nodes[1]->Location);
}

void BarElementBase::AssembleMassMatrix(Eigen::SparseMatrix<double> &mat)
{
    int indices[6];
    Eigen::VectorXd mass = NodeLumpedMass();
    for (size_t i = 0; i < 3; i++)
    {
        int idx = Nodes[0]->id * 6 + i;
        mat.coeffRef(idx, idx) += mass[0];
    }
    for (size_t i = 0; i < 3; i++)
    {
        int idx = Nodes[1]->id * 6 + i;
        mat.coeffRef(idx, idx) += mass[1];
    }
}

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

void TrussElement::AssembleGeometricStiffMatrix(
    Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp)
{
    AssembleMatrix(mat, geometric_local_stiffness_matrix(disp));
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

BeamElement::BeamElement(Node *n0, Node *n1, Section *sec, Material mat, double beta) : BeamElement()
{
    Nodes[0] = n0;
    Nodes[1] = n1;

    Sec = sec;
    Mat = mat;

    Beta = beta;
}

Eigen::MatrixXd BeamElement::StiffnessMatrix()
{
    Eigen::MatrixXd tr = trans_matrix();
    Eigen::MatrixXd k = stiffness_matrix_local();
    return tr.transpose() * k * tr;
}

void BeamElement::AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K)
{
    int indices[12];
    for (size_t i = 0; i < 6; i++)
        indices[i] = Nodes[0]->id * 6 + i;
    for (size_t i = 0; i < 6; i++)
        indices[i + 6] = Nodes[1]->id * 6 + i;

    // Eigen::MatrixXd smat = StiffnessMatrix();
    for (size_t i = 0; i < 12; i++)
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

void BeamElement::AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat)
{
    AssembleMatrix(mat, StiffnessMatrix());
    // int indices[12];
    // for (size_t i = 0; i < 6; i++)
    // 	indices[i] = Nodes[0]->id * 6 + i;
    // for (size_t i = 0; i < 6; i++)
    // 	indices[i + 6] = Nodes[1]->id * 6 + i;

    // Eigen::MatrixXd smat = StiffnessMatrix();
    // for (size_t i = 0; i < 12; i++)
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
    // 		mat.coeffRef(rj,ci) += smat(i, j);
    // 	}
    // }
}

void BeamElement::AssembleGeometricStiffMatrix(
    Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp)
{
    AssembleMatrix(mat, geometric_local_stiffness_matrix(disp));
}

std::vector<NodeLoadData> BeamElement::InertialForceToNodeLoadData(Eigen::Vector3d accel_vec)
{
    Eigen::Matrix<double, total_dof, total_dof> m;
    m = NodeConsistentMass();
    Eigen::Vector<double, total_dof> f;
    f.setZero();

    f << accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0,
        accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0;

    f = m * f;

    std::cout << "Mass Matrix: " << m.diagonal() << std::endl;
    // std::cout << "BodyforceToNodeLoadData: " << f.transpose() << std::endl;

    std::vector<NodeLoadData> loads;
    loads.push_back(NodeLoadData(Nodes[0]->id, f(0), f(1), f(2), f(3), f(4), f(5)));
    loads.push_back(NodeLoadData(Nodes[1]->id, f(6), f(7), f(8), f(9), f(10), f(11)));
    return loads;
}

Eigen::MatrixXd BeamElement::NodeConsistentMass()
{
    Eigen::MatrixXd tr = trans_matrix();
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(total_dof, total_dof);
    double l = element_length();
    double l2 = l * l;
    double l3 = l2 * l;
    // m << l / 3, 0, 0, 0, 0, 0, l / 6, 0, 0, 0, 0, 0,
    //	0, 13.0 / 35.0 * l, 0, 0, 0, 11 / 210 * l2, 0, 9 / 70 * l, 0, 0, 0, -13 / 420 * l2,
    //	0, 0, 13 / 35 * l, 0, -11 / 210 * l2, 0, 0, 0, 9 / 70 * l, 0, 13 / 420 * l2, 0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	0, 0, -11 / 210 * l2, 0, 1 / 105 * l3, 0, 0, 0, -13 / 420 * l2, 0, -1 / 140 * l3, 0,
    //	0, 11 / 210 * l2, 0, 0, 0, 1 / 105 * l3, 0, 13 / 420 * l2, 0, 0, 0, -1 / 140 * l3,
    //	l / 6, 0, 0, 0, 0, 0, l / 3, 0, 0, 0, 0, 0,
    //	0, 9 / 70 * l, 0, 0, 0, 13 / 420 * l2, 0, 13 / 35 * l, 0, 0, 0, -11 / 210 * l2,
    //	0, 0, 9 / 70 * l, 0, -13 / 420 * l2, 0, 0, 0, 13 / 35 * l, 0, 11 / 210 * l2, 0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	0, 0, 13 / 420 * l2, 0, -1 / 140 * l3, 0, 0, 0, 11 / 210 * l2, 0, 1 / 105 * l3, 0,
    //	0, -13 / 420 * l2, 0, 0, 0, -1 / 140 * l3, 0, -11 / 210 * l2, 0, 0, 0, 1 / 105 * l3;

    m << l / 3.0, 0, 0, 0, 0, 0, l / 6.0, 0, 0, 0, 0, 0,
        0, 13.0 / 35.0 * l, 0, 0, 0, 11.0 / 210.0 * l2, 0, 9.0 / 70.0 * l, 0, 0, 0, -13.0 / 420.0 * l2,
        0, 0, 13.0 / 35.0 * l, 0, -11.0 / 210.0 * l2, 0, 0, 0, 9.0 / 70.0 * l, 0, 13.0 / 420.0 * l2, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, -11.0 / 210.0 * l2, 0, 1.0 / 105.0 * l3, 0, 0, 0, -13.0 / 420.0 * l2, 0, -1.0 / 140.0 * l3, 0,
        0, 11.0 / 210.0 * l2, 0, 0, 0, 1.0 / 105.0 * l3, 0, 13.0 / 420.0 * l2, 0, 0, 0, -1.0 / 140.0 * l3,
        l / 6.0, 0, 0, 0, 0, 0, l / 3.0, 0, 0, 0, 0, 0,
        0, 9.0 / 70.0 * l, 0, 0, 0, 13.0 / 420.0 * l2, 0, 13.0 / 35.0 * l, 0, 0, 0, -11.0 / 210.0 * l2,
        0, 0, 9.0 / 70.0 * l, 0, -13.0 / 420.0 * l2, 0, 0, 0, 13.0 / 35.0 * l, 0, 11.0 / 210.0 * l2, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 13.0 / 420.0 * l2, 0, -1.0 / 140.0 * l3, 0, 0, 0, 11.0 / 210.0 * l2, 0, 1.0 / 105.0 * l3, 0,
        0, -13.0 / 420.0 * l2, 0, 0, 0, -1.0 / 140.0 * l3, 0, -11.0 / 210.0 * l2, 0, 0, 0, 1.0 / 105.0 * l3;

    return tr.transpose() * (Sec->A * Mat.dense) * m * tr;
}

BeamStress BeamElement::stress(Displacement d0, Displacement d1)
{
    Eigen::VectorXd wvec(12);
    wvec.segment(0, 6) = Eigen::Map<Eigen::VectorXd>(d0.displace, 6);
    wvec.segment(6, 6) = Eigen::Map<Eigen::VectorXd>(d1.displace, 6);

    wvec = trans_matrix() * wvec;
    wvec = stiffness_matrix_local() * wvec;
    // wvec = StiffnessMatrix() * wvec;

    // N, Qy, Qz, Mx, My, Mz
    BeamStressData s0(-wvec(0), -wvec(1), -wvec(2), -wvec(3), -wvec(4), wvec(5));
    BeamStressData s1(wvec(6), wvec(7), wvec(8), wvec(9), wvec(10), -wvec(11));
    return BeamStress(s0, s1);
}

Displacement BeamElement::DisplaceAt(Displacement d0, Displacement d1, double p)
{
    Eigen::VectorXd wvec(12);
    wvec.segment(0, 6) = Eigen::Map<Eigen::VectorXd>(d0.displace, 6);
    wvec.segment(6, 6) = Eigen::Map<Eigen::VectorXd>(d1.displace, 6);

    Eigen::MatrixXd T = trans_matrix();
    wvec = T * wvec;
    Eigen::Matrix<double, 6, total_dof> m;
    double l = length();
    double x = p * l;
    double p2 = p * p;
    double p3 = p * p * p;

    m << 1 - p, 0, 0, 0, 0, 0, p, 0, 0, 0, 0, 0,                                                                         // dx
        0, 1 - 3 * p2 + 2 * p3, 0, 0, 0, (1 - 2 * p + p2) * x, 0, 3 * p2 - 2 * p3, 0, 0, 0, (-p + p2) * x,               // dy v
        0, 0, 1 - 3 * p2 + 2 * p3, 0, (-1 + 2 * p - p2) * x, 0, 0, 0, 3 * p2 - 2 * p3, 0, (p - p2) * x, 0,               // dz w
        0, 0, 0, 1 - p, 0, 0, 0, 0, 0, p, 0, 0,                                                                          // rx
        0, 0, -6.0 * (p2 - p) / l, 0, 1.0 - 4.0 * p + 3.0 * p2, 0, 0, 0, -6.0 * (p - p2) / l, 0, -2.0 * p + 3.0 * p2, 0, // ry
        0, 6.0 * (p2 - p) / l, 0, 0, 0, 1.0 - 4.0 * p + 3.0 * p2, 0, 6.0 * (p - p2) / l, 0, 0, 0, 3.0 * p2 - 2.0 * p;    // rz

    Eigen::VectorXd d = T.block(0, 0, 6, 6).transpose() * m * wvec;

    return Displacement(d(0), d(1), d(2), d(3), d(4), d(5));
}

double BeamElement::element_length()
{
    Point p0 = Nodes[0]->Location;
    Point p1 = Nodes[1]->Location;
    double l2 = (p1.x - p0.x) * (p1.x - p0.x) + (p1.y - p0.y) * (p1.y - p0.y) + (p1.z - p0.z) * (p1.z - p0.z);
    return sqrt(l2);
}

Eigen::MatrixXd BeamElement::stiffness_matrix_local()
{
    double E = Mat.Young;
    double G = Mat.G();
    double l = element_length();
    double A = Sec->A;
    double Iy = Sec->Iy;
    double Iz = Sec->Iz;
    double K = Sec->K;

    double EAl = E * A / l;
    double GKl = G * K / l;

    double EIz12 = 12 * E * Iz / (l * l * l);
    double EIz6 = 6 * E * Iz / (l * l);
    double EIz2 = 2 * E * Iz / l;
    double EIz4 = EIz2 * 2;

    double EIy12 = 12 * E * Iy / (l * l * l);
    double EIy6 = 6 * E * Iy / (l * l);
    double EIy2 = 2 * E * Iy / l;
    double EIy4 = EIy2 * 2;

    Eigen::Matrix<double, total_dof, total_dof> m;
    m << EAl, 0, 0, 0, 0, 0, -EAl, 0, 0, 0, 0, 0,
        0, EIz12, 0, 0, 0, EIz6, 0, -EIz12, 0, 0, 0, EIz6,
        0, 0, EIy12, 0, -EIy6, 0, 0, 0, -EIy12, 0, -EIy6, 0,
        0, 0, 0, GKl, 0, 0, 0, 0, 0, -GKl, 0, 0,
        0, 0, -EIy6, 0, EIy4, 0, 0, 0, EIy6, 0, EIy2, 0,
        0, EIz6, 0, 0, 0, EIz4, 0, -EIz6, 0, 0, 0, EIz2,
        -EAl, 0, 0, 0, 0, 0, EAl, 0, 0, 0, 0, 0,
        0, -EIz12, 0, 0, 0, -EIz6, 0, EIz12, 0, 0, 0, -EIz6,
        0, 0, -EIy12, 0, EIy6, 0, 0, 0, EIy12, 0, EIy6, 0,
        0, 0, 0, -GKl, 0, 0, 0, 0, 0, GKl, 0, 0,
        0, 0, -EIy6, 0, EIy2, 0, 0, 0, EIy6, 0, EIy4, 0,
        0, EIz6, 0, 0, 0, EIz2, 0, -EIz6, 0, 0, 0, EIz4;

    return m;
}

Eigen::MatrixXd BeamElement::geometric_local_stiffness_matrix(const std::vector<Displacement> &disp)
{
    Eigen::MatrixXd Kg_mat = Eigen::MatrixXd::Zero(total_dof, total_dof);
    double l = element_length();
    BeamStress s = stress(disp[0], disp[1]);
    double Nx = s.S0.Nx;
    double Qy = s.S0.Qy;
    double Qz = s.S0.Qz;

    int indices_b[8]{1, 2, 4, 5, 7, 8, 10, 11};
    int indices_x[2]{0, 6};

    // 曲げによる幾何剛性
    Eigen::MatrixXd Kg_b = Eigen::MatrixXd::Zero(8, 8);
    Kg_b << 6.0 / (5.0 * l), 0, 0, 1.0 / 10.0, -6.0 / (5.0 * l), 0, 0, 1.0 / 10.0,
        0, 6.0 / (5.0 * l), -1.0 / 10.0, 0, 0, -6.0 / (5.0 * l), -1.0 / 10.0, 0,
        0, -1.0 / 10.0, 2.0 * l / 15.0, 0, 0, 1.0 / 10.0, -l / 30.0, 0,
        1.0 / 10.0, 0, 0, 2.0 * l / 15.0, -1.0 / 10.0, 0, 0, -l / 30.0,
        -6.0 / (5.0 * l), 0, 0, -1.0 / 10.0, 6.0 / (5.0 * l), 0, 0, -1.0 / 10.0,
        0, -6.0 / (5.0 * l), 1.0 / 10.0, 0, 0, 6.0 / (5.0 * l), 1.0 / 10.0, 0,
        0, -1.0 / 10.0, -l / 30.0, 0, 0, 1.0 / 10.0, 2.0 * l / 15.0, 0,
        1.0 / 10.0, 0, 0, -l / 30.0, -1.0 / 10.0, 0, 0, 2.0 * l / 15.0;
    Kg_b *= Nx;

    // 軸方向幾何剛性
    Eigen::MatrixXd Kg_x = Eigen::MatrixXd::Zero(2, 2);
    Kg_x << 1, -1,
        -1, 1;
    Kg_x *= Nx / l;

    for (size_t i = 0; i < 8; i++)
    {
        int ir = indices_b[i];
        for (size_t j = 0; j < 8; j++)
        {
            int ic = indices_b[j];
            Kg_mat(ir, ic) += Kg_b(i, j);
        }
    }

    for (size_t i = 0; i < 2; i++)
    {
        int ir = indices_x[i];
        for (size_t j = 0; j < 2; j++)
        {
            int ic = indices_x[j];
            Kg_mat(ir, ic) += Kg_x(i, j);
        }
    }

    // Kg_mat << 1.0 / l, 0, 0, 0, 0, 0, -1.0 / l, 0, 0, 0, 0, 0,
    // 	0, 6.0 / 5.0 / l, 0, 0, 0, 1.0 / 10.0, 0, -6.0 / 15.0 / l, 0, 0, 0, 1.0 / 10.0,
    // 	0, 0, 6.0 / 5.0 / l, 0, -1.0 / 10.0, 0, 0, 0, -6.0 / 15.0 / l, 0, -1.0 / 10.0, 0,
    // 	0, 0, 0, 1.0 / l, 0, 0, 0, 0, 0, -1.0 / l, 0, 0,
    // 	0, 0, -1.0 / 10.0, 0, 2.0 * l / 15.0, 0, 0, 0, 1.0 / 10.0, 0, -l / 30.0, 0,
    // 	0, 1.0 / 10.0, 0, 0, 0, 2.0 * l / 15.0, 0, -1.0 / 10.0, 0, 0, 0, -l / 30.0,
    // 	-1.0 / l, 0, 0, 0, 0, 0, 1.0 / l, 0, 0, 0, 0, 0,
    // 	0, -6.0 / 15.0 / l, 0, 0, 0, -1.0 / 10.0, 0, 6.0 / 5.0 / l, 0, 0, 0, -1.0 / 10.0,
    // 	0, 0, -6.0 / 15.0 / l, 0, 1.0 / 10.0, 0, 0, 0, 6.0 / 5.0 / l, 0, 1.0 / 10.0, 0,
    // 	0, 0, 0, -1.0 / l, 0, 0, 0, 0, 0, 1.0 / l, 0, 0,
    // 	0, 0, -1.0 / 10.0, 0, -l / 30.0, 0, 0, 0, 1.0 / 10.0, 0, 2.0 * l / 15.0, 0,
    // 	0, 1.0 / 10.0, 0, 0, 0, -l / 30.0, 0, -1 / 10.0, 0, 0, 0, 2.0 * l / 15.0;
    // Kg_mat *= Nx;

    // A1 A1^t
    // Eigen::MatrixXd Kg_mat1 = Eigen::MatrixXd::Zero(total_dof, total_dof);
    // Kg_mat1 << Nx / l, -(s.S0.Mz - s.S1.Mz) / (l * l), (s.S0.My - s.S1.My) / (l * l), 0, s.S0.My / l, -s.S0.Mz / l, -Nx / l, (s.S1.Mz - s.S0.Mz) / (l * l), (s.S1.My - s.S0.My) / (l * l), 0, -s.S1.My / l, s.S1.Mz / l,
    //	-(s.S0.Mz - s.S1.Mz) / (l * l), (Nx * Sec->Iz * 12) / (Sec->A * l * l * l), 0, 0, 0, -(Nx * Sec->Iz * 6) / (Sec->A * l * l), -(s.S1.Mz - s.S0.Mz) / (l * l), -(-Nx * 12 * Sec->Iz) / (Sec->A * l * l * l), 0, 0, 0, -(Nx * 6 * Sec->Iz) / (Sec->A * l * l),
    //	(s.S0.My - s.S1.My) / (l * l), 0, (Nx * Sec->Iy * 12) / (Sec->A * l * l * l), 0, (-Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (s.S1.My - s.S0.My) / (l * l), 0, (-Nx * 12 * Sec->Iy) / (Sec->A * l * l * l), 0, (-Nx * 6 * Sec->Iy) / (Sec->A * l * l), 0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	s.S0.My / l, 0, (-Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (Nx * Sec->Iy * 4) / (Sec->A * l), 0, -s.S0.My / l, 0, (Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (Nx * Sec->Iy * 2) / (Sec->A * l), 0,
    //	-s.S0.Mz / l, -(Nx * Sec->Iz * 6) / (Sec->A * l * l), 0, 0, 0, (Nx * Sec->Iz * 4) / (Sec->A * l), s.S0.Mz / l, (-Nx * Sec->Iz * 6) / (Sec->A * l * l), 0, 0, 0, (Nx * Sec->Iz * 2) / (Sec->A * l),
    //	-Nx / l, -(s.S1.Mz - s.S0.Mz) / (l * l), (s.S1.My - s.S0.My) / (l * l), 0, -s.S0.My / l, s.S0.Mz / l, Nx / l, (s.S0.Mz - s.S1.Mz) / (l * l), (s.S0.My - s.S1.My) / (l * l), 0, s.S1.My / l, -s.S1.Mz / l,
    //	(s.S1.Mz - s.S0.Mz) / (l * l), (Nx * 12 * Sec->Iz) / (Sec->A * l * l * l), 0, 0, 0, (-Nx * Sec->Iz * 6) / (Sec->A * l * l), (s.S0.Mz - s.S1.Mz) / (l * l), (Nx * 12 * Sec->Iz) / (Sec->A * l * l * l), 0, 0, 0, (-Nx * 6 * Sec->Iz) / (Sec->A * l * l),
    //	(s.S1.My - s.S0.My) / (l * l), 0, (-Nx * 12 * Sec->Iy) / (Sec->A * l * l * l), 0, (Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (s.S0.My - s.S1.My) / (l * l), 0, (Nx * 12 * Sec->Iy) / (Sec->A * l * l * l), 0, (Nx * 6 * Sec->Iy) / (Sec->A * l * l), 0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	-s.S1.My / l, 0, (-Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (Nx * Sec->Iy * 2) / (Sec->A * l), 0, s.S1.My / l, 0, (Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (Nx * Sec->Iy * 4) / (Sec->A * l), 0,
    //	s.S1.Mz / l, -(Nx * Sec->Iz * 6) / (Sec->A * l * l), 0, 0, 0, (Nx * Sec->Iz * 2) / (Sec->A * l), -s.S1.Mz / l, (-Nx * Sec->Iz * 6) / (Sec->A * l * l), 0, 0, 0, (Nx * Sec->Iz * 4) / (Sec->A * l);

    // Eigen::MatrixXd Kg_mat1 = Eigen::MatrixXd::Zero(total_dof, total_dof);
    // Kg_mat1 << Nx / l, (s.S0.Mz - s.S1.Mz) / (l * l), (s.S0.My - s.S1.My) / (l * l), 0, s.S0.My / l, -s.S0.Mz / l, -Nx / l, (s.S1.Mz - s.S0.Mz) / (l * l), (s.S1.My - s.S0.My) / (l * l), 0, -s.S1.My / l, s.S1.Mz / l,
    //	(s.S0.Mz - s.S1.Mz) / (l * l), (Nx * Sec->Iz * 12) / (Sec->A * l * l * l), 0, 0, 0, (Nx * Sec->Iz * 6) / (Sec->A * l * l), (s.S1.Mz - s.S0.Mz) / (l * l), (-Nx * 12 * Sec->Iz) / (Sec->A * l * l * l), 0, 0, 0, (Nx * 6 * Sec->Iz) / (Sec->A * l * l),
    //	(s.S0.My - s.S1.My) / (l * l), 0, (Nx * Sec->Iy * 12) / (Sec->A * l * l * l), 0, (-Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (s.S1.My - s.S0.My) / (l * l), 0, (-Nx * 12 * Sec->Iy) / (Sec->A * l * l * l), 0, (-Nx * 6 * Sec->Iy) / (Sec->A * l * l), 0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	s.S0.My / l, 0, (-Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (Nx * Sec->Iy * 4) / (Sec->A * l), 0, -s.S0.My / l, 0, (Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (Nx * Sec->Iy * 2) / (Sec->A * l), 0,
    //	-s.S0.Mz / l, (Nx * Sec->Iz * 6) / (Sec->A * l * l), 0, 0, 0, (Nx * Sec->Iz * 4) / (Sec->A * l), s.S0.Mz / l, (-Nx * Sec->Iz * 6) / (Sec->A * l * l), 0, 0, 0, (Nx * Sec->Iz * 2) / (Sec->A * l),
    //	-Nx / l, (s.S1.Mz - s.S0.Mz) / (l * l), (s.S1.My - s.S0.My) / (l * l), 0, -s.S0.My / l, s.S0.Mz / l, Nx / l, (s.S0.Mz - s.S1.Mz) / (l * l), (s.S0.My - s.S1.My) / (l * l), 0, s.S1.My / l, -s.S1.Mz / l,
    //	(s.S1.Mz - s.S0.Mz) / (l * l), (-Nx * 12 * Sec->Iz) / (Sec->A * l * l * l), 0, 0, 0, (-Nx * Sec->Iz * 6) / (Sec->A * l * l), (s.S0.Mz - s.S1.Mz) / (l * l), (Nx * 12 * Sec->Iz) / (Sec->A * l * l * l), 0, 0, 0, (-Nx * 6 * Sec->Iz) / (Sec->A * l * l),
    //	(s.S1.My - s.S0.My) / (l * l), 0, (-Nx * 12 * Sec->Iy) / (Sec->A * l * l * l), 0, (Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (s.S0.My - s.S1.My) / (l * l), 0, (Nx * 12 * Sec->Iy) / (Sec->A * l * l * l), 0, (Nx * 6 * Sec->Iy) / (Sec->A * l * l), 0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	-s.S1.My / l, 0, (-Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (Nx * Sec->Iy * 2) / (Sec->A * l), 0, s.S1.My / l, 0, (Nx * Sec->Iy * 6) / (Sec->A * l * l), 0, (Nx * Sec->Iy * 4) / (Sec->A * l), 0,
    //	s.S1.Mz / l, (Nx * Sec->Iz * 6) / (Sec->A * l * l), 0, 0, 0, (Nx * Sec->Iz * 2) / (Sec->A * l), -s.S1.Mz / l, (-Nx * Sec->Iz * 6) / (Sec->A * l * l), 0, 0, 0, (Nx * Sec->Iz * 4) / (Sec->A * l);

    // A2 A2^t + A3 A3^t
    // Eigen::MatrixXd Kg_mat2 = Eigen::MatrixXd::Zero(total_dof, total_dof);
    // Kg_mat2 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	0, (6.0 * Nx) / (5.0 * l), 0, -(s.S0.My + s.S1.My) / (2.0 * l), 0, Nx / 10.0, 0, -(6.0 * Nx) / (5.0 * l), 0, (s.S0.My + s.S1.My) / (2.0 * l), 0, Nx / 10,
    //	0, 0, (6.0 * Nx) / (5.0 * l), (s.S0.Mz + s.S1.Mz) / (2.0 * l), -Nx / 10.0, 0, 0, 0, -(6.0 * Nx) / (5.0 * l), -(s.S0.Mz + s.S1.Mz) / (2.0 * l), -Nx / 10, 0,
    //	0, -(s.S0.My + s.S1.My) / (2.0 * l), (s.S0.Mz + s.S1.Mz) / (2.0 * l), Nx * (Sec->Iy + Sec->Iz) / (Sec->A * l * l), 0, 0, 0, (s.S0.My + s.S1.My) / (2.0 * l), -(s.S0.Mz + s.S1.Mz) / (2.0 * l), -Nx * (Sec->Iy + Sec->Iz) / (Sec->A * l * l), 0, 0,
    //	0, 0, -Nx / 10, 0, 2 * Nx * l / 15, 0, 0, 0, Nx / 10, 0, -Nx * l / 30, 0,
    //	0, Nx / 10, 0, 0, 0, 2.0 * Nx * l / 15, 0, -Nx / 10, 0, 0, 0, -Nx * l / 30,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	0, -6.0 * Nx / (5 * l), 0, (s.S0.My + s.S1.My) / (2 * l), 0, -Nx / 10, 0, 6.0 * Nx / (5.0 * l), 0, -(s.S0.My + s.S1.My) / (2.0 * l), 0, -Nx / 10,
    //	0, 0, -6.0 * Nx / (5.0 * l), -(s.S0.Mz + s.S1.Mz) / (2.0 * l), Nx / 10, 0, 0, 0, 6.0 * Nx / (5.0 * l), (s.S0.Mz + s.S1.Mz) / (2.0 * l), Nx / 10, 0,
    //	0, (s.S0.My + s.S1.My) / (2.0 * l), -(s.S0.Mz + s.S1.Mz) / (2.0 * l), -Nx * (Sec->Iy + Sec->Iz) / (Sec->A * l * l), 0, 0, 0, -(s.S0.My + s.S1.My) / (2.0 * l), (s.S0.Mz + s.S1.Mz) / (2.0 * l), Nx * (Sec->Iy + Sec->Iz) / (Sec->A * l * l), 0, 0,
    //	0, 0, -Nx / 10, 0, -Nx * l / 30.0, 0, 0, 0, Nx / 10.0, 0, 2.0 * Nx * l / 15.0, 0,
    //	0, -Nx / 10.0, 0, 0, 0, -Nx * l / 30.0, 0, -Nx / 10.0, 0, 0, 0, 2.0 * l * Nx / 15.0;

    // A6 A1^t + A1 A6^t + A7 A2^t + A2 A7^t
    // Eigen::MatrixXd Kg_mat3 = Eigen::MatrixXd::Zero(total_dof, total_dof);
    // Kg_mat3 << 0, 0, -1.0/l, 0, 0, 0, 0, 0, 1.0/l, 0, 0, 0,
    //	0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0.5, 0, 0,
    //	-1.0/l, 0, 0, 0, 0, 0, 1.0/l, 0, 0, 0, 0, 0,
    //	0, 0.5, 0, 0, 0, -l/12.0, 0, -0.5, 0, 0, 0, l/12.0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	0, 0, 0, -l/12.0, 0, 0, 0, 0, 0, l/12.0, 0, 0,
    //	0, 0, 1.0/l, 0, 0, 0, 0, 0, -1.0/l, 0, 0, 0,
    //	0, 0, 0, -0.5, 0, 0, 0, 0, 0, -0.5, 0, 0,
    //	1.0 / l, 0, 0, 0, 0, 0, -1.0 / l, 0, 0, 0, 0, 0,
    //	0, 0.5, 0, 0, 0, l / 12.0, 0, -0.5, 0, 0, 0, -l / 12.0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	0, 0, 0, l / 12.0, 0, 0, 0, 0, 0, -l / 12.0, 0, 0;
    // Kg_mat3 *= Qz;
    //
    // A1 A4^t + A4 A1^t + A3 A5^t + A5 A3^t
    // Eigen::MatrixXd Kg_mat4 = Eigen::MatrixXd::Zero(total_dof, total_dof);
    // Kg_mat4 << 0, -1.0/l, 0, 0, 0, 0, 0, 1.0/l, 0, 0, 0, 0,
    //	-1.0/l, 0, 0, 0, 0, 0, 1.0/l, 0, 0, 0, 0, 0,
    //	0, 0, 0, -0.5, 0, 0, 0, 0, 0, -0.5, 0, 0,
    //	0, 0, -0.5, 0, -l/12.0, 0, 0, 0, 0.5, 0, l/12.0, 0,
    //	0, 0, 0, -l/12.0, 0, 0, 0, 0, 0, l/12.0, 0, 0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    //	0, 1.0/l, 0, 0, 0, 0, 0, -1.0/l, 0, 0, 0, 0,
    //	1.0/l, 0, 0, 0, 0, 0, -1.0/l, 0, 0, 0, 0, 0,
    //	0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0.5, 0, 0,
    //	0, 0, -0.5, 0, l / 12.0, 0, 0, 0, 0.5, 0, -l / 12.0, 0,
    //	0, 0, 0, l / 12.0, 0, 0, 0, 0, 0, -l / 12.0, 0, 0,
    //	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    // Kg_mat4 *= Qy;

    // Kg_mat = Kg_mat1 + Kg_mat2 + Kg_mat3 + Kg_mat4;
    // Kg_mat *= 0.5;

    // Kg_mat << 1.0 / l, 0, 0, 0, 0, 0, -1.0 / l, 0, 0, 0, 0, 0,
    // 	0, 6.0 / 5.0 / l, 0, 0, 0, 1.0 / 10.0, 0, -6.0 / 5.0 / l, 0, 0, 0, 1.0 / 10.0,
    // 	0, 0, 6.0 / 5.0 / l, 0, -1.0 / 10.0, 0, 0, 0, -6.0 / 5.0 / l, 0, -1.0 / 10.0, 0,
    // 	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    // 	0, 0, -1.0 / 10.0, 0, 2.0 * l / 15.0, 0, 0, 0, 1.0 / 10.0, 0, -l / 30.0, 0,
    // 	0, 1.0 / 10.0, 0, 0, 0, 2.0 * l / 15.0, 0, -1.0 / 10.0, 0, 0, 0, -l / 30.0,
    // 	-1.0 / l, 0, 0, 0, 0, 0, 1.0 / l, 0, 0, 0, 0, 0,
    // 	0, -6.0 / 5.0 / l, 0, 0, 0, -1.0 / 10.0, 0, 6.0 / 5.0 / l, 0, 0, 0, -1.0 / 10.0,
    // 	0, 0, -6.0 / 5.0 / l, 0, 1.0 / 10.0, 0, 0, 0, 6.0 / 5.0 / l, 0, 1.0 / 10.0, 0,
    // 	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    // 	0, 0, -1.0 / 10.0, 0, -l / 30.0, 0, 0, 0, 1.0 / 10.0, 0, 2.0 * l / 15.0, 0,
    // 	0, 1.0 / 10.0, 0, 0, 0, -l / 30.0, 0, -1 / 10.0, 0, 0, 0, 2.0 * l / 15.0;

    Eigen::MatrixXd tr = trans_matrix();
    return tr.transpose() * Kg_mat * tr;
}

Eigen::MatrixXd BeamElement::trans_matrix()
{

    Eigen::Matrix3d tr0 = trans_matrix3(Nodes[0]->Location, Nodes[1]->Location, Beta);

    // ブロックの数を指定
    const int numBlocks = 4;

    // 大きな行列を作成し、対角ブロック行列として同じ行列を配置
    Eigen::MatrixXd matrix(total_dof, total_dof);
    matrix.setZero(); // 行列を0で初期化

    for (int i = 0; i < numBlocks; ++i)
    {
        matrix.block(i * tr0.rows(), i * tr0.cols(), tr0.rows(), tr0.cols()) = tr0;
    }

    return matrix;
}

// Axis
Eigen::Matrix2d stiffness_matrix_truss(double E, double A, double L)
{
    Eigen::Matrix2d matrix;
    double EA_L = E * A / L;

    matrix << EA_L, -EA_L,
        -EA_L, EA_L;
    return matrix;
}

// Beam Torsion
Eigen::Matrix2d stiffness_matrix_beam_rot_x(double G, double K, double L)
{
    Eigen::Matrix2d matrix;
    double GKl = G * K / L;

    matrix << GKl, -GKl,
        -GKl, GKl;
    return matrix;
}

// Beam Rotation around z
Eigen::Matrix4d stiffness_matrix_beam_rot_z(double E, double Iz, double L)
{
    Eigen::Matrix4d matrix;
    double EI = E * Iz;
    double L2 = L * L;
    double L3 = L2 * L;

    matrix << 12.0 / L3, 6.0 / L2, -12.0 / L3, 6.0 / L2,
        6.0 / L2, 4.0 / L, -6.0 / L2, 2.0 / L,
        -12.0 / L3, -6.0 / L2, 12.0 / L3, -6.0 / L2,
        6.0 / L2, 2.0 / L, -6.0 / L2, 4.0 / L;
    return EI * matrix;
}

// Beam Rotation around y
Eigen::Matrix4d stiffness_matrix_beam_rot_y(double E, double Iy, double L)
{
    Eigen::Matrix4d matrix;
    double EI = E * Iy;
    double L2 = L * L;
    double L3 = L2 * L;

    matrix << 12.0 / L3, -6.0 / L2, -12.0 / L3, -6.0 / L2,
        -6.0 / L2, 4.0 / L, 6.0 / L2, 2.0 / L,
        -12.0 / L3, 6.0 / L2, 12.0 / L3, 6.0 / L2,
        -6.0 / L2, 2.0 / L, 6.0 / L2, 4.0 / L;

    return EI * matrix;
}


// compute_Kprime_partial 関数
// 
// 引数
//   Kz        : 4x4 の対称行列 (Eigen::Matrix4d)
//   lambda_s  : float (double)
//   lambda_z  : float (double)
//   lambda_s_ : float (double)  // λ_s'
//   lambda_z_ : float (double)  // λ_z'
//
// 戻り値
//   Kp        : 4x4 の対称行列 (Eigen::Matrix4d)
//
Eigen::Matrix4d compute_Kprime_partial(const Eigen::Matrix4d& Kz,
	double lambda_s, double lambda_z,
	double lambda_s_, double lambda_z_)
{
	// --- 1) Kz の各要素を取り出す ---
	double k11 = Kz(0, 0);
	double k12 = Kz(0, 1);
	double k13 = Kz(0, 2);
	double k14 = Kz(0, 3);
	double k22 = Kz(1, 1);
	double k23 = Kz(1, 2);
	double k24 = Kz(1, 3);
	double k33 = Kz(2, 2);
	double k34 = Kz(2, 3);
	double k44 = Kz(3, 3);

	// --- 2) Lambda|D| の計算 ---
	double LambdaD =
		k11 * k22 * k33 * k44
		+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
		- k11 * (k24 * k24) * k33 * (1.0 - lambda_z) * (1.0 - lambda_z_)
		- k11 * (k23 * k23) * k44 * (1.0 - lambda_s_) * (1.0 - lambda_z)
		- k11 * k22 * (k34 * k34) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
		- (k12 * k12) * k33 * k44 * (1.0 - lambda_s) * (1.0 - lambda_z)
		- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
		- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
		+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_s) * (1.0 - lambda_z) * (1.0 - lambda_z_)
		+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z)
		+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
		+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
		+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
		- (k13 * k13) * k22 * k44 * (1.0 - lambda_s) * (1.0 - lambda_s_)
		- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
		- (k14 * k14) * k22 * k33 * (1.0 - lambda_s) * (1.0 - lambda_z_)
		+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_);

	// --- 3) K' 用の 4x4 行列を 0 で初期化 ---
	Eigen::Matrix4d Kp = Eigen::Matrix4d::Zero();

	// ---------------------------------------------------------------------
	// 上三角 (i <= j) の要素を計算し，下三角へコピー (対称行列を構築)
	// ---------------------------------------------------------------------

	// === (A) K_{11}' ===
	Kp(0, 0) =
		((lambda_s * k11) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_z) * (1.0 - lambda_z_)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_s_) * (1.0 - lambda_z)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_z)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_z) * (1.0 - lambda_z_)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_s_) * (1.0 - lambda_z)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_s_)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_z_)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_s_) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			);

	// === (B) K_{12}' ===
	Kp(0, 1) =
		-(lambda_s * lambda_z * k11 * k22 / LambdaD) * (
			-k12 * k33 * k44
			- k13 * k34 * k24 * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			- k14 * k23 * k34 * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			+ k14 * k33 * k24 * (1.0 - lambda_z_)
			+ k13 * k23 * k44 * (1.0 - lambda_s_)
			+ k12 * (k34 * k34) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			);

	// === (C) K_{13}' ===
	Kp(0, 2) =
		-(lambda_s * lambda_s_ * k11 * k33 / LambdaD) * (
			(1.0 - lambda_z) * k12 * k23 * k44
			+ (1.0 - lambda_z) * (1.0 - lambda_z_) * k13 * (k24 * k24)
			+ (1.0 - lambda_z_) * k14 * k22 * k34
			- (1.0 - lambda_z) * (1.0 - lambda_z_) * k14 * k23 * k24
			- k13 * k22 * k44
			- (1.0 - lambda_z) * (1.0 - lambda_z_) * k12 * k24 * k34
			);

	// === (D) K_{14}' ===
	Kp(0, 3) =
		-(lambda_s * lambda_z_ * k11 * k44 / LambdaD) * (
			-(1.0 - lambda_s_) * (1.0 - lambda_z) * k12 * k23 * k34
			- (1.0 - lambda_s_) * (1.0 - lambda_z) * k13 * k24 * k23
			- k14 * k22 * k33
			+ (1.0 - lambda_s_) * (1.0 - lambda_z) * k14 * (k23 * k23)
			+ (1.0 - lambda_s_) * k13 * k22 * k34
			+ (1.0 - lambda_z) * k12 * k24 * k33
			);

	// === (E) K_{22}' ===
	Kp(1, 1) =
		((lambda_z * k22) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_z_)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_s_)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_s)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_s) * (1.0 - lambda_z_)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_s) * (1.0 - lambda_s_)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_s) * (1.0 - lambda_s_)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_s) * (1.0 - lambda_z_)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z_)
			);

	// === (F) K_{23}' ===
	Kp(1, 2) =
		-(lambda_z * lambda_s_ * k22 * k33 / LambdaD) * (
			-k11 * k23 * k44
			- (1.0 - lambda_s) * (1.0 - lambda_z_) * k13 * k24 * k14
			- (1.0 - lambda_s) * (1.0 - lambda_z_) * k14 * k12 * k34
			+ (1.0 - lambda_s) * (1.0 - lambda_z_) * (k14 * k14) * k23
			+ (1.0 - lambda_s) * k13 * k12 * k44
			+ (1.0 - lambda_z_) * k11 * k24 * k34
			);

	// === (G) K_{24}' ===
	Kp(1, 3) =
		-(lambda_z * lambda_z_ * k22 * k44 / LambdaD) * (
			(1.0 - lambda_s_) * k11 * k23 * k34
			+ (1.0 - lambda_s) * (1.0 - lambda_s_) * (k13 * k13) * k24
			+ (1.0 - lambda_s) * k14 * k12 * k33
			- (1.0 - lambda_s) * (1.0 - lambda_s_) * k14 * k23 * k13
			- (1.0 - lambda_s) * (1.0 - lambda_s_) * k13 * k12 * k34
			- k11 * k24 * k33
			);

	// === (H) K_{33}' ===
	Kp(2, 2) =
		((lambda_s_ * k33) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_z) * (1.0 - lambda_z_)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_z) * (1.0 - lambda_z_)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_z)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_z_)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_s) * (1.0 - lambda_z)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_s) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_s) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_s) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_s) * (1.0 - lambda_z)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_s) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_s) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_s) * (1.0 - lambda_z_)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_s)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_s) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_s) * (1.0 - lambda_z_)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_s) * (1.0 - lambda_z) * (1.0 - lambda_z_)
			);

	// === (I) K_{34}' ===
	Kp(2, 3) =
		-(lambda_s_ * lambda_z_ * k33 * k44 / LambdaD) * (
			-k11 * k22 * k34
			- (1.0 - lambda_s) * (1.0 - lambda_z) * k12 * k24 * k13
			- (1.0 - lambda_s) * (1.0 - lambda_z) * k14 * k12 * k23
			+ (1.0 - lambda_s) * k14 * k22 * k13
			+ (1.0 - lambda_s) * (1.0 - lambda_z) * (k12 * k12) * k34
			+ (1.0 - lambda_z) * k11 * k24 * k23
			);

	// === (J) K_{44}' ===
	Kp(3, 3) =
		((lambda_z_ * k44) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_s_) * (1.0 - lambda_z)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_z)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_s_) * (1.0 - lambda_z)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_s_)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_s) * (1.0 - lambda_z)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_s) * (1.0 - lambda_z)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_s) * (1.0 - lambda_s_)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_s) * (1.0 - lambda_s_)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_s)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_s) * (1.0 - lambda_s_) * (1.0 - lambda_z)
			);

	// --- 4) 上三角を計算したので，下三角へコピーして対称行列を完成 ---
	for (int i = 0; i < 4; i++)
	{
		for (int j = i + 1; j < 4; j++)
		{
			Kp(j, i) = Kp(i, j);
		}
	}

	return Kp;
}


Eigen::MatrixXd ComplexBeamElement::stiffness_matrix_local()
{
    double l = length();
    double lz_ = l - lzi - lzj;
    double ly_ = l - lyi - lyj;

    Eigen::Matrix4d Kbz = stiffness_matrix_beam_rot_z(Mat.Young, Sec->Iz, lz_);
    Eigen::Matrix4d Tbz = Eigen::Matrix4d::Identity();
    Tbz(0, 1) = lzi;
    Tbz(2, 3) = -lzj;

    Eigen::Matrix4d Kby = stiffness_matrix_beam_rot_y(Mat.Young, Sec->Iy, ly_);
    Eigen::Matrix4d Tby = Eigen::Matrix4d::Identity();
    Tby(0, 1) = -lzi;
    Tby(2, 3) = lzj;

    Eigen::Matrix2d Kx = stiffness_matrix_truss(Mat.Young, Sec->A, l);

    Eigen::Matrix2d Kt = stiffness_matrix_beam_rot_x(Mat.G(), Sec->K, l);

    // 端部バネの考慮
    Kbz = compute_Kprime_partial(Kbz, Lambda_sy, Lambda_bz, Lambda_sy_, Lambda_bz_);
    Kby = compute_Kprime_partial(Kby, Lambda_sz, Lambda_by, Lambda_sz_, Lambda_by_);

    // 剛域の考慮
    Kbz = Tbz.transpose() * Kbz * Tbz;
    Kby = Tby.transpose() * Kby * Tby;

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(12, 12);
    K(0, 0) = Kx(0, 0);
    K(0, 6) = Kx(0, 1);
    K(6, 0) = Kx(1, 0);
    K(6, 6) = Kx(1, 1);

    K(1, 1) = Kbz(0, 0);
    K(1, 5) = Kbz(0, 1);
    K(1, 7) = Kbz(0, 2);
    K(1, 11) = Kbz(0, 3);
    K(5, 1) = Kbz(1, 0);
    K(5, 5) = Kbz(1, 1);
    K(5, 7) = Kbz(1, 2);
    K(5, 11) = Kbz(1, 3);
    K(7, 1) = Kbz(2, 0);
    K(7, 5) = Kbz(2, 1);
    K(7, 7) = Kbz(2, 2);
    K(7, 11) = Kbz(2, 3);
    K(11, 1) = Kbz(3, 0);
    K(11, 5) = Kbz(3, 1);
    K(11, 7) = Kbz(3, 2);
    K(11, 11) = Kbz(3, 3);

    K(2, 2) = Kby(0, 0);
    K(2, 4) = Kby(0, 1);
    K(2, 8) = Kby(0, 2);
    K(1, 10) = Kby(0, 3);
    K(4, 2) = Kby(1, 0);
    K(4, 4) = Kby(1, 1);
    K(4, 8) = Kby(1, 2);
    K(4, 10) = Kby(1, 3);
    K(8, 2) = Kby(2, 0);
    K(8, 4) = Kby(2, 1);
    K(8, 8) = Kby(2, 2);
    K(8, 10) = Kby(2, 3);
    K(10, 2) = Kby(3, 0);
    K(10, 4) = Kby(3, 1);
    K(10, 8) = Kby(3, 2);
    K(10, 10) = Kby(3, 3);

    K(3, 3) = Kt(0, 0);
    K(3, 9) = Kt(0, 1);
    K(9, 3) = Kt(1, 0);
    K(9, 9) = Kt(1, 1);

    return K;
}

ComplexBeamElement::ComplexBeamElement()
{
    Lambda_bz = 1;
    Lambda_bz_ = 1;
    Lambda_sy = 1;
    Lambda_sy_ = 1;
    Lambda_by = 1;
    Lambda_by_ = 1;
    Lambda_sz = 1;
    Lambda_sz_ = 1;

    lzi = 0;
    lzj = 0;
    lyi = 0;
    lyj = 0;
}

ComplexBeamElement::ComplexBeamElement(Node *n0, Node *n1, Section *sec, Material mat, double beta)
    : BeamElement(n0, n1, sec, mat, beta)
{
    Lambda_bz = 1;
    Lambda_bz_ = 1;
    Lambda_sy = 1;
    Lambda_sy_ = 1;
    Lambda_by = 1;
    Lambda_by_ = 1;
    Lambda_sz = 1;
    Lambda_sz_ = 1;

    lzi = 0;
    lzj = 0;
    lyi = 0;
    lyj = 0;
}

ComplexBeamElement::ComplexBeamElement(int _id, Node *n0, Node *n1, Section *sec, Material mat, double beta)
    : ComplexBeamElement(n0, n1, sec, mat, beta)
{
    id = _id;

    Lambda_bz = 1;
    Lambda_bz_ = 1;
    Lambda_sy = 1;
    Lambda_sy_ = 1;
    Lambda_by = 1;
    Lambda_by_ = 1;
    Lambda_sz = 1;
    Lambda_sz_ = 1;

    lzi = 0;
    lzj = 0;
    lyi = 0;
    lyj = 0;
}

Eigen::MatrixXd ComplexBeamElement::StiffnessMatrix()
{
    Eigen::MatrixXd tr = trans_matrix();
    Eigen::MatrixXd k = stiffness_matrix_local();
    return tr.transpose() * k * tr;
}

// Eigen::MatrixXd ComplexBeamElement::NodeConsistentMass()
//{
//	Eigen::MatrixXd tr = trans_matrix();
//	Eigen::MatrixXd m = Eigen::MatrixXd::Zero(total_dof, total_dof);
//	double l = element_length();
//	double l2 = l * l;
//	double l3 = l2 * l;
//	m << l / 3, 0, 0, 0, 0, 0, l / 6, 0, 0, 0, 0, 0,
//		0, 13 / 35 * l, 0, 0, 0, 11 / 210 * l2, 0, 9 / 70 * l, 0, 0, 0, -13 / 420 * l2,
//		0, 0, 13 / 35 * l, 0, -11 / 210 * l2, 0, 0, 0, 9 / 70 * l, 0, 13 / 420 * l2, 0,
//		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//		0, 0, -11 / 210 * l2, 0, 1 / 105 * l3, 0, 0, 0, -13 / 420 * l2, 0, -1 / 140 * l3, 0,
//		0, 11 / 210 * l2, 0, 0, 0, 1 / 105 * l3, 0, 13 / 420 * l2, 0, 0, 0, -1 / 140 * l3,
//		l / 6, 0, 0, 0, 0, 0, l / 3, 0, 0, 0, 0, 0,
//		0, 9 / 70 * l, 0, 0, 0, 13 / 420 * l2, 0, 13 / 35 * l, 0, 0, 0, -11 / 210 * l2,
//		0, 0, 9 / 70 * l, 0, -13 / 420 * l2, 0, 0, 0, 13 / 35 * l, 0, 11 / 210 * l2, 0,
//		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//		0, 0, 13 / 420 * l2, 0, -1 / 140 * l3, 0, 0, 0, 11 / 210 * l2, 0, 1 / 105 * l3, 0,
//		0, -13 / 420 * l2, 0, 0, 0, -1 / 140 * l3, 0, -11 / 210 * l2, 0, 0, 0, 1 / 105 * l3;
//
//	return tr.transpose() * (Sec->A * Mat.dense) * m * tr;
// }

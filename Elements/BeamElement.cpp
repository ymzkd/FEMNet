#include "BeamElement.h"

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
    double N = s.S0.Nx;

	double Mzi = s.S0.Mz;
	double Mzj = s.S1.Mz;
	double Myi = s.S0.My;
	double Myj = s.S1.My;

    double Qz = (Myj + Myi) / l;
    double Qy = -(Mzj + Mzi) / l;

	double R = (Sec->Iy + Sec->Iz) / Sec->A;

    int indices_b[10]{1, 2, 3, 4, 5, 7, 8, 9, 10, 11};
    int indices_x[2]{0, 6};

    // 曲げによる幾何剛性
    Eigen::MatrixXd Kg_b = Eigen::MatrixXd::Zero(10, 10);
    Kg_b << 6.0*N/(5.0*l), 0, Myi/l, 0, -N/10.0, -6.0*N/(5.0*l), 0, Myj/l, 0, N/10.0,
        0, 6.0*N/(5.0*l), Mzi/l, N/10.0, 0, 0, -6.0*N/(5.0*l), Mzj/l, N*l/30.0, 0,
        Myi/l, Mzi/l, N*R/l, -Qy*l/6.0, -Qz*l/6.0, Myj/l, -Mzi/l, -N*R/l, Qy*l/6.0, Qz*l/6.0,
        0, N/10.0, -Qy*l/6.0, 2.0*l*N/15.0, 0, 0, -N/10.0, Qy*l/6.0, l*N/30.0, 0,
        -N/10.0, 0, -Qz*l/6.0, 0, 2.0*l*N/15.0, N/10.0, 0, Qz*l/6.0, 0, N*l/30.0,
        -6.0*N/(5.0*l), 0, Myj/l, 0, N/10.0, 6.0*N/(5.0*l), 0, Myi/l, 0, -N/10.0,
        0, -6.0*N/(5.0*l), -Mzi/l, -N/10.0, 0, 0, 6.0*N/(5.0*l), -Mzj/l, N/10.0, 0,
        Myj/l, Mzj/l, -N*R/l, Qy*l/6.0, Qz*l/6.0, Myi/l, -Mzj/l, N*R/l, -Qy*l/6.0, -Qz*l/6.0,
        0, N*l/30.0, Qy*l/6.0, l*N/30.0, 0, 0, N/10.0, -Qy*l/6.0, 2.0*l*N/15.0, 0,
        N/10.0, 0, Qz*l/6.0, 0, N*l/30.0, -N/10.0, 0, -Qz*l/6.0, 0, 2.0*l*N/15.0;
    
    for (size_t i = 0; i < 8; i++)
        for (size_t j = 0; j < 8; j++)
            Kg_mat(indices_b[i], indices_b[j]) += Kg_b(i, j);
    
    // 軸方向幾何剛性
    Eigen::MatrixXd Kg_x = Eigen::MatrixXd::Zero(2, 2);
    Kg_x << 1, -1,
        -1, 1;
    Kg_x *= N / l;

    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < 2; j++)
            Kg_mat(indices_x[i], indices_x[j]) += Kg_x(i, j);

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

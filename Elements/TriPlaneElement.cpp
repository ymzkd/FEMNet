#include "TriPlaneElement.h"

TriPlaneElement::TriPlaneElement(Node *n0, Node *n1, Node *n2,
                                 double t, Material mat, double beta)
{
    Nodes[0] = n0;
    Nodes[1] = n1;
    Nodes[2] = n2;

    thickness = Thickness(t);
    Mat = mat;
    plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
    plane.Rotate(beta, plane.ez);
}

TriPlaneElement::TriPlaneElement(Node *n0, Node *n1, Node *n2, Thickness t, Material mat, double beta)
{
    Nodes[0] = n0;
    Nodes[1] = n1;
    Nodes[2] = n2;

    thickness = t;
    Mat = mat;
    plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
    plane.Rotate(beta, plane.ez);
}

// Eigen::Matrix2d TriPlaneElement::invJMatrix()
//{
//	Point p1 = Nodes[0]->Location;
//	Point p2 = Nodes[1]->Location;
//	Point p3 = Nodes[2]->Location;
//	Eigen::Matrix2d mat;
//	double jd = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
//	mat(0, 0) = (p3.y - p1.y) / jd;
//	mat(0, 1) = (p1.y - p2.y) / jd;
//	mat(1, 0) = (p1.x - p3.x) / jd;
//	mat(1, 1) = (p2.x - p1.x) / jd;
//	return mat;
// }

std::vector<NodeLoadData> TriPlaneElement::AreaForceToNodeLoadData(std::vector<Vector> load_vecs)
{
    Eigen::Vector3d p1 = load_vecs[0].toEigen();
    Eigen::Vector3d p2 = load_vecs[1].toEigen();
    Eigen::Vector3d p3 = load_vecs[2].toEigen();

    double area = Area();
    NodeLoadData nl1(Nodes[0]->id,
                     area * (p1.x() / 6 + p2.x() / 12 + p3.x() / 12),
                     area * (p1.y() / 6 + p2.y() / 12 + p3.y() / 12),
                     area * (p1.z() / 6 + p2.z() / 12 + p3.z() / 12));

    NodeLoadData nl2(Nodes[1]->id,
                     area * (p1.x() / 12 + p2.x() / 6 + p3.x() / 12),
                     area * (p1.y() / 12 + p2.y() / 6 + p3.y() / 12),
                     area * (p1.z() / 12 + p2.z() / 6 + p3.z() / 12));

    NodeLoadData nl3(Nodes[2]->id,
                     area * (p1.x() / 12 + p2.x() / 12 + p3.x() / 6),
                     area * (p1.y() / 12 + p2.y() / 12 + p3.y() / 6),
                     area * (p1.z() / 12 + p2.z() / 12 + p3.z() / 6));

    std::vector<NodeLoadData> loads;
    loads.push_back(nl1);
    loads.push_back(nl2);
    loads.push_back(nl3);

    return loads;
}

std::vector<NodeLoadData> TriPlaneElement::InertialForceToNodeLoadData(Eigen::Vector3d accel_vec)
{
    Eigen::Matrix<double, total_dof, total_dof> m;
    m = NodeConsistentMass();
    Eigen::Vector<double, total_dof> f;
    f.setZero();

    f << accel_vec.x(), accel_vec.y(), accel_vec.z(),
        accel_vec.x(), accel_vec.y(), accel_vec.z(),
        accel_vec.x(), accel_vec.y(), accel_vec.z();

    f = m * f;
    std::vector<NodeLoadData> loads;
    loads.push_back(NodeLoadData(Nodes[0]->id, f(0), f(1), f(2)));
    loads.push_back(NodeLoadData(Nodes[1]->id, f(3), f(4), f(5)));
    loads.push_back(NodeLoadData(Nodes[2]->id, f(6), f(7), f(8)));
    return loads;
}

Eigen::MatrixXd TriPlaneElement::NodeConsistentMass()
{
    return Eigen::MatrixXd::Identity(total_dof, total_dof);
}

double TriPlaneElement::Area()
{
    Vector v01 = Nodes[1]->Location - Nodes[0]->Location;
    Vector v02 = Nodes[2]->Location - Nodes[0]->Location;
    Vector vn = Vector::cross(v01, v02);
    return vn.norm() / 2;
}

Eigen::MatrixXd TriPlaneElement::StiffnessMatrix()
{
    Eigen::MatrixXd K = localStiffnessMatrix();
    Eigen::MatrixXd TrMat = trans_matrix();
    return TrMat.transpose() * K * TrMat;
}

void TriPlaneElement::AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K)
{
    int indices[total_dof];
    for (size_t i = 0; i < node_dof; i++)
        indices[i] = Nodes[0]->id * 6 + i;
    for (size_t i = 0; i < node_dof; i++)
        indices[i + node_dof] = Nodes[1]->id * 6 + i;
    for (size_t i = 0; i < node_dof; i++)
        indices[i + node_dof * 2] = Nodes[2]->id * 6 + i;

    // Eigen::MatrixXd smat = StiffnessMatrix();
    for (size_t i = 0; i < total_dof; i++)
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

void TriPlaneElement::AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat)
{
    AssembleMatrix(mat, StiffnessMatrix());
    // int indices[total_dof];
    // for (size_t i = 0; i < node_dof; i++)
    // 	indices[i] = Nodes[0]->id * 6 + i;
    // for (size_t i = 0; i < node_dof; i++)
    // 	indices[i + node_dof] = Nodes[1]->id * 6 + i;
    // for (size_t i = 0; i < node_dof; i++)
    // 	indices[i + node_dof * 2] = Nodes[2]->id * 6 + i;

    // Eigen::MatrixXd smat = StiffnessMatrix();
    // for (size_t i = 0; i < total_dof; i++)
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

void TriPlaneElement::AssembleGeometricStiffMatrix(
    Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp)
{
    AssembleMatrix(mat, geometric_local_stiffness_matrix(disp));
}

void TriPlaneElement::AssembleMassMatrix(Eigen::SparseMatrix<double> &mat)
{
    Eigen::VectorXd mass = NodeLumpedMass();
    for (size_t ni = 0; ni < node_num; ni++)
    {
        for (size_t i = 0; i < 3; i++)
        {
            int idx = Nodes[ni]->id * 6 + i;
            mat.coeffRef(idx, idx) += mass[ni];
        }
    }
}

MembraneStressData TriPlaneElement::stress(
    Displacement d0, Displacement d1, Displacement d2)
{
    Eigen::VectorXd wvec(total_dof);
    wvec.segment(0, node_dof) = Eigen::Map<Eigen::VectorXd>(d0.displace, node_dof);
    wvec.segment(3, node_dof) = Eigen::Map<Eigen::VectorXd>(d1.displace, node_dof);
    wvec.segment(6, node_dof) = Eigen::Map<Eigen::VectorXd>(d2.displace, node_dof);

    wvec = trans_matrix() * wvec;
    Eigen::Vector3d strs = DMatrix() * BMatrix() * wvec;

    return MembraneStressData(strs(0), strs(1), strs(2));
}

Eigen::MatrixXd TriPlaneElement::BMatrix()
{
    Eigen::MatrixXd mat(3, 6);
    Point p1 = plane.PointToCoord(Nodes[0]->Location);
    Point p2 = plane.PointToCoord(Nodes[1]->Location);
    Point p3 = plane.PointToCoord(Nodes[2]->Location);

    mat << p2.y - p3.y, 0, p3.y - p1.y, 0, p1.y - p2.y, 0,
        0, p3.x - p2.x, 0, p1.x - p3.x, 0, p2.x - p1.x,
        p3.x - p2.x, p2.y - p3.y, p1.x - p3.x, p3.y - p1.y, p2.x - p1.x, p1.y - p2.y;
    return mat / (2 * Area());
}

Eigen::Matrix3d TriPlaneElement::DMatrix()
{
    Eigen::Matrix3d mat;
    mat << 1, Mat.Poisson, 0,
        Mat.Poisson, 1, 0,
        0, 0, (1 - Mat.Poisson) / 2;
    return Mat.Young / (1 - Mat.Poisson * Mat.Poisson) * mat;
}

Eigen::MatrixXd TriPlaneElement::localStiffnessMatrix()
{
    Eigen::MatrixXd BMat = BMatrix();
    Eigen::MatrixXd DMat = DMatrix();

    Eigen::MatrixXd K = thickness.plane_thick * Area() * BMat.transpose() * DMat * BMat;
    return K;
}

Eigen::MatrixXd TriPlaneElement::geometric_local_stiffness_matrix(const std::vector<Displacement> &disp)
{
    // Eigen::MatrixXd Kg_mat = Eigen::MatrixXd::Zero(total_dof, total_dof);
    double area = Area();

    Eigen::MatrixXd Gmat(6, 9);
    Point p1 = plane.PointToCoord(Nodes[0]->Location);
    Point p2 = plane.PointToCoord(Nodes[1]->Location);
    Point p3 = plane.PointToCoord(Nodes[2]->Location);

    Gmat << p2.y - p3.y, 0, 0, p3.y - p1.y, 0, 0, p1.y - p2.y, 0, 0,
        p3.x - p2.x, 0, 0, p1.x - p3.x, 0, 0, p2.x - p1.x, 0, 0,
        0, p2.y - p3.y, 0, 0, p3.y - p1.y, 0, 0, p1.y - p2.y, 0,
        0, p3.x - p2.x, 0, 0, p1.x - p3.x, 0, 0, p2.x - p1.x, 0,
        0, 0, p2.y - p3.y, 0, 0, p3.y - p1.y, 0, 0, p1.y - p2.y,
        0, 0, p3.x - p2.x, 0, 0, p1.x - p3.x, 0, 0, p2.x - p1.x;
    Gmat /= (2 * area);

    MembraneStressData strs = stress(disp[0], disp[1], disp[2]);
    Eigen::Matrix2d SigMat;
    SigMat << strs.sigx, strs.sigxy,
        strs.sigxy, strs.sigy;

    // 2x2ブロックを3つの対角位置に配置
    Eigen::MatrixXd SigMat3 = Eigen::MatrixXd::Zero(6, 6);
    SigMat3.block<2, 2>(0, 0) = SigMat;
    SigMat3.block<2, 2>(2, 2) = SigMat;
    SigMat3.block<2, 2>(4, 4) = SigMat;

    Eigen::MatrixXd Kg = Gmat.transpose() * SigMat3 * Gmat * area * thickness.plane_thick;
    // Eigen::MatrixXd tr = trans_matrix();
    // return tr.transpose() * Kg * tr;
    return Kg;
}

Eigen::MatrixXd TriPlaneElement::trans_matrix()
{
    Eigen::Matrix3d tr0 = trans_matrix3(plane);

    // ブロックの数を指定
    const int numBlocks = node_num;

    // 大きな行列を作成し、対角ブロック行列として同じ行列を配置
    Eigen::MatrixXd matrix(3 * 2, total_dof);
    matrix.setZero(); // 行列を0で初期化

    for (int i = 0; i < numBlocks; ++i)
    {
        matrix.block(i * 2, i * node_dof, 2, node_dof) = tr0.block(0, 0, 2, node_dof);
    }

    return matrix;
}

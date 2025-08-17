#include "QuadPlaneElement.h"

Eigen::Matrix2d QuadPlaneElement::JMatrix(double xi, double eta)
{
    Point p1 = plane.PointToCoord(Nodes[0]->Location);
    Point p2 = plane.PointToCoord(Nodes[1]->Location);
    Point p3 = plane.PointToCoord(Nodes[2]->Location);
    Point p4 = plane.PointToCoord(Nodes[3]->Location);

    Eigen::VectorXd dndxi(4), dndeta(4);
    dndxi << -0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1 + eta);
    dndeta << -0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1 - xi);
    // double dxdxi, dxdeta, dydxi, dydeta;

    Eigen::Matrix2d JMat;
    JMat << dndxi(0) * p1.x + dndxi(1) * p2.x + dndxi(2) * p3.x + dndxi(3) * p4.x,
        dndxi(0) * p1.y + dndxi(1) * p2.y + dndxi(2) * p3.y + dndxi(3) * p4.y,
        dndeta(0) * p1.x + dndeta(1) * p2.x + dndeta(2) * p3.x + dndeta(3) * p4.x,
        dndeta(0) * p1.y + dndeta(1) * p2.y + dndeta(2) * p3.y + dndeta(3) * p4.y;

    return JMat;
}

Eigen::MatrixXd QuadPlaneElement::BMatrix(double xi, double eta)
{
    Eigen::VectorXd dndxi(4), dndeta(4);
    dndxi << -0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1 + eta);
    dndeta << -0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1 - xi);

    // Eigen::Matrix2d JMat = JMatrix(xi, eta);
    // double detJ = JMat.determinant();
    // JMat = JMat.inverse();

    Eigen::MatrixXd dndx(2, 4);
    dndx.row(0) = dndxi;
    dndx.row(1) = dndeta;
    dndx = JMatrix(xi, eta).inverse() * dndx;

    Eigen::MatrixXd BMat(3, 8);
    BMat << dndx(0, 0), 0, dndx(0, 1), 0, dndx(0, 2), 0, dndx(0, 3), 0,
        0, dndx(1, 0), 0, dndx(1, 1), 0, dndx(1, 2), 0, dndx(1, 3),
        dndx(1, 0), dndx(0, 0), dndx(1, 1), dndx(0, 1), dndx(1, 2), dndx(0, 2), dndx(1, 3), dndx(0, 3);

    return BMat;
}

Eigen::Matrix3d QuadPlaneElement::DMatrix()
{
    Eigen::Matrix3d mat;
    mat << 1, Mat.Poisson, 0,
        Mat.Poisson, 1, 0,
        0, 0, (1 - Mat.Poisson) / 2.0;
    return Mat.Young / (1.0 - Mat.Poisson * Mat.Poisson) * mat;
}

Eigen::MatrixXd QuadPlaneElement::localStiffnessMatrix()
{
    const double intgp = 1.0 / sqrt(3.0);
    Eigen::Vector4d xi_list(-intgp, intgp, intgp, -intgp);
    Eigen::Vector4d eta_list(-intgp, -intgp, intgp, intgp);

    Eigen::MatrixXd K(total_dof_local, total_dof_local);
    K.setZero();

    Eigen::Matrix3d DMat = DMatrix();
    for (size_t i = 0; i < node_num; i++)
    {
        double detJ = JMatrix(xi_list(i), eta_list(i)).determinant();
        Eigen::MatrixXd BMat = BMatrix(xi_list(i), eta_list(i));
        K += thickness.plane_thick * BMat.transpose() * DMat * BMat * detJ;
    }
    return K;
}

Eigen::MatrixXd QuadPlaneElement::geometric_local_stiffness_matrix(const std::vector<Displacement> &disp)
{
    const double intgp = 1.0 / sqrt(3.0);
    Eigen::Vector4d xi_list(-intgp, intgp, intgp, -intgp);
    Eigen::Vector4d eta_list(-intgp, -intgp, intgp, intgp);
    Eigen::MatrixXd Kg = Eigen::MatrixXd::Zero(total_dof, total_dof);

    for (size_t i = 0; i < 4; i++)
    {
        double xi = xi_list(i);
        double eta = eta_list(i);

        Eigen::VectorXd dndxi(4), dndeta(4);
        dndxi << -0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1 + eta);
        dndeta << -0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1 - xi);

        Eigen::Matrix2d JMat = JMatrix(xi, eta);
        double detJ = JMat.determinant();
        // JMat = JMat.inverse().eval();

        Eigen::MatrixXd dndx(2, 4);
        dndx.row(0) = dndxi;
        dndx.row(1) = dndeta;
        dndx = JMat.inverse() * dndx;
        // dndx = JMat * dndx;

        Eigen::MatrixXd Gmat = Eigen::MatrixXd::Zero(6, 12);
        Gmat << dndx(0, 0), 0, 0, dndx(0, 1), 0, 0, dndx(0, 2), 0, 0, dndx(0, 3), 0, 0,
            dndx(1, 0), 0, 0, dndx(1, 1), 0, 0, dndx(1, 2), 0, 0, dndx(1, 3), 0, 0,
            0, dndx(0, 0), 0, 0, dndx(0, 1), 0, 0, dndx(0, 2), 0, 0, dndx(0, 3), 0,
            0, dndx(1, 0), 0, 0, dndx(1, 1), 0, 0, dndx(1, 2), 0, 0, dndx(1, 3), 0,
            0, 0, dndx(0, 0), 0, 0, dndx(0, 1), 0, 0, dndx(0, 2), 0, 0, dndx(0, 3),
            0, 0, dndx(1, 0), 0, 0, dndx(1, 1), 0, 0, dndx(1, 2), 0, 0, dndx(1, 3);

        MembraneStressData strs = stress(disp[0], disp[1], disp[2], disp[3], xi, eta);
        Eigen::Matrix2d SigMat;
        SigMat << strs.sigx, strs.sigxy,
            strs.sigxy, strs.sigy;

        // 2x2ブロックを3つの対角位置に配置
        Eigen::MatrixXd SigMat3 = Eigen::MatrixXd::Zero(6, 6);
        SigMat3.block<2, 2>(0, 0) = SigMat;
        SigMat3.block<2, 2>(2, 2) = SigMat;
        SigMat3.block<2, 2>(4, 4) = SigMat;

        Kg += thickness.plane_thick * Gmat.transpose() * SigMat3 * Gmat * detJ;
        // Kg += thickness.plane_thick * intgp * Gmat.transpose() * SigMat3 * Gmat * detJ;
    }

    // Eigen::MatrixXd tr = trans_matrix();
    // return tr.transpose() * Kg * tr;
    return Kg;
}

QuadPlaneElement::QuadPlaneElement(Node *n0, Node *n1, Node *n2, Node *n3, double t, Material mat)
{
    Nodes[0] = n0;
    Nodes[1] = n1;
    Nodes[2] = n2;
    Nodes[3] = n3;

    thickness = Thickness(t);
    Mat = mat;

    // Point p01 = (n0->Location + n1->Location) / 2;
    Point p12 = (n1->Location + n2->Location) / 2;
    Point p23 = (n2->Location + n3->Location) / 2;
    Point p30 = (n3->Location + n0->Location) / 2;
    plane = Plane::CreateFromPoints(p30, p12, p23);
    // plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
}

QuadPlaneElement::QuadPlaneElement(Node *n0, Node *n1, Node *n2, Node *n3, Thickness t, Material mat)
{
    Nodes[0] = n0;
    Nodes[1] = n1;
    Nodes[2] = n2;
    Nodes[3] = n3;

    thickness = t;
    Mat = mat;

    // Point p01 = (n0->Location + n1->Location) / 2;
    Point p12 = (n1->Location + n2->Location) / 2;
    Point p23 = (n2->Location + n3->Location) / 2;
    Point p30 = (n3->Location + n0->Location) / 2;
    plane = Plane::CreateFromPoints(p30, p12, p23);
}

Eigen::MatrixXd QuadPlaneElement::trans_matrix()
{
    Eigen::Matrix3d tr0 = trans_matrix3(plane);

    // ブロックの数を指定
    const int numBlocks = node_num;

    // 大きな行列を作成し、対角ブロック行列として同じ行列を配置
    Eigen::MatrixXd matrix(node_num * node_dof_local, total_dof);
    matrix.setZero(); // 行列を0で初期化

    for (int i = 0; i < numBlocks; ++i)
        matrix.block(i * 2, i * node_dof, 2, node_dof) = tr0.block(0, 0, 2, node_dof);

    return matrix;
}

Eigen::MatrixXd QuadPlaneElement::NodeConsistentMass()
{
    return Eigen::MatrixXd::Identity(total_dof, total_dof);
}

// std::vector<NodeLoadData> QuadPlaneElement::BodyforceToNodeLoadData(Eigen::Vector3d accel_vec)
//{
//	return std::vector<NodeLoadData>();
// }

std::vector<NodeLoadData> QuadPlaneElement::InertialForceToNodeLoadData(Eigen::Vector3d accel_vec)
{
    Eigen::Matrix<double, total_dof, total_dof> m = NodeConsistentMass();
    Eigen::Vector<double, total_dof> f;
    f << accel_vec.x(), accel_vec.y(), accel_vec.z(),
        accel_vec.x(), accel_vec.y(), accel_vec.z(),
        accel_vec.x(), accel_vec.y(), accel_vec.z(),
        accel_vec.x(), accel_vec.y(), accel_vec.z();

    f = m * f;
    std::vector<NodeLoadData> loads;
    loads.push_back(NodeLoadData(Nodes[0]->id, f(0), f(1), f(2)));
    loads.push_back(NodeLoadData(Nodes[1]->id, f(3), f(4), f(5)));
    loads.push_back(NodeLoadData(Nodes[2]->id, f(6), f(7), f(8)));
    loads.push_back(NodeLoadData(Nodes[3]->id, f(9), f(10), f(11)));
    return loads;
}


std::vector<NodeLoadData> QuadPlaneElement::AreaForceToNodeLoadData(std::vector<Vector> load_vecs)
{
    // 2-Point Gauss Quadrature
    // const Eigen::Vector2d intg_weights(1, 1);
    // const Eigen::Vector2d intg_params(-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0));

    // 3-Point Gauss Quadrature
    const Eigen::Vector3d intg_weights(5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0);
    const Eigen::Vector3d intg_params(-sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0));

    // 4-Point Gauss Quadrature
    // const Eigen::Vector4d intg_weights(0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451);
    // const Eigen::Vector4d intg_params(-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116);

    // 5-Point Gauss Quadrature
    // Eigen::VectorXd intg_weights(5), intg_params(5);
    // intg_weights << 0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851;
    // intg_params << -0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459;

    double t3 = thickness.plane_thick * thickness.plane_thick * thickness.plane_thick;

    Eigen::MatrixXd M(total_dof, total_dof);
    M.setZero();

    Eigen::Vector<double, 12> p_vec;
    p_vec << load_vecs[0].x, load_vecs[0].y, load_vecs[0].z,
        load_vecs[1].x, load_vecs[1].y, load_vecs[1].z,
        load_vecs[2].x, load_vecs[2].y, load_vecs[2].z,
        load_vecs[3].x, load_vecs[3].y, load_vecs[3].z;
    // Eigen::Vector<int, 12> shape_indices;
    // const int shape_indices[12] = { 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20 };
    //  shape_indices << 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20;

    // Eigen::Vector<int, 12> HVecs_indices;
    // HVecs_indices << 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22;
    // const int HVecs_indices[12] = { 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22 };

    for (size_t ix = 0; ix < intg_params.size(); ix++)
    {
        for (size_t iy = 0; iy < intg_params.size(); iy++)
        {
            double detJ = JMatrix(intg_params(ix), intg_params(iy)).determinant();

            // ui, vi, wi...
            Eigen::Vector4d spf = ShapeFunction4(intg_params(ix), intg_params(iy));
            Eigen::MatrixXd shape_mat = Eigen::MatrixXd::Zero(3, 12);
            shape_mat << spf(0), 0, 0, spf(1), 0, 0, spf(2), 0, 0, spf(3), 0, 0,
                0, spf(0), 0, 0, spf(1), 0, 0, spf(2), 0, 0, spf(3), 0,
                0, 0, spf(0), 0, 0, spf(1), 0, 0, spf(2), 0, 0, spf(3);

            M += shape_mat.transpose() * shape_mat * intg_weights(ix) * intg_weights(iy) * detJ;
        }
    }

    p_vec = M * p_vec;
    std::vector<NodeLoadData> node_loads;
    node_loads.push_back(NodeLoadData(Nodes[0]->id, p_vec(0), p_vec(1), p_vec(2)));
    node_loads.push_back(NodeLoadData(Nodes[1]->id, p_vec(3), p_vec(4), p_vec(5)));
    node_loads.push_back(NodeLoadData(Nodes[2]->id, p_vec(6), p_vec(7), p_vec(8)));
    node_loads.push_back(NodeLoadData(Nodes[3]->id, p_vec(9), p_vec(10), p_vec(11)));

    return node_loads;
}

double QuadPlaneElement::Area()
{
    Vector v01 = Nodes[1]->Location - Nodes[0]->Location;
    Vector v03 = Nodes[3]->Location - Nodes[0]->Location;
    Vector vn13 = Vector::cross(v01, v03);

    Vector v23 = Nodes[3]->Location - Nodes[2]->Location;
    Vector v21 = Nodes[1]->Location - Nodes[2]->Location;
    Vector vn31 = Vector::cross(v01, v03);
    return (vn13.norm() + vn31.norm()) / 2;
}

Eigen::MatrixXd QuadPlaneElement::StiffnessMatrix()
{
    Eigen::MatrixXd K = localStiffnessMatrix();
    Eigen::MatrixXd TrMat = trans_matrix();
    return TrMat.transpose() * K * TrMat;
}

void QuadPlaneElement::AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K)
{
    int indices[total_dof];
    for (size_t i = 0; i < node_num; i++)
        for (size_t j = 0; j < node_dof; j++)
            indices[node_dof * i + j] = Nodes[i]->id * 6 + j;

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

void QuadPlaneElement::AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat)
{
    AssembleMatrix(mat, StiffnessMatrix());
    // int indices[total_dof];
    // for (size_t i = 0; i < node_num; i++)
    // 	for (size_t j = 0; j < node_dof; j++)
    // 		indices[node_dof * i + j] = Nodes[i]->id * 6 + j;

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

void QuadPlaneElement::AssembleGeometricStiffMatrix(
    Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp)
{
    AssembleMatrix(mat, geometric_local_stiffness_matrix(disp));
}

void QuadPlaneElement::AssembleMassMatrix(Eigen::SparseMatrix<double> &mat)
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

MembraneStressData QuadPlaneElement::stress(Displacement d0, Displacement d1,
                                            Displacement d2, Displacement d3, double xi, double eta)
{
    Eigen::VectorXd wvec(total_dof);
    wvec.segment(0, node_dof) = Eigen::Map<Eigen::VectorXd>(d0.displace, node_dof);
    wvec.segment(3, node_dof) = Eigen::Map<Eigen::VectorXd>(d1.displace, node_dof);
    wvec.segment(6, node_dof) = Eigen::Map<Eigen::VectorXd>(d2.displace, node_dof);
    wvec.segment(9, node_dof) = Eigen::Map<Eigen::VectorXd>(d3.displace, node_dof);

    wvec = trans_matrix() * wvec;
    Eigen::Vector3d strs = DMatrix() * BMatrix(xi, eta) * wvec;

    return MembraneStressData(strs(0), strs(1), strs(2));
    // return MembraneStressData();
}

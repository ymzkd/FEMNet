#include "TriPlateElement.h"

// Eigen::MatrixXd TriPlateElement::BMatrix_Save(double xi, double eta)
//{
//	Point p1 = plane.PointToCoord(Nodes[0]->Location);
//	Point p2 = plane.PointToCoord(Nodes[1]->Location);
//	Point p3 = plane.PointToCoord(Nodes[2]->Location);
//
//	Vector v23 = p2 - p3;
//	Vector v31 = p3 - p1;
//	Vector v12 = p1 - p2;
//
//	double p4 = -6 * v23.x / v23.squared_norm();
//	double p5 = -6 * v31.x / v31.squared_norm();
//	double p6 = -6 * v12.x / v12.squared_norm();
//
//	double q4 = 3 * v23.x * v23.y / v23.squared_norm();
//	double q5 = 3 * v31.x * v31.y / v31.squared_norm();
//	double q6 = 3 * v12.x * v12.y / v12.squared_norm();
//
//	double r4 = 3 * v23.y * v23.y / v23.squared_norm();
//	double r5 = 3 * v31.y * v31.y / v31.squared_norm();
//	double r6 = 3 * v12.y * v12.y / v12.squared_norm();
//
//	double t4 = -6 * v23.y / v23.squared_norm();
//	double t5 = -6 * v31.y / v31.squared_norm();
//	double t6 = -6 * v12.y / v12.squared_norm();
//
//	Eigen::Vector<double, 9> Hx_xi, Hy_xi, Hx_eta, Hy_eta;
//	Hx_xi <<
//		p6 * (1 - 2 * xi) + (p5 - p6) * eta,
//		q6 * (1 - 2 * xi) - (q5 + q6) * eta,
//		-4 + 6 * (xi + eta) + r6 * (1 - 2 * xi) - eta * (r5 + r6),
//		-p6 * (1 - 2 * xi) + eta * (p4 + p6),
//		q6 * (1 - 2 * xi) - eta * (q6 - q4),
//		-2 + 6 * xi + r6 * (1 - 2 * xi) + eta * (r4 - r6),
//		-eta * (p5 + p4),
//		eta * (q4 - q5),
//		-eta * (r5 - r4);
//
//	// 元論文
//	Hy_xi <<
//		t6 * (1 - 2 * xi) + (t5 - t6) * eta,
//		1 + r6 * (1 - 2 * xi) - (r5 + r6) * eta,
//		-q6 * (1 - 2 * xi) + eta * (q5 + q6),
//		-t6 * (1 - 2 * xi) + eta * (t4 + t6),
//		-1 + r6 * (1 - 2 * xi) + eta * (r4 - r6),
//		-q6 * (1 - 2 * xi) - eta * (q4 - q6),
//		-eta * (t4 + t5),
//		eta * (r4 - r5),
//		-eta * (q4 - q5);
//
//	Hx_eta <<
//		-p5 * (1 - 2 * eta) - xi * (p6 - p5),
//		q5 * (1 - 2 * eta) - xi * (q5 + q6),
//		-4 + 6 * (xi + eta) + r5 * (1 - 2 * eta) - xi * (r6 + r5),
//		xi * (p4 + p6),
//		xi * (q4 - q6),
//		-xi * (r6 - r4),
//		p5 * (1 - 2 * eta) - xi * (p4 + p5),
//		q5 * (1 - 2 * eta) + xi * (q4 - q5),
//		-2 + 6 * eta + r5 * (1 - 2 * eta) + xi * (r4 - r5);
//
//	// 元論文
//	Hy_eta <<
//		-t5 * (1 - 2 * eta) - xi * (t6 - t5),
//		1 + r5 * (1 - 2 * eta) - xi * (r5 + r6),
//		-q5 * (1 - 2 * eta) + xi * (q5 + q6),
//		xi * (t4 + t6),
//		xi * (r4 - r6),
//		-xi * (q4 - q6),
//		t5 * (1 - 2 * eta) - xi * (t4 + t5),
//		-1 + r5 * (1 - 2 * eta) + xi * (r4 - r5),
//		-q5 * (1 - 2 * eta) - xi * (q4 - q5);
//
//	Eigen::MatrixXd mat(3, 9);
//	mat.row(0) = v31.y * Hx_xi + v12.y * Hx_eta;
//	mat.row(1) = -v31.x * Hy_xi - v12.x * Hy_eta;
//	mat.row(2) = -v31.x * Hx_xi - v12.x * Hx_eta + v31.y * Hy_xi + v12.y * Hy_eta;
//	mat /= (2.0 * Area());
//	return mat;
// }

// 入力形状関数に応じたx,yそれぞれの方向の列ベクトルによるHベクトル行列
Eigen::MatrixXd TriPlateElement::HVecs(Eigen::Vector<double, 6> shape_funcs)
{
    Point p1 = plane.PointToCoord(Nodes[0]->Location);
    Point p2 = plane.PointToCoord(Nodes[1]->Location);
    Point p3 = plane.PointToCoord(Nodes[2]->Location);

    Vector v23 = p2 - p3;
    Vector v31 = p3 - p1;
    Vector v12 = p1 - p2;

    double l12 = v12.norm();
    double l23 = v23.norm();
    double l31 = v31.norm();

    double c4 = -v23.y / l23;
    double c5 = -v31.y / l31;
    double c6 = -v12.y / l12;

    double s4 = v23.x / l23;
    double s5 = v31.x / l31;
    double s6 = v12.x / l12;

    // -------------------------------------------------------
    // double L1 = 1 - L2 - L3;
    // double area = Area();

    // Eigen::VectorXd dLdx(3), dLdy(3);
    // dLdx << v23.y * 0.5 / area, v31.y * 0.5 / area, v12.y * 0.5 / area;
    // dLdy << -v23.x * 0.5 / area, -v31.x * 0.5 / area, -v12.x * 0.5 / area;

    // Eigen::VectorXd dndx(6), dndy(6);
    // dndx << dLdx(0) * (4 * L1 - 1), dLdx(1)* (4 * L2 - 1), dLdx(2)* (4 * L3 - 1),
    //	4 * (dLdx(1) * L3 + L2 * dLdx(2)), 4 * (dLdx(2) * L1 + L3 * dLdx(0)),
    //	4 * (dLdx(0) * L2 + L1 * dLdx(1));
    // dndy << dLdy(0) * (4 * L1 - 1), dLdy(1)* (4 * L2 - 1), dLdy(2)* (4 * L3 - 1),
    //	4 * (dLdy(1) * L3 + L2 * dLdy(2)), 4 * (dLdy(2) * L1 + L3 * dLdy(0)),
    //	4 * (dLdy(0) * L2 + L1 * dLdy(1));
    //  -------------------------------------------------------

    Eigen::Matrix<double, 2, 9> HVecs;

    // Hx
    HVecs.row(0) << 1.5 * s5 / l31 * shape_funcs(4) - 1.5 * s6 / l12 * shape_funcs(5),
        -3 * s5 * c5 / 4 * shape_funcs(4) - 3 * s6 * c6 / 4 * shape_funcs(5),
        shape_funcs(0) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * shape_funcs(4) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * shape_funcs(5),
        1.5 * s6 / l12 * shape_funcs(5) - 1.5 * s4 / l23 * shape_funcs(3),
        -3 * s4 * c4 / 4 * shape_funcs(3) - 3 * s6 * c6 / 4 * shape_funcs(5),
        shape_funcs(1) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * shape_funcs(3) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * shape_funcs(5),
        1.5 * s4 / l23 * shape_funcs(3) - 1.5 * s5 / l31 * shape_funcs(4),
        -3 * s4 * c4 / 4 * shape_funcs(3) - 3 * s5 * c5 / 4 * shape_funcs(4),
        shape_funcs(2) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * shape_funcs(3) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * shape_funcs(4);

    // Hy
    HVecs.row(1) << 1.5 * c6 / l12 * shape_funcs(5) - 1.5 * c5 / l31 * shape_funcs(4),
        -shape_funcs(0) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * shape_funcs(4) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * shape_funcs(5),
        3 * s5 * c5 / 4 * shape_funcs(4) + 3 * s6 * c6 / 4 * shape_funcs(5),
        1.5 * c4 / l23 * shape_funcs(3) - 1.5 * c6 / l12 * shape_funcs(5),
        -shape_funcs(1) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * shape_funcs(3) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * shape_funcs(5),
        3 * s4 * c4 / 4 * shape_funcs(3) + 3 * s6 * c6 / 4 * shape_funcs(5),
        1.5 * c5 / l31 * shape_funcs(4) - 1.5 * c4 / l23 * shape_funcs(3),
        -shape_funcs(2) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * shape_funcs(3) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * shape_funcs(4),
        3 * s4 * c4 / 4 * shape_funcs(3) + 3 * s5 * c5 / 4 * shape_funcs(4);

    return HVecs;
}

Eigen::MatrixXd TriPlateElement::BMatrix(double L2, double L3)
{
    Point p1 = plane.PointToCoord(Nodes[0]->Location);
    Point p2 = plane.PointToCoord(Nodes[1]->Location);
    Point p3 = plane.PointToCoord(Nodes[2]->Location);

    Vector v23 = p2 - p3;
    Vector v31 = p3 - p1;
    Vector v12 = p1 - p2;

    double L1 = 1 - L2 - L3;
    double area = Area();

    Eigen::VectorXd dLdx(3), dLdy(3);
    dLdx << v23.y * 0.5 / area, v31.y * 0.5 / area, v12.y * 0.5 / area;
    dLdy << -v23.x * 0.5 / area, -v31.x * 0.5 / area, -v12.x * 0.5 / area;

    Eigen::VectorXd dndx(6), dndy(6);
    dndx << dLdx(0) * (4 * L1 - 1), dLdx(1) * (4 * L2 - 1), dLdx(2) * (4 * L3 - 1),
        4 * (dLdx(1) * L3 + L2 * dLdx(2)), 4 * (dLdx(2) * L1 + L3 * dLdx(0)),
        4 * (dLdx(0) * L2 + L1 * dLdx(1));
    dndy << dLdy(0) * (4 * L1 - 1), dLdy(1) * (4 * L2 - 1), dLdy(2) * (4 * L3 - 1),
        4 * (dLdy(1) * L3 + L2 * dLdy(2)), 4 * (dLdy(2) * L1 + L3 * dLdy(0)),
        4 * (dLdy(0) * L2 + L1 * dLdy(1));

    double l12 = v12.norm();
    double l23 = v23.norm();
    double l31 = v31.norm();

    double c4 = -v23.y / l23;
    double c5 = -v31.y / l31;
    double c6 = -v12.y / l12;

    double s4 = v23.x / l23;
    double s5 = v31.x / l31;
    double s6 = v12.x / l12;

    Eigen::Vector<double, 9> dHx_dx, dHx_dy, dHy_dx, dHy_dy;
    dHx_dx << 1.5 * s5 / l31 * dndx(4) - 1.5 * s6 / l12 * dndx(5),
        -3 * s5 * c5 / 4 * dndx(4) - 3 * s6 * c6 / 4 * dndx(5),
        dndx(0) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * dndx(4) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * dndx(5),
        1.5 * s6 / l12 * dndx(5) - 1.5 * s4 / l23 * dndx(3),
        -3 * s4 * c4 / 4 * dndx(3) - 3 * s6 * c6 / 4 * dndx(5),
        dndx(1) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * dndx(3) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * dndx(5),
        1.5 * s4 / l23 * dndx(3) - 1.5 * s5 / l31 * dndx(4),
        -3 * s4 * c4 / 4 * dndx(3) - 3 * s5 * c5 / 4 * dndx(4),
        dndx(2) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * dndx(3) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * dndx(4);

    dHx_dy << 1.5 * s5 / l31 * dndy(4) - 1.5 * s6 / l12 * dndy(5),
        -3 * s5 * c5 / 4 * dndy(4) - 3 * s6 * c6 / 4 * dndy(5),
        dndy(0) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * dndy(4) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * dndy(5),
        1.5 * s6 / l12 * dndy(5) - 1.5 * s4 / l23 * dndy(3),
        -3 * s4 * c4 / 4 * dndy(3) - 3 * s6 * c6 / 4 * dndy(5),
        dndy(1) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * dndy(3) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * dndy(5),
        1.5 * s4 / l23 * dndy(3) - 1.5 * s5 / l31 * dndy(4),
        -3 * s4 * c4 / 4 * dndy(3) - 3 * s5 * c5 / 4 * dndy(4),
        dndy(2) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * dndy(3) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * dndy(4);

    dHy_dx << 1.5 * c6 / l12 * dndx(5) - 1.5 * c5 / l31 * dndx(4),
        -dndx(0) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * dndx(4) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * dndx(5),
        3 * s5 * c5 / 4 * dndx(4) + 3 * s6 * c6 / 4 * dndx(5),
        1.5 * c4 / l23 * dndx(3) - 1.5 * c6 / l12 * dndx(5),
        -dndx(1) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * dndx(3) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * dndx(5),
        3 * s4 * c4 / 4 * dndx(3) + 3 * s6 * c6 / 4 * dndx(5),
        1.5 * c5 / l31 * dndx(4) - 1.5 * c4 / l23 * dndx(3),
        -dndx(2) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * dndx(3) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * dndx(4),
        3 * s4 * c4 / 4 * dndx(3) + 3 * s5 * c5 / 4 * dndx(4);

    dHy_dy << 1.5 * c6 / l12 * dndy(5) - 1.5 * c5 / l31 * dndy(4),
        -dndy(0) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * dndy(4) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * dndy(5),
        3 * s5 * c5 / 4 * dndy(4) + 3 * s6 * c6 / 4 * dndy(5),
        1.5 * c4 / l23 * dndy(3) - 1.5 * c6 / l12 * dndy(5),
        -dndy(1) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * dndy(3) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * dndy(5),
        3 * s4 * c4 / 4 * dndy(3) + 3 * s6 * c6 / 4 * dndy(5),
        1.5 * c5 / l31 * dndy(4) - 1.5 * c4 / l23 * dndy(3),
        -dndy(2) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * dndy(3) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * dndy(4),
        3 * s4 * c4 / 4 * dndy(3) + 3 * s5 * c5 / 4 * dndy(4);

    Eigen::MatrixXd mat(3, 9);
    mat.row(0) = dHx_dx;
    mat.row(1) = dHy_dy;
    mat.row(2) = dHx_dy + dHy_dx;
    // mat /= (2.0 * Area());
    return mat;
}

// void TriPlateElement::shearstress(Displacement d0, Displacement d1, Displacement d2)
//{
//	Point p1 = plane.PointToCoord(Nodes[0]->Location);
//	Point p2 = plane.PointToCoord(Nodes[1]->Location);
//	Point p3 = plane.PointToCoord(Nodes[2]->Location);
//
//	Vector v23 = p2 - p3;
//	Vector v31 = p3 - p1;
//	Vector v12 = p1 - p2;
//
//	double p4 = -6 * v23.x / v23.squared_norm();
//	double p5 = -6 * v31.x / v31.squared_norm();
//	double p6 = -6 * v12.x / v12.squared_norm();
//
//	double q4 = 3 * v23.x * v23.y / v23.squared_norm();
//	double q5 = 3 * v31.x * v31.y / v31.squared_norm();
//	double q6 = 3 * v12.x * v12.y / v12.squared_norm();
//
//	double r4 = 3 * v23.y * v23.y / v23.squared_norm();
//	double r5 = 3 * v31.y * v31.y / v31.squared_norm();
//	double r6 = 3 * v12.y * v12.y / v12.squared_norm();
//
//	double t4 = -6 * v23.y / v23.squared_norm();
//	double t5 = -6 * v31.y / v31.squared_norm();
//	double t6 = -6 * v12.y / v12.squared_norm();
//
//	double Jl = v12.y * v31.x - v12.x * v31.y;
//	Eigen::Matrix2d Jinv;
//	Jinv << v31.y / Jl, v12.y / Jl, -v31.x / Jl, -v12.x / Jl;
//
//	Eigen::Vector<double, 9> Hy_xi_xi, Hy_xi_eta, Hy_eta_xi, Hy_eta_eta;
//	Hy_xi_xi << -2.0 * t6, -2.0 * r6, 2.0 * q6, 2.0 * t6, -2.0 * r6, 2.0 * q6, 0, 0, 0;
//	Hy_xi_eta << t5 - t6, -(r5 + r6), q5 + q6, t4 + t6, r4 - r6, -(q4 - q6), -(t4 + t5), r4 - r5, q5 - q4;
//	Hy_eta_xi << t5 - t6, -(r5 + r6), q5 + q6, t4 + t6, r4 - r6, -(q4 - q6), -(t4 + t5), r4 - r5, q5 - q4;
//	Hy_eta_eta << 2.0 * t5, -2.0 * r5, 2.0 * q5, 0, 0, 0, -2.0 * t5, -2.0 * r5, 2.0 * q5;
//
//	Eigen::Vector<double, 9> Hx_xi_xi, Hx_xi_eta, Hx_eta_xi, Hx_eta_eta;
//	Hx_xi_xi << -2.0 * p6, -2.0 * q6, 6.0 - 2.0 * r6, 2.0 * p6, -2.0 * q6, 6.0 - 2.0 * r6, 0, 0, 0;
//	Hx_xi_eta << p5 - p6, -(q5 + q6), 6.0 - (r5 + r6), p4 + p6, q4 - q6, r4 - r6, -(p5 + p4), (q4 - q5), -(r5 - r4);
//	Hx_eta_xi << p5 - p6, -(q5 + q6), 6.0 - (r5 + r6), p4 + p6, q4 - q6, r4 - r6, -(p4 + p5), q4 - q5, r4 - r5;
//	Hx_eta_eta << 2.0 * p5, -2.0 * q5, 6.0 - 2.0 * r5, 0, 0, 0, -2.0 * p5, -2.0 * q5, 6.0 - 2.0 * r5;
//
//	Eigen::Vector<double, 9> Hy_xi_x, Hy_xi_y, Hy_eta_x, Hy_eta_y;
//	Hy_xi_x = Jinv(0, 0) * Hy_xi_xi + Jinv(0, 1) * Hy_xi_eta;
//	Hy_xi_y = Jinv(1, 0) * Hy_xi_xi + Jinv(1, 1) * Hy_xi_eta;
//	Hy_eta_x = Jinv(0, 0) * Hy_eta_xi + Jinv(0, 1) * Hy_eta_eta;
//	Hy_eta_y = Jinv(1, 0) * Hy_eta_xi + Jinv(1, 1) * Hy_eta_eta;
//
//	Eigen::Vector<double, 9> Hx_xi_x, Hx_xi_y, Hx_eta_x, Hx_eta_y;
//	Hx_xi_x = Jinv(0, 0) * Hx_xi_xi + Jinv(0, 1) * Hx_xi_eta;
//	Hx_xi_y = Jinv(1, 0) * Hx_xi_xi + Jinv(1, 1) * Hx_xi_eta;
//	Hx_eta_x = Jinv(0, 0) * Hx_eta_xi + Jinv(0, 1) * Hx_eta_eta;
//	Hx_eta_y = Jinv(1, 0) * Hx_eta_xi + Jinv(1, 1) * Hx_eta_eta;
//
//	Eigen::MatrixXd Bmat_y(3, 9);
//	Bmat_y.row(0) = v31.y * Hx_xi_y + v12.y * Hx_eta_y;
//	Bmat_y.row(1) = -v31.x * Hy_xi_y - v12.x * Hy_eta_y;
//	Bmat_y.row(2) = -v31.x * Hx_xi_y - v12.x * Hx_eta_y + v31.y * Hy_xi_y + v12.y * Hy_eta_y;
//	Bmat_y /= (2.0 * Area());
//
//	Eigen::MatrixXd Bmat_x(3, 9);
//	Bmat_x.row(0) = v31.y * Hx_xi_x + v12.y * Hx_eta_x;
//	Bmat_x.row(1) = -v31.x * Hy_xi_x - v12.x * Hy_eta_x;
//	Bmat_x.row(2) = -v31.x * Hx_xi_x - v12.x * Hx_eta_x + v31.y * Hy_xi_x + v12.y * Hy_eta_x;
//	Bmat_x /= (2.0 * Area());
//
//	Eigen::Matrix3d tr = trans_matrix3(plane);
//	Displacement d0t = d0.translate(tr);
//	Displacement d1t = d1.translate(tr);
//	Displacement d2t = d2.translate(tr);
//
//	Eigen::VectorXd wvec(9);
//	wvec << d0t.Dz(), d0t.Rx(), d0t.Ry(), d1t.Dz(), d1t.Rx(), d1t.Ry(), d2t.Dz(), d2t.Rx(), d2t.Ry();
//	Eigen::Matrix3d Dmat = DMatrix(); // .topLeftCorner<2, 2>();
//
//	Eigen::Vector3d dMdx = Dmat * (Bmat_x * wvec);
//	Eigen::Vector3d dMdy = Dmat * (Bmat_y * wvec);
//
//	double qx = dMdx(0) + dMdy(2);
//	double qy = dMdy(1) + dMdx(2);
//	// double qy = Dmat.row(1) * (Bmat_y * wvec);
//
//	std::cout << "qx: " << qx << std::endl;
//	std::cout << "qy: " << qy << std::endl;
// }

Eigen::Matrix3d TriPlateElement::DMatrix()
{
    Eigen::Matrix3d mat;
    mat << 1, Mat.Poisson, 0,
        Mat.Poisson, 1, 0,
        0, 0, (1 - Mat.Poisson) / 2;
    return Mat.Young * pow(thickness.plate_thick, 3) / 12 / (1 - Mat.Poisson * Mat.Poisson) * mat;
}

// Eigen::MatrixXd TriPlateElement::localStiffnessMatrix_Save()
//{
//	double weight = 1.0 / 3.0;
//	double xi_list[3]{ 0.5 , 0.5 , 0 };
//	double eta_list[3]{ 0, 0.5, 0.5 };
//	double a2 = Area();
//	//double a2 = 2.0 * Area();
//
//	Eigen::Matrix3d DMat = DMatrix();
//	Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(9, 9);
//	for (size_t i = 0; i < 3; i++)
//	{
//		Eigen::MatrixXd BMat = BMatrix_Save(xi_list[i], eta_list[i]);
//		mat += BMat.transpose() * DMat * BMat;
//	}
//	return a2 * weight * mat;
// }

Eigen::MatrixXd TriPlateElement::localStiffnessMatrix()
{
    double weight = 1.0 / 3.0;
    double xi_list[3]{0.5, 0.5, 0};
    double eta_list[3]{0, 0.5, 0.5};

    Eigen::Matrix3d DMat = DMatrix();
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(9, 9);
    for (size_t i = 0; i < 3; i++)
    {
        Eigen::MatrixXd BMat = BMatrix(xi_list[i], eta_list[i]);
        mat += BMat.transpose() * DMat * BMat * weight;
    }
    return mat * Area();
}

Eigen::MatrixXd TriPlateElement::geometric_local_stiffness_matrix(const std::vector<Displacement> &disp)
{
    Eigen::MatrixXd Kg = Eigen::MatrixXd::Zero(9, 9);

    double weight = 1.0 / 3.0;
    double xi_list[3]{0.5, 0.5, 0};
    double eta_list[3]{0, 0.5, 0.5};

    Point p1 = plane.PointToCoord(Nodes[0]->Location);
    Point p2 = plane.PointToCoord(Nodes[1]->Location);
    Point p3 = plane.PointToCoord(Nodes[2]->Location);

    Vector v23 = p2 - p3;
    Vector v31 = p3 - p1;
    Vector v12 = p1 - p2;

    double l12 = v12.norm();
    double l23 = v23.norm();
    double l31 = v31.norm();

    double c4 = -v23.y / l23;
    double c5 = -v31.y / l31;
    double c6 = -v12.y / l12;

    double s4 = v23.x / l23;
    double s5 = v31.x / l31;
    double s6 = v12.x / l12;
    double area = Area();

    for (size_t i = 0; i < 3; i++)
    {
        double L2 = xi_list[i];
        double L3 = eta_list[i];
        double L1 = 1 - L2 - L3;

        Eigen::VectorXd shape_func(6);
        shape_func << 2 * L1 * L1 - L1, 2 * L2 * L2 - L2, 2 * L3 * L3 - L3,
            4 * L2 * L3, 4 * L3 * L1, 4 * L1 * L2;

        // Hx, Hyを計算
        Eigen::Vector<double, 9> Hx, Hy;
        Hx << 1.5 * s5 / l31 * shape_func(4) - 1.5 * s6 / l12 * shape_func(5),
            -3 * s5 * c5 / 4 * shape_func(4) - 3 * s6 * c6 / 4 * shape_func(5),
            shape_func(0) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * shape_func(4) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * shape_func(5),
            1.5 * s6 / l12 * shape_func(5) - 1.5 * s4 / l23 * shape_func(3),
            -3 * s4 * c4 / 4 * shape_func(3) - 3 * s6 * c6 / 4 * shape_func(5),
            shape_func(1) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * shape_func(3) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * shape_func(5),
            1.5 * s4 / l23 * shape_func(3) - 1.5 * s5 / l31 * shape_func(4),
            -3 * s4 * c4 / 4 * shape_func(3) - 3 * s5 * c5 / 4 * shape_func(4),
            shape_func(2) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * shape_func(3) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * shape_func(4);

        Hy << 1.5 * c6 / l12 * shape_func(5) - 1.5 * c5 / l31 * shape_func(4),
            -shape_func(0) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * shape_func(4) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * shape_func(5),
            3 * s5 * c5 / 4 * shape_func(4) + 3 * s6 * c6 / 4 * shape_func(5),
            1.5 * c4 / l23 * shape_func(3) - 1.5 * c6 / l12 * shape_func(5),
            -shape_func(1) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * shape_func(3) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * shape_func(5),
            3 * s4 * c4 / 4 * shape_func(3) + 3 * s6 * c6 / 4 * shape_func(5),
            1.5 * c5 / l31 * shape_func(4) - 1.5 * c4 / l23 * shape_func(3),
            -shape_func(2) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * shape_func(3) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * shape_func(4),
            3 * s4 * c4 / 4 * shape_func(3) + 3 * s5 * c5 / 4 * shape_func(4);

        PlateStressData strs = stress(disp[0], disp[1], disp[2], L2, L3);
        Eigen::Matrix2d SigMat;
        SigMat << strs.My, strs.Mxy, strs.Mxy, strs.Mx;

        Eigen::MatrixXd Gmat = Eigen::MatrixXd::Zero(2, 9);
        Gmat.row(0) = Hx;
        Gmat.row(1) = Hy;

        Kg += Gmat.transpose() * SigMat * Gmat * weight * area;
    }

    // Eigen::MatrixXd trMat = trans_matrix();
    // return trMat.transpose() * Kg * trMat;
    return Kg;
}

TriPlateElement::TriPlateElement(Node *n0, Node *n1, Node *n2, Thickness t, Material mat, double beta)
    : plane_element(n0, n1, n2, t, mat)
{
    Nodes[0] = n0;
    Nodes[1] = n1;
    Nodes[2] = n2;

    Mat = mat;
    thickness = t;
    plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
    plane.Rotate(beta, plane.ez);
}

TriPlateElement::TriPlateElement(Node *n0, Node *n1, Node *n2, double t, Material mat, double beta)
    : plane_element(n0, n1, n2, t, mat)
{
    Nodes[0] = n0;
    Nodes[1] = n1;
    Nodes[2] = n2;

    Mat = mat;
    thickness = Thickness(t);
    plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
    plane.Rotate(beta, plane.ez);
}

/// <summary>
/// 18 x 18 blocked translate matrix.
/// </summary>
TriPlateElement::LocalMatrixd
TriPlateElement::trans_matrix()
{
    Eigen::Matrix3d tr0 = trans_matrix3(plane);

    // ブロックの数を指定
    const int numBlocks = 6;

    // 大きな行列を作成し、対角ブロック行列として同じ行列を配置
    LocalMatrixd matrix;
    matrix.setZero();

    for (int i = 0; i < numBlocks; ++i)
        matrix.block(i * tr0.rows(), i * tr0.cols(), tr0.rows(), tr0.cols()) = tr0;

    return matrix;
}

// Eigen::MatrixXd TriPlateElement::NodeConsistentMass()
//{
//	return Eigen::MatrixXd::Identity(total_dof, total_dof);
// }

std::vector<NodeLoadData> TriPlateElement::InertialForceToNodeLoadData(Eigen::Vector3d accel_vec)
{
    Eigen::Matrix<double, total_dof, total_dof> m;
    m = NodeConsistentMass();
    Eigen::Vector<double, total_dof> f;
    f.setZero();

    f << accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0,
        accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0,
        accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0;

    f = m * f;
    std::vector<NodeLoadData> loads;
    loads.push_back(NodeLoadData(Nodes[0]->id, f(0), f(1), f(2), f(3), f(4), f(5)));
    loads.push_back(NodeLoadData(Nodes[1]->id, f(6), f(7), f(8), f(9), f(10), f(11)));
    loads.push_back(NodeLoadData(Nodes[2]->id, f(12), f(13), f(14), f(15), f(16), f(17)));
    return loads;
}

std::vector<NodeLoadData> TriPlateElement::AreaForceToNodeLoadData(std::vector<Vector> load_vecs)
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

Eigen::MatrixXd TriPlateElement::NodeConsistentMass()
{
    // 1-Point Gauss Quadrature
    // Eigen::Vector<double, 1> intg_weights; intg_weights(0) = 1.0;
    // Eigen::Vector<double, 1> xi_params; xi_params(0) = 1.0 / 3.0;
    // Eigen::Vector<double, 1> eta_params; eta_params(0) = 1.0 / 3.0;

    // 3-Point Gauss Quadrature
    // const Eigen::Vector3d intg_weights(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    // const Eigen::Vector3d xi_params(0.5, 0.5, 0);
    // const Eigen::Vector3d eta_params(0, 0.5, 0.5);

    // 4-Point Gauss Quadrature
    // const Eigen::Vector4d intg_weights(-27.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0);
    // const Eigen::Vector4d xi_params(1.0 / 3.0, 0.6, 0.2, 0.2);
    // const Eigen::Vector4d eta_params(1.0 / 3.0, 0.2, 0.6, 0.2);

    // 7-Point Gauss Quadrature(この精度が必要になりそう)
    Eigen::Vector<double, 7> intg_weights;
    Eigen::Vector<double, 7> xi_params;
    Eigen::Vector<double, 7> eta_params;
    intg_weights << 0.225, 0.1323941527, 0.1323941527, 0.1323941527, 0.1259391805, 0.1259391805, 0.1259391805;
    xi_params << 1.0 / 3.0, 0.0597158717, 0.4701420641, 0.4701420641, 0.7974269853, 0.1012865073, 0.1012865073;
    eta_params << 1.0 / 3.0, 0.4701420641, 0.0597158717, 0.4701420641, 0.1012865073, 0.7974269853, 0.1012865073;

    //// 9-Point Gauss Quadrature
    // Eigen::VectorXd intg_weights(9);
    // intg_weights << 0.2059505047608870, 0.2059505047608870, 0.2059505047608870,
    //	0.0636914142862230, 0.0636914142862230, 0.0636914142862230,
    //	0.0636914142862230, 0.0636914142862230, 0.0636914142862230;
    // Eigen::VectorXd xi_params(9);
    // xi_params << 0.1249495032332320, 0.4375252483838400, 0.4375252483838400,
    //	0.7971126518600710, 0.1654099273984100, 0.0374774207500880,
    //	0.7971126518600710, 0.0374774207500880, 0.1654099273984100;
    // Eigen::VectorXd eta_params(9);
    // eta_params << 0.4375252483838400, 0.1249495032332320, 0.4375252483838400,
    //	0.1654099273984100, 0.0374774207500880, 0.7971126518600710,
    //	0.0374774207500880, 0.1654099273984100, 0.7971126518600710;

    double t3 = thickness.plane_thick * thickness.plane_thick * thickness.plane_thick;
    double area = Area();

    Eigen::MatrixXd M(total_dof, total_dof);
    M.setZero();

    // inplane mass matrix
    Eigen::Matrix<double, 9, 9> Mp;
    Mp << 1.0 / 6, 0, 0, 1.0 / 12, 0, 0, 1.0 / 12, 0, 0,
        0, 1.0 / 6, 0, 0, 1.0 / 12, 0, 0, 1.0 / 12, 0,
        0, 0, 1.0 / 6, 0, 0, 1.0 / 12, 0, 0, 1.0 / 12,
        1.0 / 12, 0, 0, 1.0 / 6, 0, 0, 1.0 / 12, 0, 0,
        0, 1.0 / 12, 0, 0, 1.0 / 6, 0, 0, 1.0 / 12, 0,
        0, 0, 1.0 / 12, 0, 0, 1.0 / 6, 0, 0, 1.0 / 12,
        1.0 / 12, 0, 0, 1.0 / 12, 0, 0, 1.0 / 6, 0, 0,
        0, 1.0 / 12, 0, 0, 1.0 / 12, 0, 0, 1.0 / 6, 0,
        0, 0, 1.0 / 12, 0, 0, 1.0 / 12, 0, 0, 1.0 / 6;

    Mp *= Mat.dense * thickness.plane_thick * area;
    // std::cout << "Mp: \n" << Mp << std::endl;
    // std::cout << "area: " << area << ", " << Mat.dense << ", " << thickness.plane_thick << std::endl;
    // Eigen::Vector<int, 9> shape_indices;
    // shape_indices << 0, 1, 2, 6, 7, 8, 12, 13, 14;
    const int shape_indices[9]{0, 1, 2, 6, 7, 8, 12, 13, 14};

    for (size_t j = 0; j < 9; j++)
        for (size_t k = 0; k < 9; k++)
            M(shape_indices[j], shape_indices[k]) += Mp(j, k);

    // Eigen::Vector<int, 9> HVecs_indices;
    // HVecs_indices << 2, 3, 4, 8, 9, 10, 14, 15, 16;
    const int HVecs_indices[9]{2, 3, 4, 8, 9, 10, 14, 15, 16};

    for (size_t ix = 0; ix < xi_params.size(); ix++)
    {
        // wi, theta_xi, theta_yi...
        Eigen::MatrixXd H_mat = HVecs(ShapeFunctionTriangle6(xi_params(ix), eta_params(ix)));
        // std::cout << "sum of shape func: " << ShapeFunctionSerendipity8(intg_params(ix), intg_params(iy)).sum() << std::endl;
        Eigen::Matrix<double, 9, 9> MassComp;
        MassComp = H_mat.transpose() * H_mat;
        MassComp *= area * Mat.dense * t3 / 12.0 * intg_weights(ix);
        for (size_t j = 0; j < 9; j++)
            for (size_t k = 0; k < 9; k++)
                M(HVecs_indices[j], HVecs_indices[k]) += MassComp(j, k);
    }

    Eigen::MatrixXd trMat = trans_matrix();
    return trMat.transpose() * M * trMat;
}

double TriPlateElement::Area()
{
    Vector v01 = Nodes[1]->Location - Nodes[0]->Location;
    Vector v02 = Nodes[2]->Location - Nodes[0]->Location;
    Vector vn = Vector::cross(v01, v02);
    return vn.norm() / 2;
}

Eigen::MatrixXd TriPlateElement::StiffnessMatrix()
{
    Eigen::MatrixXd Kpln = plane_element.localStiffnessMatrix();
    Eigen::MatrixXd Kplt = localStiffnessMatrix();
    // Eigen::MatrixXd Kplt = localStiffnessMatrix();
    // std::cout << "K plane: \n" << Kpln << std::endl;
    // std::cout << "K plate: \n" << Kplt << std::endl;
    double max_diag = std::max(Kpln.diagonal().maxCoeff(), Kplt.diagonal().maxCoeff());

    Eigen::Matrix3d Krotz = Eigen::Matrix3d::Identity() * max_diag / 1000;

    // Krotz << 1, -0.5, -0.5,
    //	-0.5, 1, -0.5,
    //	-0.5, -0.5, 1;
    // Krotz *=  0.03 * Mat->Young * thickness * Area() / 1000;

    // Krotz.fill(100);
    // std::cout << "K rotz: \n" << Krotz << std::endl;
    // std::cout << "-------------------" << std::endl;
    int indices_rotz[3]{5, 11, 17};
    int indices_pln[6]{0, 1, 6, 7, 12, 13};
    int indices_plt[9]{2, 3, 4, 8, 9, 10, 14, 15, 16};

    Eigen::MatrixXd mat(total_dof, total_dof);
    mat.setZero();
    for (size_t i = 0; i < 6; i++)
    {
        int ir = indices_pln[i];
        for (size_t j = 0; j < 6; j++)
        {
            int ic = indices_pln[j];
            mat(ir, ic) = Kpln(i, j);
        }
    }

    for (size_t i = 0; i < 9; i++)
    {
        int ir = indices_plt[i];
        for (size_t j = 0; j < 9; j++)
        {
            int ic = indices_plt[j];
            mat(ir, ic) = Kplt(i, j);
        }
    }

    for (size_t i = 0; i < 3; i++)
    {
        int ir = indices_rotz[i];
        for (size_t j = 0; j < 3; j++)
        {
            int ic = indices_rotz[j];
            mat(ir, ic) = Krotz(i, j);
        }
    }

    // std::cout << "K: \n" << mat << std::endl;
    // std::cout << "-------------------" << std::endl;
    // std::cout << "rows: " << mat.rows() << "cols: " << mat.cols() << std::endl;
    Eigen::MatrixXd trMat = trans_matrix();
    return trMat.transpose() * mat * trMat;
}

Eigen::MatrixXd TriPlateElement::GeometricStiffnessMatrix(const std::vector<Displacement> &disp)
{
    Eigen::MatrixXd KGpln = plane_element.geometric_local_stiffness_matrix(disp);
    Eigen::MatrixXd KGplt = geometric_local_stiffness_matrix(disp);

    // int indices_rotz[3]{ 5,11,17 };
    int indices_pln[9]{0, 1, 2, 6, 7, 8, 12, 13, 14};
    // int indices_pln[9]{ 0,1,-1,6,7,-1,12,13,-1 };
    int indices_plt[9]{2, 3, 4, 8, 9, 10, 14, 15, 16};

    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(total_dof, total_dof);

    for (size_t i = 0; i < 9; i++)
    {
        int ir = indices_pln[i];
        if (ir < 0)
            continue; // -1は無視
        for (size_t j = 0; j < 9; j++)
        {
            int ic = indices_pln[j];
            if (ic < 0)
                continue; // -1は無視
            mat(ir, ic) += KGpln(i, j);
        }
    }

    for (size_t i = 0; i < 9; i++)
    {
        int ir = indices_plt[i];
        for (size_t j = 0; j < 9; j++)
        {
            int ic = indices_plt[j];
            mat(ir, ic) += KGplt(i, j);
        }
    }

    Eigen::MatrixXd trMat = trans_matrix();
    return trMat.transpose() * mat * trMat;
}

void TriPlateElement::AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K)
{
    int indices[total_dof];
    for (size_t i = 0; i < node_dof; i++)
        indices[i] = Nodes[0]->id * 6 + i;
    for (size_t i = 0; i < node_dof; i++)
        indices[i + node_dof] = Nodes[1]->id * 6 + i;
    for (size_t i = 0; i < node_dof; i++)
        indices[i + node_dof * 2] = Nodes[2]->id * 6 + i;

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

void TriPlateElement::AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat)
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

void TriPlateElement::AssembleGeometricStiffMatrix(
    Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp)
{
    AssembleMatrix(mat, GeometricStiffnessMatrix(disp));
}

void TriPlateElement::AssembleMassMatrix(Eigen::SparseMatrix<double> &mat)
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

/// @brief DKT要素向けのHxベクトルを計算する。
/// @param nfuncs 形状関数やその微分
/// @return 入力形状関数の階数に応じたHxベクトル
Eigen::Vector<double, 9> TriPlateElement_Hx_vector(
    double s4, double s5, double s6,
    double c4, double c5, double c6,
    double l12, double l23, double l31,
    const Eigen::VectorXd &nfuncs)
{
    Eigen::Vector<double, 9> Hx_vector;
    Hx_vector << 1.5 * s5 / l31 * nfuncs(4) - 1.5 * s6 / l12 * nfuncs(5),
        -3 * s5 * c5 / 4 * nfuncs(4) - 3 * s6 * c6 / 4 * nfuncs(5),
        nfuncs(0) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * nfuncs(4) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * nfuncs(5),
        1.5 * s6 / l12 * nfuncs(5) - 1.5 * s4 / l23 * nfuncs(3),
        -3 * s4 * c4 / 4 * nfuncs(3) - 3 * s6 * c6 / 4 * nfuncs(5),
        nfuncs(1) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * nfuncs(3) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * nfuncs(5),
        1.5 * s4 / l23 * nfuncs(3) - 1.5 * s5 / l31 * nfuncs(4),
        -3 * s4 * c4 / 4 * nfuncs(3) - 3 * s5 * c5 / 4 * nfuncs(4),
        nfuncs(2) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * nfuncs(3) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * nfuncs(4);
    return Hx_vector;
}

/// @brief DKT要素向けのHyベクトルを計算する。
/// @param nfuncs 形状関数やその微分
/// @return 入力形状関数の階数に応じたHyベクトル
Eigen::Vector<double, 9> TriPlateElement_Hy_vector(
    double s4, double s5, double s6,
    double c4, double c5, double c6,
    double l12, double l23, double l31,
    const Eigen::VectorXd &nfuncs)
{
    Eigen::Vector<double, 9> Hy_vector;
    Hy_vector << 1.5 * c6 / l12 * nfuncs(5) - 1.5 * c5 / l31 * nfuncs(4),
        -nfuncs(0) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * nfuncs(4) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * nfuncs(5),
        3 * s5 * c5 / 4 * nfuncs(4) + 3 * s6 * c6 / 4 * nfuncs(5),
        1.5 * c4 / l23 * nfuncs(3) - 1.5 * c6 / l12 * nfuncs(5),
        -nfuncs(1) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * nfuncs(3) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * nfuncs(5),
        3 * s4 * c4 / 4 * nfuncs(3) + 3 * s6 * c6 / 4 * nfuncs(5),
        1.5 * c5 / l31 * nfuncs(4) - 1.5 * c4 / l23 * nfuncs(3),
        -nfuncs(2) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * nfuncs(3) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * nfuncs(4),
        3 * s4 * c4 / 4 * nfuncs(3) + 3 * s5 * c5 / 4 * nfuncs(4);
    return Hy_vector;
}

PlateStressData TriPlateElement::stress(
    Displacement d0, Displacement d1, Displacement d2, double L2, double L3)
{
    Eigen::Matrix3d tr = trans_matrix3(plane);

    Displacement d0t = d0.translate(tr);
    Displacement d1t = d1.translate(tr);
    Displacement d2t = d2.translate(tr);

    // Moments
    Eigen::Matrix3d Dmat = DMatrix();
    Eigen::VectorXd wvec(9);
    wvec << d0t.Dz(), d0t.Rx(), d0t.Ry(), d1t.Dz(), d1t.Rx(), d1t.Ry(), d2t.Dz(), d2t.Rx(), d2t.Ry();
    Eigen::Vector3d strs = Dmat * BMatrix(L2, L3) * wvec;

    // Shear
    Point p1 = plane.PointToCoord(Nodes[0]->Location);
    Point p2 = plane.PointToCoord(Nodes[1]->Location);
    Point p3 = plane.PointToCoord(Nodes[2]->Location);

    Vector v23 = p2 - p3;
    Vector v31 = p3 - p1;
    Vector v12 = p1 - p2;

    double L1 = 1 - L2 - L3;
    double area = Area();

    Eigen::VectorXd dLdx(3), dLdy(3);
    dLdx << v23.y * 0.5 / area, v31.y * 0.5 / area, v12.y * 0.5 / area;
    dLdy << -v23.x * 0.5 / area, -v31.x * 0.5 / area, -v12.x * 0.5 / area;

    // x方向の二階導関数
    Eigen::VectorXd d2ndx2(6);
    d2ndx2 << 4 * dLdx(0) * dLdx(0), // d²N₁/dx²
        4 * dLdx(1) * dLdx(1),       // d²N₂/dx²
        4 * dLdx(2) * dLdx(2),       // d²N₃/dx²
        8 * dLdx(2) * dLdx(1),       // d²N₄/dx²
        8 * dLdx(0) * dLdx(2),       // d²N₅/dx²
        8 * dLdx(1) * dLdx(0);       // d²N₆/dx²

    // y方向の二階導関数
    Eigen::VectorXd d2ndy2(6);
    d2ndy2 << 4 * dLdy(0) * dLdy(0), // d²N₁/dy²
        4 * dLdy(1) * dLdy(1),       // d²N₂/dy²
        4 * dLdy(2) * dLdy(2),       // d²N₃/dy²
        8 * dLdy(1) * dLdy(2),       // d²N₄/dy²
        8 * dLdy(0) * dLdy(2),       // d²N₅/dy²
        8 * dLdy(1) * dLdy(0);       // d²N₆/dy²

    // x,y方向の二階導関数
    Eigen::VectorXd d2ndxdy(6);
    d2ndxdy << 4 * dLdx(0) * dLdy(0),                // d²N₁/dxdy
        4 * dLdx(1) * dLdy(1),                       // d²N₂/dxdy
        4 * dLdx(2) * dLdy(2),                       // d²N₃/dxdy
        4 * (dLdx(1) * dLdy(2) + dLdy(1) * dLdx(2)), // d²N₄/dxdy
        4 * (dLdx(2) * dLdy(0) + dLdy(2) * dLdx(0)), // d²N₅/dxdy
        4 * (dLdx(0) * dLdy(1) + dLdy(0) * dLdx(1)); // d²N₆/dxdy

    double l12 = v12.norm();
    double l23 = v23.norm();
    double l31 = v31.norm();

    double c4 = -v23.y / l23;
    double c5 = -v31.y / l31;
    double c6 = -v12.y / l12;

    double s4 = v23.x / l23;
    double s5 = v31.x / l31;
    double s6 = v12.x / l12;

    Eigen::Vector<double, 9> d2Hx_dx2, d2Hx_dy2, d2Hx_dxdy, d2Hy_dx2, d2Hy_dy2, d2Hy_dxdy;
    d2Hx_dx2 = TriPlateElement_Hx_vector(s4, s5, s6, c4, c5, c6, l12, l23, l31, d2ndx2);
    d2Hx_dy2 = TriPlateElement_Hx_vector(s4, s5, s6, c4, c5, c6, l12, l23, l31, d2ndy2);
    d2Hx_dxdy = TriPlateElement_Hx_vector(s4, s5, s6, c4, c5, c6, l12, l23, l31, d2ndxdy);

    d2Hy_dx2 = TriPlateElement_Hy_vector(s4, s5, s6, c4, c5, c6, l12, l23, l31, d2ndx2);
    d2Hy_dy2 = TriPlateElement_Hy_vector(s4, s5, s6, c4, c5, c6, l12, l23, l31, d2ndy2);
    d2Hy_dxdy = TriPlateElement_Hy_vector(s4, s5, s6, c4, c5, c6, l12, l23, l31, d2ndxdy);

    Eigen::MatrixXd dBdx(3, 9), dBdy(3, 9);
    dBdx.row(0) = d2Hx_dx2;
    dBdx.row(1) = d2Hy_dxdy;
    dBdx.row(2) = d2Hx_dxdy + d2Hy_dx2;

    dBdy.row(0) = d2Hx_dxdy;
    dBdy.row(1) = d2Hy_dy2;
    dBdy.row(2) = d2Hx_dy2 + d2Hy_dxdy;

    Eigen::Vector3d dMdx = Dmat * (dBdx * wvec);
    Eigen::Vector3d dMdy = Dmat * (dBdy * wvec);

    double qx = dMdx(0) + dMdy(2);
    double qy = dMdy(1) + dMdx(2);

    // Plane Composition
    MembraneStressData mstr = plane_element.stress(d0, d1, d2);

    return PlateStressData(strs[0], strs[1], strs[2], qx, qy,
                           mstr.sigx * thickness.plane_thick, mstr.sigy * thickness.plane_thick, mstr.sigxy * thickness.plane_thick);
}

// PlateStressData TriPlateElement::stress_save(
//	Displacement d0, Displacement d1, Displacement d2, double xi, double eta)
//{
//	Eigen::Matrix3d tr = trans_matrix3(plane);
//
//	Displacement d0t = d0.translate(tr);
//	Displacement d1t = d1.translate(tr);
//	Displacement d2t = d2.translate(tr);
//
//	// Moments
//	Eigen::Matrix3d Dmat = DMatrix();
//	Eigen::VectorXd wvec(9);
//	wvec << d0t.Dz(), d0t.Rx(), d0t.Ry(), d1t.Dz(), d1t.Rx(), d1t.Ry(), d2t.Dz(), d2t.Rx(), d2t.Ry();
//	Eigen::Vector3d strs = Dmat * BMatrix(xi, eta) * wvec;
//
//	// Shear
//	Point p1 = plane.PointToCoord(Nodes[0]->Location);
//	Point p2 = plane.PointToCoord(Nodes[1]->Location);
//	Point p3 = plane.PointToCoord(Nodes[2]->Location);
//
//	Vector v23 = p2 - p3;
//	Vector v31 = p3 - p1;
//	Vector v12 = p1 - p2;
//
//	double p4 = -6 * v23.x / v23.squared_norm();
//	double p5 = -6 * v31.x / v31.squared_norm();
//	double p6 = -6 * v12.x / v12.squared_norm();
//
//	double q4 = 3 * v23.x * v23.y / v23.squared_norm();
//	double q5 = 3 * v31.x * v31.y / v31.squared_norm();
//	double q6 = 3 * v12.x * v12.y / v12.squared_norm();
//
//	double r4 = 3 * v23.y * v23.y / v23.squared_norm();
//	double r5 = 3 * v31.y * v31.y / v31.squared_norm();
//	double r6 = 3 * v12.y * v12.y / v12.squared_norm();
//
//	double t4 = -6 * v23.y / v23.squared_norm();
//	double t5 = -6 * v31.y / v31.squared_norm();
//	double t6 = -6 * v12.y / v12.squared_norm();
//
//	double Jl = v12.y * v31.x - v12.x * v31.y;
//	Eigen::Matrix2d Jinv;
//	Jinv << v31.y / Jl, v12.y / Jl, -v31.x / Jl, -v12.x / Jl;
//
//	Eigen::Vector<double, 9> Hy_xi_xi, Hy_xi_eta, Hy_eta_xi, Hy_eta_eta;
//	Hy_xi_xi << -2.0 * t6, -2.0 * r6, 2.0 * q6, 2.0 * t6, -2.0 * r6, 2.0 * q6, 0, 0, 0;
//	Hy_xi_eta << t5 - t6, -(r5 + r6), q5 + q6, t4 + t6, r4 - r6, -(q4 - q6), -(t4 + t5), r4 - r5, q5 - q4;
//	Hy_eta_xi << t5 - t6, -(r5 + r6), q5 + q6, t4 + t6, r4 - r6, -(q4 - q6), -(t4 + t5), r4 - r5, q5 - q4;
//	Hy_eta_eta << 2.0 * t5, -2.0 * r5, 2.0 * q5, 0, 0, 0, -2.0 * t5, -2.0 * r5, 2.0 * q5;
//
//	Eigen::Vector<double, 9> Hx_xi_xi, Hx_xi_eta, Hx_eta_xi, Hx_eta_eta;
//	Hx_xi_xi << -2.0 * p6, -2.0 * q6, 6.0 - 2.0 * r6, 2.0 * p6, -2.0 * q6, 6.0 - 2.0 * r6, 0, 0, 0;
//	Hx_xi_eta << p5 - p6, -(q5 + q6), 6.0 - (r5 + r6), p4 + p6, q4 - q6, r4 - r6, -(p5 + p4), (q4 - q5), -(r5 - r4);
//	Hx_eta_xi << p5 - p6, -(q5 + q6), 6.0 - (r5 + r6), p4 + p6, q4 - q6, r4 - r6, -(p4 + p5), q4 - q5, r4 - r5;
//	Hx_eta_eta << 2.0 * p5, -2.0 * q5, 6.0 - 2.0 * r5, 0, 0, 0, -2.0 * p5, -2.0 * q5, 6.0 - 2.0 * r5;
//
//	Eigen::Vector<double, 9> Hy_xi_x, Hy_xi_y, Hy_eta_x, Hy_eta_y;
//	Hy_xi_x = Jinv(0, 0) * Hy_xi_xi + Jinv(0, 1) * Hy_xi_eta;
//	Hy_xi_y = Jinv(1, 0) * Hy_xi_xi + Jinv(1, 1) * Hy_xi_eta;
//	Hy_eta_x = Jinv(0, 0) * Hy_eta_xi + Jinv(0, 1) * Hy_eta_eta;
//	Hy_eta_y = Jinv(1, 0) * Hy_eta_xi + Jinv(1, 1) * Hy_eta_eta;
//
//	Eigen::Vector<double, 9> Hx_xi_x, Hx_xi_y, Hx_eta_x, Hx_eta_y;
//	Hx_xi_x = Jinv(0, 0) * Hx_xi_xi + Jinv(0, 1) * Hx_xi_eta;
//	Hx_xi_y = Jinv(1, 0) * Hx_xi_xi + Jinv(1, 1) * Hx_xi_eta;
//	Hx_eta_x = Jinv(0, 0) * Hx_eta_xi + Jinv(0, 1) * Hx_eta_eta;
//	Hx_eta_y = Jinv(1, 0) * Hx_eta_xi + Jinv(1, 1) * Hx_eta_eta;
//
//	Eigen::MatrixXd Bmat_y(3, 9);
//	Bmat_y.row(0) = v31.y * Hx_xi_y + v12.y * Hx_eta_y;
//	Bmat_y.row(1) = -v31.x * Hy_xi_y - v12.x * Hy_eta_y;
//	Bmat_y.row(2) = -v31.x * Hx_xi_y - v12.x * Hx_eta_y + v31.y * Hy_xi_y + v12.y * Hy_eta_y;
//	Bmat_y /= (2.0 * Area());
//
//	Eigen::MatrixXd Bmat_x(3, 9);
//	Bmat_x.row(0) = v31.y * Hx_xi_x + v12.y * Hx_eta_x;
//	Bmat_x.row(1) = -v31.x * Hy_xi_x - v12.x * Hy_eta_x;
//	Bmat_x.row(2) = -v31.x * Hx_xi_x - v12.x * Hx_eta_x + v31.y * Hy_xi_x + v12.y * Hy_eta_x;
//	Bmat_x /= (2.0 * Area());
//
//	Eigen::Vector3d dMdx = Dmat * (Bmat_x * wvec);
//	Eigen::Vector3d dMdy = Dmat * (Bmat_y * wvec);
//
//	double qx = dMdx(0) + dMdy(2);
//	double qy = dMdy(1) + dMdx(2);
//
//	// Plane Composition
//	MembraneStressData mstr = plane_element.stress(d0, d1, d2);
//
//	return PlateStressData(strs[0], strs[1], strs[2], qx, qy,
//		mstr.sigx * thickness.plane_thick, mstr.sigy * thickness.plane_thick, mstr.sigxy * thickness.plane_thick);
// }
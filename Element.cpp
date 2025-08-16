//#ifdef USE_MKL
//#define EIGEN_USE_MKL_ALL
//#endif
#include <Eigen/Sparse>
#include <Eigen/LU>

#include "Components.h"
#include "Element.h"

//Eigen::Matrix3d ElementBase::trans_matrix3(const Plane &plane)
//{
//	Eigen::Matrix3d mat;
//	mat << plane.ex.x, plane.ex.y, plane.ex.z,
//		plane.ey.x, plane.ey.y, plane.ey.z,
//		plane.ez.x, plane.ez.y, plane.ez.z;
//
//	return mat;
//}

Eigen::Matrix3d trans_matrix3(const Plane plane)
{
	Eigen::Matrix3d m;
	m << plane.ex.x, plane.ex.y, plane.ex.z,
		plane.ey.x, plane.ey.y, plane.ey.z,
		plane.ez.x, plane.ez.y, plane.ez.z;
	return m;
}


Eigen::Matrix3d trans_matrix3(const Point p0, const Point p1, const double beta)
{
	const double TOL = 0.0001;
	Vector ex, ey, ez;
	//Point p0, p1;
	//p0 = Nodes[0]->Location;
	//p1 = Nodes[1]->Location;

	ex = Point::subtract(p1, p0);
	// ex = p1 - p0;
	ex = ex * (1 / ex.norm());

	// 材軸方向と全体座標系Z軸が並行のとき
	double dot = ex * Vector::ZAxis();
	if (std::abs(1 - std::abs(dot)) < TOL)
	{
		if (dot < 0)
			ey = Vector::YAxis();
		else
			ey = Vector::YAxis() * -1.0;
		
		// Ref Book
		//ey = Vector::XAxis();
		//ez = Vector::YAxis();

		// Midas仕様
		ez = Vector::XAxis();
	}
	// 材軸方向と全体座標系Z軸が並行ではないとき
	else {
		ez = Vector::ZAxis() - ex * (ex * Vector::ZAxis());
		ez = ez * (1 / ez.norm());
		ey = Vector::cross(ez, ex);
	}

	// Rotate
	Vector ey2 = ey * std::cos(beta) + ez * std::sin(beta);
	Vector ez2 = ez * std::cos(beta) - ey * std::sin(beta);

	Eigen::Matrix3d m;
	m << ex.x, ex.y, ex.z,
		ey2.x, ey2.y, ey2.z,
		ez2.x, ez2.y, ez2.z;
	return m;
}

//Eigen::Matrix3d trans_matrix3(const Plane plane)
//{
//	Eigen::Matrix3d m;
//	m << plane.ex.x, plane.ex.y, plane.ex.z,
//		plane.ey.x, plane.ey.y, plane.ey.z,
//		plane.ez.x, plane.ez.y, plane.ez.z;
//	return m;
//}

//Eigen::MatrixXd stiff_matrix_local(
//	double E, double G, double l, double A, double Iy, double Iz, double K)
//{
//	double EAl = E * A / l;
//	double GKl = G * K / l;
//
//	double EIz12 = 12 * E * Iz / (l * l * l);
//	double EIz6 = 6 * E * Iz / (l * l);
//	double EIz2 = 2 * E * Iz / l;
//	double EIz4 = EIz2 * 2;
//
//	double EIy12 = 12 * E * Iy / (l * l * l);
//	double EIy6 = 6 * E * Iy / (l * l);
//	double EIy2 = 2 * E * Iy / l;
//	double EIy4 = EIy2 * 2;
//
//	Eigen::MatrixXd X(12, 12);
//	X << EAl, 0, 0, 0, 0, 0, -EAl, 0, 0, 0, 0, 0,
//		0, EIz12, 0, 0, 0, EIz6, 0, -EIz12, 0, 0, 0, EIz6,
//		0, 0, EIy12, 0, EIy6, 0, 0, 0, -EIy12, 0, EIy6, 0,
//		0, 0, 0, GKl, 0, 0, 0, 0, 0, -GKl, 0, 0,
//		0, 0, EIy6, 0, EIy4, 0, 0, 0, -EIy6, 0, EIy2, 0,
//		0, EIz6, 0, 0, 0, EIz4, 0, -EIz6, 0, 0, 0, EIz2,
//		-EAl, 0, 0, 0, 0, 0, EAl, 0, 0, 0, 0, 0,
//		0, -EIz12, 0, 0, 0, -EIz6, 0, EIz12, 0, 0, 0, -EIz6,
//		0, 0, -EIy12, 0, -EIy6, 0, 0, 0, EIy12, 0, -EIy6, 0,
//		0, 0, 0, -GKl, 0, 0, 0, 0, 0, GKl, 0, 0,
//		0, 0, EIy6, 0, EIy2, 0, 0, 0, -EIy6, 0, EIy4, 0,
//		0, EIz6, 0, 0, 0, EIz2, 0, -EIz6, 0, 0, 0, EIz4;
//
//	return X;
//}





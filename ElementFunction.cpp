#include <Eigen/Dense>

#include "ElementFunction.h"



Eigen::MatrixXd stiff_matrix_local(
	double E, double G, double l, double A, double Iy, double Iz, double K)
{
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

	Eigen::MatrixXd X(12, 12);
	X << EAl, 0, 0, 0, 0, 0, -EAl, 0, 0, 0, 0, 0,
		 0, EIz12, 0, 0, 0, EIz6, 0, -EIz12, 0, 0, 0, EIz6,
		 0, 0, EIy12, 0, EIy6, 0, 0, 0, -EIy12, 0, EIy6, 0,
		 0, 0, 0, GKl, 0, 0, 0, 0, 0, -GKl, 0, 0,
		 0, 0, EIy6, 0, EIy4, 0, 0, 0, -EIy6, 0, EIy2, 0,
		 0, EIz6, 0, 0, 0, EIz4, 0, -EIz6, 0, 0, 0, EIz2,
		 -EAl, 0, 0, 0, 0, 0, EAl, 0, 0, 0, 0, 0,
		 0, -EIz12, 0, 0, 0, -EIz6, 0, EIz12, 0, 0, 0, -EIz6,
		 0, 0, -EIy12, 0, -EIy6, 0, 0, 0, EIy12, 0, -EIy6, 0,
		 0, 0, 0, -GKl, 0, 0, 0, 0, 0, GKl, 0, 0,
		 0, 0, EIy6, 0, EIy2, 0, 0, 0, -EIy6, 0, EIy4, 0,
		 0, EIz6, 0, 0, 0, EIz2, 0, -EIz6, 0, 0, 0, EIz4;
	
	return X;
}
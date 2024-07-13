#ifndef _ELEMENT_FUNCTION_
#define _ELEMENT_FUNCTION_

Eigen::MatrixXd stiff_matrix_local(
	double E, double G, double l, double A, double Iy, double Iz, double K);

#endif
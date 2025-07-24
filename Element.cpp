//#ifdef USE_MKL
//#define EIGEN_USE_MKL_ALL
//#endif
#include <Eigen/Sparse>
#include <Eigen/LU>

#include "Components.h"
#include "Element.h"


Eigen::Matrix3d trans_matrix_from_plane(const Plane& plane)
{
	Eigen::Matrix3d mat;
	mat << plane.ex.x, plane.ex.y, plane.ex.z,
		plane.ey.x, plane.ey.y, plane.ey.z,
		plane.ez.x, plane.ez.y, plane.ez.z;

	return mat;
}

Eigen::Vector<double, 4> ShapeFunction4(double xi, double eta)
{
	Eigen::Vector<double, 4> Nfunc;
	Nfunc << 0.25 * (1 - xi) * (1 - eta),
		0.25 * (1 + xi) * (1 - eta),
		0.25 * (1 + xi) * (1 + eta),
		0.25 * (1 - xi) * (1 + eta);
	return Nfunc;
}

Eigen::Vector<double, 8> ShapeFunctionSerendipity8(double xi, double eta)
{
	Eigen::Vector<double, 8> Nfunc;
	Nfunc << 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1),
		0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1),
		0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1),
		0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1),
		0.5 * (1 - xi * xi) * (1 - eta),
		0.5 * (1 + xi) * (1 - eta * eta),
		0.5 * (1 - xi * xi) * (1 + eta),
		0.5 * (1 - xi) * (1 - eta * eta);
	return Nfunc;
}

Eigen::Vector<double, 6> ShapeFunctionTriangle6(double xi, double eta)
{
	double gamma = 1.0 - xi - eta;
	Eigen::Vector<double, 6> Nfunc;
	Nfunc << 2.0 * xi * xi - xi,
		2.0 * eta * eta - eta,
		2.0 * gamma * gamma - gamma,
		4.0 * eta * gamma,
		4.0 * gamma * xi,
		4.0 * xi * eta;
	return Nfunc;
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
	//double alpha = vij * Vector::XAxis() / l;
	//double beta = vij * Vector::YAxis() / l;
	//double gamma = vij * Vector::ZAxis() / l;

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

	Kg_mat << 
		1, 0, 0, -1, 0, 0,
		0, 1, 0, 0, -1, 0,
		0, 0, 1, 0, 0, -1,
		-1, 0, 0, 1, 0, 0,
		0, -1, 0, 0, 1, 0,
		0, 0, -1, 0, 0, 1;

	Kg_mat = Kg_mat * Nx / length;
	//Eigen::MatrixXd transMat = trans_matrix();
	//return transMat.transpose() * Kg_mat * transMat;
	return Kg_mat;
}

TrussElement::TrussElement(Node* n0, Node* n1, Section* sec, Material mat)
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
	//Vector vij = ;
	//double l = (Nodes[1]->Location - Nodes[0]->Location).norm();
	//double EAl = Mat->Young * Sec->A / l;
	//double alpha = vij * Vector::XAxis() / l;
	//double beta = vij * Vector::YAxis() / l;
	//double gamma = vij * Vector::ZAxis() / l;

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


BeamElement::BeamElement(Node* n0, Node* n1, Section* sec, Material mat, double beta):
	BeamElement()
{
	Nodes[0] = n0;
	Nodes[1] = n1;

	Sec = sec;
	Mat = mat;

	Beta = beta;
}

Eigen::MatrixXd BeamElement::StiffnessMatrix() {
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
	Eigen::Vector<double, total_dof> f; f.setZero();

	f << accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0,
		accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0;

	f = m * f;

	std::cout << "Mass Matrix: " << m.diagonal() << std::endl;
	//std::cout << "BodyforceToNodeLoadData: " << f.transpose() << std::endl;

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
	//m << l / 3, 0, 0, 0, 0, 0, l / 6, 0, 0, 0, 0, 0,
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
	//wvec = StiffnessMatrix() * wvec;

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
	
	m << 1 - p, 0, 0, 0, 0, 0, p, 0, 0, 0, 0, 0, // dx
		0, 1-3*p2+2*p3, 0, 0, 0, (1-2*p+p2)*x, 0, 3*p2-2*p3, 0, 0, 0, (-p+p2)*x, // dy v
		0, 0, 1-3*p2+2*p3, 0, (-1+2*p-p2)*x, 0, 0, 0, 3*p2-2*p3, 0, (p-p2)*x, 0, // dz w
		0, 0, 0, 1 - p, 0, 0, 0, 0, 0, p, 0, 0, // rx
		0, 0, -6.0 * (p2 - p) / l, 0, 1.0 - 4.0 * p + 3.0 * p2, 0, 0, 0, -6.0 * (p - p2) / l, 0, -2.0 * p + 3.0 * p2, 0, // ry
		0, 6.0 * (p2 - p)/l, 0, 0, 0, 1.0-4.0*p+3.0*p2, 0, 6.0 * (p-p2)/l, 0, 0, 0, 3.0 * p2 - 2.0 * p; // rz

	Eigen::VectorXd d = T.block(0, 0, 6, 6).transpose() * m * wvec;

	return Displacement(d(0), d(1), d(2), d(3), d(4), d(5));
}

double BeamElement::element_length()
{
	Point p0 = Nodes[0]->Location;
	Point p1 = Nodes[1]->Location;
	double l2 = (p1.x - p0.x)* (p1.x - p0.x) + (p1.y - p0.y) * (p1.y - p0.y) + (p1.z - p0.z) * (p1.z - p0.z);
	return sqrt(l2);
}

Eigen::MatrixXd BeamElement::stiffness_matrix_local() {
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

Eigen::MatrixXd BeamElement::geometric_local_stiffness_matrix(const std::vector<Displacement> &disp){
	Eigen::MatrixXd Kg_mat = Eigen::MatrixXd::Zero(total_dof, total_dof);
	double l = element_length();
	BeamStress s = stress(disp[0], disp[1]);
	double Nx = s.S0.Nx;
	double Qy = s.S0.Qy;
	double Qz = s.S0.Qz;

	int indices_b[8]{ 1,2,4,5,7,8,10,11 };
	int indices_x[6]{ 0,1,2,6,7,8 };

	Eigen::MatrixXd Kg_b = Eigen::MatrixXd::Zero(8, 8);
	Kg_b << 6.0/(5.0*l), 0, 0, 1.0/10.0, -6.0/(5.0*l), 0, 0, 1.0/10.0,
		0, 6.0/(5.0*l), -1.0/10.0, 0, 0, -6.0/(5.0*l), -1.0/10.0, 0,
		0, -1.0/10.0, 2.0*l/15.0, 0, 0, 1.0/10.0, -l/30.0, 0,
		1.0/10.0, 0, 0, 2.0*l/15.0, -1.0/10.0, 0, 0, -l/30.0,
		-6.0/(5.0*l), 0, 0, -1.0/10.0, 6.0/(5.0*l), 0, 0, -1.0/10.0,
		0, -6.0/(5.0*l), 1.0/10.0, 0, 0, 6.0/(5.0*l), 1.0/10.0, 0,
		0, -1.0/10.0, -l/30.0, 0, 0, 1.0/10.0, 2.0*l/15.0, 0,
		1.0/10.0, 0, 0, -l/30.0, -1.0/10.0, 0, 0, 2.0*l/15.0;
	 Kg_b *= Nx;

	Eigen::MatrixXd Kg_x = Eigen::MatrixXd::Zero(6, 6);
	Kg_x << 1, 0, 0, -1, 0, 0,
			0, 1, 0, 0, -1, 0,
			0, 0, 1, 0, 0, -1,
			-1, 0, 0, 1, 0, 0,
			0, -1, 0, 0, 1, 0,
			0, 0, -1, 0, 0, 1;
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

	for (size_t i = 0; i < 6; i++)
	{
		int ir = indices_x[i];
		for (size_t j = 0; j < 6; j++)
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
	//Eigen::MatrixXd Kg_mat1 = Eigen::MatrixXd::Zero(total_dof, total_dof);
	//Kg_mat1 << Nx / l, -(s.S0.Mz - s.S1.Mz) / (l * l), (s.S0.My - s.S1.My) / (l * l), 0, s.S0.My / l, -s.S0.Mz / l, -Nx / l, (s.S1.Mz - s.S0.Mz) / (l * l), (s.S1.My - s.S0.My) / (l * l), 0, -s.S1.My / l, s.S1.Mz / l,
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

	//Eigen::MatrixXd Kg_mat1 = Eigen::MatrixXd::Zero(total_dof, total_dof);
	//Kg_mat1 << Nx / l, (s.S0.Mz - s.S1.Mz) / (l * l), (s.S0.My - s.S1.My) / (l * l), 0, s.S0.My / l, -s.S0.Mz / l, -Nx / l, (s.S1.Mz - s.S0.Mz) / (l * l), (s.S1.My - s.S0.My) / (l * l), 0, -s.S1.My / l, s.S1.Mz / l,
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
	//Eigen::MatrixXd Kg_mat2 = Eigen::MatrixXd::Zero(total_dof, total_dof);
	//Kg_mat2 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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
	//Eigen::MatrixXd Kg_mat3 = Eigen::MatrixXd::Zero(total_dof, total_dof);
	//Kg_mat3 << 0, 0, -1.0/l, 0, 0, 0, 0, 0, 1.0/l, 0, 0, 0,
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
	//Kg_mat3 *= Qz;
	//
	// A1 A4^t + A4 A1^t + A3 A5^t + A5 A3^t
	//Eigen::MatrixXd Kg_mat4 = Eigen::MatrixXd::Zero(total_dof, total_dof);
	//Kg_mat4 << 0, -1.0/l, 0, 0, 0, 0, 0, 1.0/l, 0, 0, 0, 0,
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
	//Kg_mat4 *= Qy;

	//Kg_mat = Kg_mat1 + Kg_mat2 + Kg_mat3 + Kg_mat4;
	//Kg_mat *= 0.5;

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

Eigen::Matrix3d trans_matrix3(const Plane plane)
{
	Eigen::Matrix3d m;
	m << plane.ex.x, plane.ex.y, plane.ex.z,
		plane.ey.x, plane.ey.y, plane.ey.z,
		plane.ez.x, plane.ez.y, plane.ez.z;
	return m;
}

Eigen::MatrixXd BeamElement::trans_matrix()
{
	
	Eigen::Matrix3d tr0 = trans_matrix3(Nodes[0]->Location, Nodes[1]->Location, Beta);
	
	// ブロックの数を指定
	const int numBlocks = 4;

	// 大きな行列を作成し、対角ブロック行列として同じ行列を配置
	Eigen::MatrixXd matrix(total_dof, total_dof);
	matrix.setZero();  // 行列を0で初期化

	for (int i = 0; i < numBlocks; ++i) {
		matrix.block(i * tr0.rows(), i * tr0.cols(), tr0.rows(), tr0.cols()) = tr0;
	}

	return matrix;
}

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

std::ostream& operator<<(std::ostream& os, const BeamStressData& bsd)
{
	/*std::string force = std::format("Nx: {}, Qy: {}, Qz: {}, Mx: {}, My: {}, Mz: {}", 
		bsd.Nx, bsd.Qy, bsd.Qz, bsd.Mx, bsd.My, bsd.Mz);
	os << force << std::endl;*/

	os << "Nx: " << bsd.Nx << ", Qy: " << bsd.Qy << ", Qz: " << bsd.Qz
		<< ", Mx: " << bsd.Mx << ", My: " << bsd.My << ", Mz: " << bsd.Mz << std::endl;
	return os;
}

std::ostream& operator<<(std::ostream& os, const BeamStress& bsd)
{
	os << "Beam Stress" << std::endl;
	os << "Start  " << bsd.S0;
	os << "End    " << bsd.S1;
	return os;
}

std::ostream& operator<<(std::ostream& os, const MembraneStressData& strs)
{
	os << "Membrane Stress  " << "sig x: " << strs.sigx << "sig y: " << strs.sigy
		<< "sig xy: " << strs.sigxy;
	return os;
}

std::ostream& operator<<(std::ostream& os, const PlateStressData& strs)
{
	os << "Plate Stress Mx: " << strs.Mx << ", My: " << strs.My
		<< ", Mxy: " << strs.Mxy << ", Qx: " << strs.Qx << ", Qy: " << strs.Qy;
	return os;
}

TriPlaneElement::TriPlaneElement(Node* n0, Node* n1, Node* n2, 
	double t, Material mat, double beta)
{
	Nodes[0] = n0;
	Nodes[1] = n1;
	Nodes[2] = n2;
	
	thickness = Thickness(t);
	Mat = mat;
	plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
	plane.Rotate(beta, Vector::ZAxis());
}

TriPlaneElement::TriPlaneElement(Node* n0, Node* n1, Node* n2, Thickness t, Material mat, double beta)
{
	Nodes[0] = n0;
	Nodes[1] = n1;
	Nodes[2] = n2;

	thickness = t;
	Mat = mat;
	plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
	plane.Rotate(beta, Vector::ZAxis());
}

//Eigen::Matrix2d TriPlaneElement::invJMatrix()
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
//}

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
	Eigen::Vector<double, total_dof> f; f.setZero();

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

void TriPlaneElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat)
{
	Eigen::VectorXd mass = NodeLumpedMass();
	for (size_t ni = 0; ni < node_num; ni++)
	{
		for (size_t i = 0; i < 3; i++) {
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

Eigen::MatrixXd TriPlaneElement::geometric_local_stiffness_matrix(const std::vector<Displacement>& disp)
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
	//Eigen::MatrixXd tr = trans_matrix();
	//return tr.transpose() * Kg * tr;
	return Kg;
}

Eigen::MatrixXd TriPlaneElement::trans_matrix()
{	
	Eigen::Matrix3d tr0 = trans_matrix3(plane);

	// ブロックの数を指定
	const int numBlocks = node_num;

	// 大きな行列を作成し、対角ブロック行列として同じ行列を配置
	Eigen::MatrixXd matrix(3*2, total_dof);
	matrix.setZero();  // 行列を0で初期化

	for (int i = 0; i < numBlocks; ++i) {
		matrix.block(i * 2, i * node_dof, 2, node_dof) = tr0.block(0,0,2,node_dof);
	}

	return matrix;
}

//Eigen::MatrixXd TriPlateElement::BMatrix_Save(double xi, double eta)
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
//}

// 入力形状関数に応じたx,yそれぞれの方向の列ベクトルによるHベクトル行列
Eigen::MatrixXd TriPlateElement::HVecs(Eigen::Vector<double, 6> shape_funcs) {
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
	//double L1 = 1 - L2 - L3;
	//double area = Area();

	//Eigen::VectorXd dLdx(3), dLdy(3);
	//dLdx << v23.y * 0.5 / area, v31.y * 0.5 / area, v12.y * 0.5 / area;
	//dLdy << -v23.x * 0.5 / area, -v31.x * 0.5 / area, -v12.x * 0.5 / area;

	//Eigen::VectorXd dndx(6), dndy(6);
	//dndx << dLdx(0) * (4 * L1 - 1), dLdx(1)* (4 * L2 - 1), dLdx(2)* (4 * L3 - 1),
	//	4 * (dLdx(1) * L3 + L2 * dLdx(2)), 4 * (dLdx(2) * L1 + L3 * dLdx(0)),
	//	4 * (dLdx(0) * L2 + L1 * dLdx(1));
	//dndy << dLdy(0) * (4 * L1 - 1), dLdy(1)* (4 * L2 - 1), dLdy(2)* (4 * L3 - 1),
	//	4 * (dLdy(1) * L3 + L2 * dLdy(2)), 4 * (dLdy(2) * L1 + L3 * dLdy(0)),
	//	4 * (dLdy(0) * L2 + L1 * dLdy(1));
	// -------------------------------------------------------

	Eigen::Matrix<double, 2, 9> HVecs;

	// Hx
	HVecs.row(0) <<
		1.5 * s5 / l31 * shape_funcs(4) - 1.5 * s6 / l12 * shape_funcs(5),
		-3 * s5 * c5 / 4 * shape_funcs(4) - 3 * s6 * c6 / 4 * shape_funcs(5),
		shape_funcs(0) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * shape_funcs(4) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * shape_funcs(5),
		1.5 * s6 / l12 * shape_funcs(5) - 1.5 * s4 / l23 * shape_funcs(3),
		-3 * s4 * c4 / 4 * shape_funcs(3) - 3 * s6 * c6 / 4 * shape_funcs(5),
		shape_funcs(1) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * shape_funcs(3) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * shape_funcs(5),
		1.5 * s4 / l23 * shape_funcs(3) - 1.5 * s5 / l31 * shape_funcs(4),
		-3 * s4 * c4 / 4 * shape_funcs(3) - 3 * s5 * c5 / 4 * shape_funcs(4),
		shape_funcs(2) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * shape_funcs(3) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * shape_funcs(4);

	// Hy
	HVecs.row(1) <<
		1.5 * c6 / l12 * shape_funcs(5) - 1.5 * c5 / l31 * shape_funcs(4),
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
	dndx << dLdx(0) * (4 * L1 - 1), dLdx(1)* (4 * L2 - 1), dLdx(2)* (4 * L3 - 1),
		4 * (dLdx(1) * L3 + L2 * dLdx(2)), 4 * (dLdx(2) * L1 + L3 * dLdx(0)),
		4 * (dLdx(0) * L2 + L1 * dLdx(1));
	dndy << dLdy(0) * (4 * L1 - 1), dLdy(1)* (4 * L2 - 1), dLdy(2)* (4 * L3 - 1),
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
	dHx_dx <<
		1.5 * s5 / l31 * dndx(4) - 1.5 * s6 / l12 * dndx(5),
		-3 * s5 * c5 / 4 * dndx(4) - 3 * s6 * c6 / 4 * dndx(5),
		dndx(0) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * dndx(4) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * dndx(5),
		1.5 * s6 / l12 * dndx(5) - 1.5 * s4 / l23 * dndx(3),
		-3 * s4 * c4 / 4 * dndx(3) - 3 * s6 * c6 / 4 * dndx(5),
		dndx(1) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * dndx(3) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * dndx(5),
		1.5 * s4 / l23 * dndx(3) - 1.5 * s5 / l31 * dndx(4),
		-3 * s4 * c4 / 4 * dndx(3) - 3 * s5 * c5 / 4 * dndx(4),
		dndx(2) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * dndx(3) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * dndx(4);

	dHx_dy <<
		1.5 * s5 / l31 * dndy(4) - 1.5 * s6 / l12 * dndy(5),
		-3 * s5 * c5 / 4 * dndy(4) - 3 * s6 * c6 / 4 * dndy(5),
		dndy(0) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * dndy(4) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * dndy(5),
		1.5 * s6 / l12 * dndy(5) - 1.5 * s4 / l23 * dndy(3),
		-3 * s4 * c4 / 4 * dndy(3) - 3 * s6 * c6 / 4 * dndy(5),
		dndy(1) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * dndy(3) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * dndy(5),
		1.5 * s4 / l23 * dndy(3) - 1.5 * s5 / l31 * dndy(4),
		-3 * s4 * c4 / 4 * dndy(3) - 3 * s5 * c5 / 4 * dndy(4),
		dndy(2) + (c4 * c4 * 0.5 - s4 * s4 * 0.25) * dndy(3) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * dndy(4);

	dHy_dx <<
		1.5 * c6 / l12 * dndx(5) - 1.5 * c5 / l31 * dndx(4),
		-dndx(0) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * dndx(4) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * dndx(5),
		3 * s5 * c5 / 4 * dndx(4) + 3 * s6 * c6 / 4 * dndx(5),
		1.5 * c4 / l23 * dndx(3) - 1.5 * c6 / l12 * dndx(5),
		-dndx(1) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * dndx(3) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * dndx(5),
		3 * s4 * c4 / 4 * dndx(3) + 3 * s6 * c6 / 4 * dndx(5),
		1.5 * c5 / l31 * dndx(4) - 1.5 * c4 / l23 * dndx(3),
		-dndx(2) + (c4 * c4 * 0.25 - s4 * s4 * 0.5) * dndx(3) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * dndx(4),
		3 * s4 * c4 / 4 * dndx(3) + 3 * s5 * c5 / 4 * dndx(4);
		
	dHy_dy <<
		1.5 * c6 / l12 * dndy(5) - 1.5 * c5 / l31 * dndy(4),
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
	//mat /= (2.0 * Area());
	return mat;
}

//void TriPlateElement::shearstress(Displacement d0, Displacement d1, Displacement d2)
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
//}


Eigen::Matrix3d TriPlateElement::DMatrix()
{
	Eigen::Matrix3d mat;
	mat << 1, Mat.Poisson, 0,
		Mat.Poisson, 1, 0,
		0, 0, (1 - Mat.Poisson) / 2;
	return Mat.Young * pow(thickness.plate_thick, 3) / 12 / (1 - Mat.Poisson * Mat.Poisson) * mat;
}

//Eigen::MatrixXd TriPlateElement::localStiffnessMatrix_Save()
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
//}

Eigen::MatrixXd TriPlateElement::localStiffnessMatrix()
{
	double weight = 1.0 / 3.0;
	double xi_list[3]{ 0.5 , 0.5 , 0 };
	double eta_list[3]{ 0, 0.5, 0.5 };

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

	//Eigen::MatrixXd trMat = trans_matrix();
	//return trMat.transpose() * Kg * trMat;
	return Kg;
}

TriPlateElement::TriPlateElement(Node* n0, Node* n1, Node* n2, Thickness t, Material mat)
	: plane_element(n0, n1, n2, t, mat)
{
	Nodes[0] = n0;
	Nodes[1] = n1;
	Nodes[2] = n2;

	Mat = mat;
	thickness = t;
	plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
}

TriPlateElement::TriPlateElement(Node* n0, Node* n1, Node* n2, double t, Material mat)
	: plane_element(n0, n1, n2, t, mat)
{
	Nodes[0] = n0;
	Nodes[1] = n1;
	Nodes[2] = n2;

	Mat = mat;
	thickness = Thickness(t);
	plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
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

//Eigen::MatrixXd TriPlateElement::NodeConsistentMass()
//{
//	return Eigen::MatrixXd::Identity(total_dof, total_dof);
//}

std::vector<NodeLoadData> TriPlateElement::InertialForceToNodeLoadData(Eigen::Vector3d accel_vec)
{
	Eigen::Matrix<double, total_dof, total_dof> m;
	m = NodeConsistentMass();
	Eigen::Vector<double, total_dof> f; f.setZero();

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
	//Eigen::Vector<double, 1> intg_weights; intg_weights(0) = 1.0;
	//Eigen::Vector<double, 1> xi_params; xi_params(0) = 1.0 / 3.0;
	//Eigen::Vector<double, 1> eta_params; eta_params(0) = 1.0 / 3.0;

	// 3-Point Gauss Quadrature
	//const Eigen::Vector3d intg_weights(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
	//const Eigen::Vector3d xi_params(0.5, 0.5, 0);
	//const Eigen::Vector3d eta_params(0, 0.5, 0.5);

	// 4-Point Gauss Quadrature
	//const Eigen::Vector4d intg_weights(-27.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0);
	//const Eigen::Vector4d xi_params(1.0 / 3.0, 0.6, 0.2, 0.2);
	//const Eigen::Vector4d eta_params(1.0 / 3.0, 0.2, 0.6, 0.2);

	// 7-Point Gauss Quadrature(この精度が必要になりそう)
	Eigen::Vector<double, 7> intg_weights;
	Eigen::Vector<double, 7> xi_params;
	Eigen::Vector<double, 7> eta_params;
	intg_weights << 0.225, 0.1323941527, 0.1323941527, 0.1323941527, 0.1259391805, 0.1259391805, 0.1259391805;
	xi_params << 1.0 / 3.0, 0.0597158717, 0.4701420641, 0.4701420641, 0.7974269853, 0.1012865073, 0.1012865073;
	eta_params << 1.0 / 3.0, 0.4701420641, 0.0597158717, 0.4701420641, 0.1012865073, 0.7974269853, 0.1012865073;
	
	//// 9-Point Gauss Quadrature
	//Eigen::VectorXd intg_weights(9);
	//intg_weights << 0.2059505047608870, 0.2059505047608870, 0.2059505047608870,
	//	0.0636914142862230, 0.0636914142862230, 0.0636914142862230,
	//	0.0636914142862230, 0.0636914142862230, 0.0636914142862230;
	//Eigen::VectorXd xi_params(9);
	//xi_params << 0.1249495032332320, 0.4375252483838400, 0.4375252483838400,
	//	0.7971126518600710, 0.1654099273984100, 0.0374774207500880,
	//	0.7971126518600710, 0.0374774207500880, 0.1654099273984100;
	//Eigen::VectorXd eta_params(9);
	//eta_params << 0.4375252483838400, 0.1249495032332320, 0.4375252483838400,
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
	//std::cout << "Mp: \n" << Mp << std::endl;
	//std::cout << "area: " << area << ", " << Mat.dense << ", " << thickness.plane_thick << std::endl;
	//Eigen::Vector<int, 9> shape_indices;
	//shape_indices << 0, 1, 2, 6, 7, 8, 12, 13, 14;
	const int shape_indices[9]{ 0, 1, 2, 6, 7, 8, 12, 13, 14 };

	for (size_t j = 0; j < 9; j++)
		for (size_t k = 0; k < 9; k++)
			M(shape_indices[j], shape_indices[k]) += Mp(j, k);

	//Eigen::Vector<int, 9> HVecs_indices;
	//HVecs_indices << 2, 3, 4, 8, 9, 10, 14, 15, 16;
	const int HVecs_indices[9]{ 2, 3, 4, 8, 9, 10, 14, 15, 16 };

	for (size_t ix = 0; ix < xi_params.size(); ix++)
	{
		// wi, theta_xi, theta_yi...
		Eigen::MatrixXd H_mat = HVecs(ShapeFunctionTriangle6(xi_params(ix), eta_params(ix)));
		//std::cout << "sum of shape func: " << ShapeFunctionSerendipity8(intg_params(ix), intg_params(iy)).sum() << std::endl;
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

double TriPlateElement::Area() {
	Vector v01 = Nodes[1]->Location - Nodes[0]->Location;
	Vector v02 = Nodes[2]->Location - Nodes[0]->Location;
	Vector vn = Vector::cross(v01, v02);
	return vn.norm() / 2;
}

Eigen::MatrixXd TriPlateElement::StiffnessMatrix()
{
	Eigen::MatrixXd Kpln = plane_element.localStiffnessMatrix();
	Eigen::MatrixXd Kplt = localStiffnessMatrix();
	//Eigen::MatrixXd Kplt = localStiffnessMatrix();
	//std::cout << "K plane: \n" << Kpln << std::endl;
	//std::cout << "K plate: \n" << Kplt << std::endl;
	double max_diag = std::max(Kpln.diagonal().maxCoeff(), Kplt.diagonal().maxCoeff());

	Eigen::Matrix3d Krotz = Eigen::Matrix3d::Identity() * max_diag / 1000;

	//Krotz << 1, -0.5, -0.5,
	//	-0.5, 1, -0.5,
	//	-0.5, -0.5, 1;
	//Krotz *=  0.03 * Mat->Young * thickness * Area() / 1000;

	//Krotz.fill(100);
	//std::cout << "K rotz: \n" << Krotz << std::endl;
	//std::cout << "-------------------" << std::endl;
	int indices_rotz[3]{ 5,11,17 };
	int indices_pln[6]{ 0,1,6,7,12,13 };
	int indices_plt[9]{ 2,3,4,8,9,10,14,15,16 };

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

	//std::cout << "K: \n" << mat << std::endl;
	//std::cout << "-------------------" << std::endl;
	//std::cout << "rows: " << mat.rows() << "cols: " << mat.cols() << std::endl;
	Eigen::MatrixXd trMat = trans_matrix();
	return trMat.transpose() * mat * trMat;
}

Eigen::MatrixXd TriPlateElement::GeometricStiffnessMatrix(const std::vector<Displacement>& disp)
{
	Eigen::MatrixXd KGpln = plane_element.geometric_local_stiffness_matrix(disp);
	Eigen::MatrixXd KGplt = geometric_local_stiffness_matrix(disp);

	//int indices_rotz[3]{ 5,11,17 };
	int indices_pln[9]{ 0,1,2,6,7,8,12,13,14 };
	//int indices_pln[9]{ 0,1,-1,6,7,-1,12,13,-1 };
	int indices_plt[9]{ 2,3,4,8,9,10,14,15,16 };

	Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(total_dof, total_dof);

	for (size_t i = 0; i < 9; i++)
	{
		int ir = indices_pln[i];
		if (ir < 0) continue; // -1は無視
		for (size_t j = 0; j < 9; j++)
		{
			int ic = indices_pln[j];
			if (ic < 0) continue; // -1は無視
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

void TriPlateElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat)
{
	Eigen::VectorXd mass = NodeLumpedMass();
	for (size_t ni = 0; ni < node_num; ni++)
	{
		for (size_t i = 0; i < 3; i++) {
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

PlateStressData	TriPlateElement::stress(
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
	d2ndx2 << 4 * dLdx(0) * dLdx(0),  // d²N₁/dx²
		4 * dLdx(1) * dLdx(1),  // d²N₂/dx²
		4 * dLdx(2) * dLdx(2),  // d²N₃/dx²
		8 * dLdx(2) * dLdx(1),  // d²N₄/dx²
		8 * dLdx(0) * dLdx(2),  // d²N₅/dx²
		8 * dLdx(1) * dLdx(0);  // d²N₆/dx²

	// y方向の二階導関数
	Eigen::VectorXd d2ndy2(6);
	d2ndy2 << 4 * dLdy(0) * dLdy(0),  // d²N₁/dy²
		4 * dLdy(1) * dLdy(1),  // d²N₂/dy²
		4 * dLdy(2) * dLdy(2),  // d²N₃/dy²
		8 * dLdy(1) * dLdy(2),  // d²N₄/dy²
		8 * dLdy(0) * dLdy(2),  // d²N₅/dy²
		8 * dLdy(1) * dLdy(0);  // d²N₆/dy²

	// x,y方向の二階導関数
	Eigen::VectorXd d2ndxdy(6);
	d2ndxdy << 4 * dLdx(0) * dLdy(0),  // d²N₁/dxdy
		4 * dLdx(1) * dLdy(1),  // d²N₂/dxdy
		4 * dLdx(2) * dLdy(2),  // d²N₃/dxdy
		4 * (dLdx(1) * dLdy(2) + dLdy(1) * dLdx(2)),  // d²N₄/dxdy
		4 * (dLdx(2) * dLdy(0) + dLdy(2) * dLdx(0)),  // d²N₅/dxdy
		4 * (dLdx(0) * dLdy(1) + dLdy(0) * dLdx(1));  // d²N₆/dxdy

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


//PlateStressData TriPlateElement::stress_save(
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
//}


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
		dndxi(0)* p1.y + dndxi(1) * p2.y + dndxi(2) * p3.y + dndxi(3) * p4.y,
		dndeta(0)* p1.x + dndeta(1) * p2.x + dndeta(2) * p3.x + dndeta(3) * p4.x,
		dndeta(0)* p1.y + dndeta(1) * p2.y + dndeta(2) * p3.y + dndeta(3) * p4.y;

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
		//JMat = JMat.inverse().eval();

		Eigen::MatrixXd dndx(2, 4);
		dndx.row(0) = dndxi;
		dndx.row(1) = dndeta;
		dndx = JMat.inverse() * dndx;
		//dndx = JMat * dndx;

		Eigen::MatrixXd Gmat = Eigen::MatrixXd::Zero(6, 12);
		Gmat << 
			dndx(0, 0), 0, 0, dndx(0, 1), 0, 0, dndx(0, 2), 0, 0, dndx(0, 3), 0, 0,
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
		//Kg += thickness.plane_thick * intgp * Gmat.transpose() * SigMat3 * Gmat * detJ;
	}

	//Eigen::MatrixXd tr = trans_matrix();
	//return tr.transpose() * Kg * tr;
	return Kg;
}

QuadPlaneElement::QuadPlaneElement(Node* n0, Node* n1, Node* n2, Node* n3, double t, Material mat)
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

QuadPlaneElement::QuadPlaneElement(Node* n0, Node* n1, Node* n2, Node* n3, Thickness t, Material mat)
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
	matrix.setZero();  // 行列を0で初期化

	for (int i = 0; i < numBlocks; ++i)
		matrix.block(i * 2, i * node_dof, 2, node_dof) = tr0.block(0, 0, 2, node_dof);

	return matrix;
}

Eigen::MatrixXd QuadPlaneElement::NodeConsistentMass()
{
	return Eigen::MatrixXd::Identity(total_dof, total_dof);
}

//std::vector<NodeLoadData> QuadPlaneElement::BodyforceToNodeLoadData(Eigen::Vector3d accel_vec)
//{
//	return std::vector<NodeLoadData>();
//}

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
	//const Eigen::Vector2d intg_weights(1, 1);
	//const Eigen::Vector2d intg_params(-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0));

	// 3-Point Gauss Quadrature
	const Eigen::Vector3d intg_weights(5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0);
	const Eigen::Vector3d intg_params(-sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0));

	// 4-Point Gauss Quadrature
	//const Eigen::Vector4d intg_weights(0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451);
	//const Eigen::Vector4d intg_params(-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116);

	// 5-Point Gauss Quadrature
	//Eigen::VectorXd intg_weights(5), intg_params(5);
	//intg_weights << 0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851;
	//intg_params << -0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459;

	double t3 = thickness.plane_thick * thickness.plane_thick * thickness.plane_thick;

	Eigen::MatrixXd M(total_dof, total_dof);
	M.setZero();

	Eigen::Vector<double, 12> p_vec;
	p_vec << load_vecs[0].x, load_vecs[0].y, load_vecs[0].z,
		load_vecs[1].x, load_vecs[1].y, load_vecs[1].z,
		load_vecs[2].x, load_vecs[2].y, load_vecs[2].z,
		load_vecs[3].x, load_vecs[3].y, load_vecs[3].z;
	//Eigen::Vector<int, 12> shape_indices;
	//const int shape_indices[12] = { 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20 };
	// shape_indices << 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20;

	//Eigen::Vector<int, 12> HVecs_indices;
	//HVecs_indices << 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22;
	//const int HVecs_indices[12] = { 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22 };

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

void QuadPlaneElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat)
{
	Eigen::VectorXd mass = NodeLumpedMass();
	for (size_t ni = 0; ni < node_num; ni++)
	{
		for (size_t i = 0; i < 3; i++) {
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
	//return MembraneStressData();
}

/// <summary>
/// Coumpute the following Jacobi matrix
/// +                     +
/// | dx / dxi   dy / dxi  |
/// | dx / deta  dy / deta |
/// +                     +
/// </summary>
/// <param name="xi">parameter xi</param>
/// <param name="eta">parameter eta</param>
/// <returns></returns>
Eigen::Matrix2d QuadPlateElement::JMatrix(double xi, double eta)
{
	Point p1 = plane.PointToCoord(Nodes[0]->Location);
	Point p2 = plane.PointToCoord(Nodes[1]->Location);
	Point p3 = plane.PointToCoord(Nodes[2]->Location);
	Point p4 = plane.PointToCoord(Nodes[3]->Location);

	Vector v12 = p1 - p2;
	Vector v23 = p2 - p3;
	Vector v34 = p3 - p4;
	Vector v41 = p4 - p1;
	
	Eigen::Matrix2d JMat;
	JMat <<	-v12.x + v34.x + eta * (v12.x + v34.x), -v12.y + v34.y + eta * (v12.y + v34.y),
		-v23.x + v41.x + xi * (v12.x + v34.x), -v23.y + v41.y + xi * (v12.y + v34.y);
	return JMat * 0.25; 
}

Eigen::Matrix2d QuadPlateElement::dJinv_dxi(double xi, double eta)
{
	Point p1 = plane.PointToCoord(Nodes[0]->Location);
	Point p2 = plane.PointToCoord(Nodes[1]->Location);
	Point p3 = plane.PointToCoord(Nodes[2]->Location);
	Point p4 = plane.PointToCoord(Nodes[3]->Location);

	Vector v12 = p1 - p2;
	Vector v23 = p2 - p3;
	Vector v34 = p3 - p4;
	Vector v41 = p4 - p1;

	Eigen::Matrix2d JMat;
	JMat << -v12.x + v34.x + eta * (v12.x + v34.x), -v12.y + v34.y + eta * (v12.y + v34.y),
		-v23.x + v41.x + xi * (v12.x + v34.x), -v23.y + v41.y + xi * (v12.y + v34.y);
	JMat *= 0.25;

	double detJ = JMat.determinant();
	double ddetJ_dxi = 0.125 * (-v34.y * v12.x + v12.y * v34.x);
	Eigen::Matrix2d Matdfdxi;
	
	Matdfdxi << v12.y + v34.y, -0,
		-(v12.x + v34.x), 0;
	Matdfdxi *= 0.25;

	Eigen::Matrix2d dJinv_dxi;
	dJinv_dxi = 1 / detJ * Matdfdxi - ddetJ_dxi / detJ * JMat.inverse();
	
	return dJinv_dxi;
}

//Eigen::Matrix2d QuadPlateElement::dJinv_dxi(double xi, double eta)
//{
//	Point p1 = plane.PointToCoord(Nodes[0]->Location);
//	Point p2 = plane.PointToCoord(Nodes[1]->Location);
//	Point p3 = plane.PointToCoord(Nodes[2]->Location);
//	Point p4 = plane.PointToCoord(Nodes[3]->Location);
//
//	Vector v12 = p1 - p2;
//	Vector v23 = p2 - p3;
//	Vector v34 = p3 - p4;
//	Vector v41 = p4 - p1;
//
//	Eigen::Matrix2d JMat;
//	JMat << -v12.x + v34.x + eta * (v12.x + v34.x), -v12.y + v34.y + eta * (v12.y + v34.y),
//		-v23.x + v41.x + xi * (v12.x + v34.x), -v23.y + v41.y + xi * (v12.y + v34.y);
//	JMat *= 0.25;
//
//	double detJ = JMat.determinant();
//	double ddetJ_dxi = 0.125 * (-v34.y * v12.x + v12.y * v34.x);
//	Eigen::Matrix2d Matdfdxi;
//
//	Matdfdxi << v12.y + v34.y, -0,
//		-(v12.x + v34.x), 0;
//
//	Eigen::Matrix2d dJinv_dxi;
//	dJinv_dxi = 0.25 / detJ * Matdfdxi - ddetJ_dxi / detJ * JMat.inverse();
//
//	return dJinv_dxi;
//}

Eigen::Matrix2d QuadPlateElement::dJinv_deta(double xi, double eta)
{
	Point p1 = plane.PointToCoord(Nodes[0]->Location);
	Point p2 = plane.PointToCoord(Nodes[1]->Location);
	Point p3 = plane.PointToCoord(Nodes[2]->Location);
	Point p4 = plane.PointToCoord(Nodes[3]->Location);

	Vector v12 = p1 - p2;
	Vector v23 = p2 - p3;
	Vector v34 = p3 - p4;
	Vector v41 = p4 - p1;

	Eigen::Matrix2d JMat;
	JMat << -v12.x + v34.x + eta * (v12.x + v34.x), -v12.y + v34.y + eta * (v12.y + v34.y),
		-v23.x + v41.x + xi * (v12.x + v34.x), -v23.y + v41.y + xi * (v12.y + v34.y);
	JMat *= 0.25;

	double detJ = JMat.determinant();
	double ddetJ_deta = 0.125 * (-v41.y * v23.x + v23.y * v41.x);
	Eigen::Matrix2d Matdfdeta;

	Matdfdeta << 0, -(v12.y + v34.y),
		-0, v12.x + v34.x;
	Matdfdeta *= 0.25;

	Eigen::Matrix2d dJinv_deta;
	dJinv_deta = 1.0 / detJ * Matdfdeta - ddetJ_deta / detJ * JMat.inverse();

	return dJinv_deta;
}

// 入力形状関数に応じたx,yそれぞれの方向の列ベクトルによるHベクトル行列
Eigen::MatrixXd QuadPlateElement::HVecs(Eigen::VectorXd shape_funcs) {
	Point p1 = plane.PointToCoord(Nodes[0]->Location);
	Point p2 = plane.PointToCoord(Nodes[1]->Location);
	Point p3 = plane.PointToCoord(Nodes[2]->Location);
	Point p4 = plane.PointToCoord(Nodes[3]->Location);

	Vector v12 = p1 - p2;
	Vector v23 = p2 - p3;
	Vector v34 = p3 - p4;
	Vector v41 = p4 - p1;

	double l12 = v12.norm();
	double l23 = v23.norm();
	double l34 = v34.norm();
	double l41 = v41.norm();

	double c5 = (p2.y - p1.y) / l12;
	double c6 = (p3.y - p2.y) / l23;
	double c7 = (p4.y - p3.y) / l34;
	double c8 = (p1.y - p4.y) / l41;

	double s5 = (p1.x - p2.x) / l12;
	double s6 = (p2.x - p3.x) / l23;
	double s7 = (p3.x - p4.x) / l34;
	double s8 = (p4.x - p1.x) / l41;

	Eigen::Matrix<double, 2, total_dof_local> HVecs;
	
	// Hx
	HVecs.row(0) << 1.5 * (s8 / l41 * shape_funcs(7) - s5 / l12 * shape_funcs(4)),
		-0.75 * (s5 * c5 * shape_funcs(4) + s8 * c8 * shape_funcs(7)),
		shape_funcs(0) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * shape_funcs(4) + (c8 * c8 * 0.5 - s8 * s8 * 0.25) * shape_funcs(7),
		1.5 * (s5 / l12 * shape_funcs(4) - s6 / l23 * shape_funcs(5)),
		-0.75 * (s5 * c5 * shape_funcs(4) + s6 * c6 * shape_funcs(5)),
		shape_funcs(1) + (c5 * c5 * 0.5 - s5 * s5 * 0.25) * shape_funcs(4) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * shape_funcs(5),
		1.5 * (s6 / l23 * shape_funcs(5) - s7 / l34 * shape_funcs(6)),
		-0.75 * (s6 * c6 * shape_funcs(5) + s7 * c7 * shape_funcs(6)),
		shape_funcs(2) + (c6 * c6 * 0.5 - s6 * s6 * 0.25) * shape_funcs(5) + (c7 * c7 * 0.5 - s7 * s7 * 0.25) * shape_funcs(6),
		1.5 * (s7 / l34 * shape_funcs(6) - s8 / l41 * shape_funcs(7)),
		-0.75 * (s7 * c7 * shape_funcs(6) + s8 * c8 * shape_funcs(7)),
		shape_funcs(3) + (c7 * c7 * 0.5 - s7 * s7 * 0.25) * shape_funcs(6) + (c8 * c8 * 0.5 - s8 * s8 * 0.25) * shape_funcs(7);
	// Hy
	HVecs.row(1) << 1.5 * (c5 / l12 * shape_funcs(4) - c8 / l41 * shape_funcs(7)),
		-shape_funcs(0) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * shape_funcs(4) + (c8 * c8 * 0.25 - s8 * s8 * 0.5) * shape_funcs(7),
		0.75 * (s5 * c5 * shape_funcs(4) + s8 * c8 * shape_funcs(7)),
		1.5 * (c6 / l23 * shape_funcs(5) - c5 / l12 * shape_funcs(4)),
		-shape_funcs(1) + (c5 * c5 * 0.25 - s5 * s5 * 0.5) * shape_funcs(4) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * shape_funcs(5),
		0.75 * (s5 * c5 * shape_funcs(4) + s6 * c6 * shape_funcs(5)),
		1.5 * (c7 / l34 * shape_funcs(6) - c6 / l23 * shape_funcs(5)),
		-shape_funcs(2) + (c6 * c6 * 0.25 - s6 * s6 * 0.5) * shape_funcs(5) + (c7 * c7 * 0.25 - s7 * s7 * 0.5) * shape_funcs(6),
		0.75 * (s6 * c6 * shape_funcs(5) + s7 * c7 * shape_funcs(6)),
		1.5 * (c8 / l41 * shape_funcs(7) - c7 / l34 * shape_funcs(6)),
		-shape_funcs(3) + (c7 * c7 * 0.25 - s7 * s7 * 0.5) * shape_funcs(6) + (c8 * c8 * 0.25 - s8 * s8 * 0.5) * shape_funcs(7),
		0.75 * (s7 * c7 * shape_funcs(6) + s8 * c8 * shape_funcs(7));

	return HVecs;
}

//Eigen::MatrixXd QuadPlateElement::HVecs(Eigen::VectorXd shape_funcs) {
//	Point p1 = plane.PointToCoord(Nodes[0]->Location);
//	Point p2 = plane.PointToCoord(Nodes[1]->Location);
//	Point p3 = plane.PointToCoord(Nodes[2]->Location);
//	Point p4 = plane.PointToCoord(Nodes[3]->Location);
//
//	Vector v12 = p1 - p2;
//	Vector v23 = p2 - p3;
//	Vector v34 = p3 - p4;
//	Vector v41 = p4 - p1;
//
//	double a5 = -v12.x / v12.squared_norm();
//	double a6 = -v23.x / v23.squared_norm();
//	double a7 = -v34.x / v34.squared_norm();
//	double a8 = -v41.x / v41.squared_norm();
//
//	double b5 = 0.75 * v12.x * v12.y / v12.squared_norm();
//	double b6 = 0.75 * v23.x * v23.y / v23.squared_norm();
//	double b7 = 0.75 * v34.x * v34.y / v34.squared_norm();
//	double b8 = 0.75 * v41.x * v41.y / v41.squared_norm();
//
//	double c5 = (0.25 * v12.x * v12.x - 0.5 * v12.y * v12.y) / v12.squared_norm();
//	double c6 = (0.25 * v23.x * v23.x - 0.5 * v23.y * v23.y) / v23.squared_norm();
//	double c7 = (0.25 * v34.x * v34.x - 0.5 * v34.y * v34.y) / v34.squared_norm();
//	double c8 = (0.25 * v41.x * v41.x - 0.5 * v41.y * v41.y) / v41.squared_norm();
//
//	double d5 = -v12.y / v12.squared_norm();
//	double d6 = -v23.y / v23.squared_norm();
//	double d7 = -v34.y / v34.squared_norm();
//	double d8 = -v41.y / v41.squared_norm();
//
//	double e5 = (-0.5 * v12.x * v12.x + 0.25 * v12.y * v12.y) / v12.squared_norm();
//	double e6 = (-0.5 * v23.x * v23.x + 0.25 * v23.y * v23.y) / v23.squared_norm();
//	double e7 = (-0.5 * v34.x * v34.x + 0.25 * v34.y * v34.y) / v34.squared_norm();
//	double e8 = (-0.5 * v41.x * v41.x + 0.25 * v41.y * v41.y) / v41.squared_norm();
//
//	// Eigen::Vector<double, 8> shape_funcs;
//
//	Eigen::Matrix<double, 2, total_dof_local> HVecs;
//	// Eigen::Vector<double, total_dof_local> Hx_x, Hx_eta, Hy_xi, Hy_eta;
//
//	HVecs.row(0) << 1.5 * (a5 * shape_funcs(4) - a8 * shape_funcs(7)),
//		b5* shape_funcs(4) + b8 * shape_funcs(7),
//		shape_funcs(0) - c5 * shape_funcs(4) - c8 * shape_funcs(7),
//		1.5 * (a6 * shape_funcs(5) - a5 * shape_funcs(4)),
//		b6* shape_funcs(5) + b5 * shape_funcs(4),
//		shape_funcs(1) - c6 * shape_funcs(5) - c5 * shape_funcs(4),
//		1.5 * (a7 * shape_funcs(6) - a6 * shape_funcs(5)),
//		b7* shape_funcs(6) + b6 * shape_funcs(5),
//		shape_funcs(2) - c7 * shape_funcs(6) - c6 * shape_funcs(5),
//		1.5 * (a8 * shape_funcs(7) - a7 * shape_funcs(6)),
//		b8* shape_funcs(7) + b7 * shape_funcs(6),
//		shape_funcs(3) - c8 * shape_funcs(7) - c7 * shape_funcs(6);
//
//	HVecs.row(1) << 1.5 * (d5 * shape_funcs(4) - d8 * shape_funcs(7)),
//		-shape_funcs(0) + e5 * shape_funcs(4) + e8 * shape_funcs(7),
//		-b5 * shape_funcs(4) - b8 * shape_funcs(7),
//		1.5 * (d6 * shape_funcs(5) - d5 * shape_funcs(4)),
//		-shape_funcs(1) + e6 * shape_funcs(5) + e5 * shape_funcs(4),
//		-b6 * shape_funcs(5) - b5 * shape_funcs(4),
//		1.5 * (d7 * shape_funcs(6) - d6 * shape_funcs(5)),
//		-shape_funcs(2) + e7 * shape_funcs(6) + e6 * shape_funcs(5),
//		-b7 * shape_funcs(6) - b6 * shape_funcs(5),
//		1.5 * (d8 * shape_funcs(7) - d7 * shape_funcs(6)),
//		-shape_funcs(3) + e8 * shape_funcs(7) + e7 * shape_funcs(6),
//		-b8 * shape_funcs(7) - b7 * shape_funcs(6);
//
//	return HVecs;
//}

Eigen::Vector<double, 8> shape_func(double xi, double eta)
{
	Eigen::Vector<double, 8> sf;

	sf(0) = 0.25 * (1 - xi) * (1 - eta) * (-1 - xi - eta); // N1
	sf(1) = 0.25 * (1 + xi) * (1 - eta) * (-1 + xi - eta); // N2
	sf(2) = 0.25 * (1 + xi) * (1 + eta) * (-1 + xi + eta); // N3
	sf(3) = 0.25 * (1 - xi) * (1 + eta) * (-1 - xi + eta); // N4
	sf(4) = 0.5 * (1 - xi * xi) * (1 - eta);			   // N5
	sf(5) = 0.5 * (1 + xi) * (1 - eta * eta);			   // N6
	sf(6) = 0.5 * (1 - xi * xi) * (1 + eta);			   // N7
	sf(7) = 0.5 * (1 - xi) * (1 - eta * eta);			   // N8

	return sf;
}

Eigen::Vector<double, 8> dndxi(double xi, double eta) {
	Eigen::Vector<double, 8> dndxi;
	dndxi << 0.25 * (1 - eta) * (2 * xi + eta),
		0.25 * (1 - eta) * (2 * xi - eta),
		0.25 * (1 + eta) * (2 * xi + eta),
		0.25 * (1 + eta) * (2 * xi - eta),
		-xi * (1 - eta),
		0.5 * (1 - eta * eta),
		-xi * (1 + eta),
		-0.5 * (1 - eta * eta);

	return dndxi;
}

Eigen::Vector<double, 8> dndeta(double xi, double eta) {
	Eigen::Vector<double, 8> dndeta;
	dndeta << 0.25 * (1 - xi) * (2 * eta + xi),
		0.25 * (1 + xi) * (2 * eta - xi),
		0.25 * (1 + xi) * (2 * eta + xi),
		0.25 * (1 - xi) * (2 * eta - xi),
		-0.5 * (1 - xi * xi),
		-eta * (1 + xi),
		0.5 * (1 - xi * xi),
		-eta * (1 - xi);

	return dndeta;
}

Eigen::Vector<double, 8> dn2_dxi2(double xi, double eta) {
	Eigen::Vector<double, 8> dn2_dxi2;
	dn2_dxi2 << 0.5 * (1 - eta),
		0.5 * (1 - eta),
		0.5 * (1 + eta),
		0.5 * (1 + eta),
		eta - 1,
		0,
		-(1.0 + eta),
		0;

	return dn2_dxi2;
}

Eigen::Vector<double, 8> dn2_deta2(double xi, double eta) {
	Eigen::Vector<double, 8> dn2_deta2;
	dn2_deta2 << 0.5 * (1 - xi),
		0.5 * (1 + xi),
		0.5 * (1 + xi),
		0.5 * (1 - xi),
		0,
		-(1.0 + xi),
		0,
		xi - 1;

	return dn2_deta2;
}

Eigen::Vector<double, 8> dn2_dxieta(double xi, double eta) {
	Eigen::Vector<double, 8> dn2_dxieta;
	dn2_dxieta << -0.25 * (2 * xi + 2 * eta - 1),
		0.25 * (-2 * xi + 2 * eta - 1),
		0.25 * (2 * xi + 2 * eta + 1),
		0.25 * (2 * xi - 2 * eta - 1),
		xi,
		-eta,
		-xi,
		eta;

	return dn2_dxieta;
}

Eigen::MatrixXd QuadPlateElement::BMatrix(double xi, double eta)
{
	Eigen::MatrixXd H_xi = HVecs(dndxi(xi, eta));
	Eigen::MatrixXd H_eta = HVecs(dndeta(xi, eta));

	Eigen::Vector<double, total_dof_local> Hx_xi, Hx_eta, Hy_xi, Hy_eta;
	Hx_xi = H_xi.row(0);
	Hx_eta = H_eta.row(0);
	Hy_xi = H_xi.row(1);
	Hy_eta = H_eta.row(1);

	Eigen::Matrix2d Jinv = JMatrix(xi, eta).inverse();
	Eigen::MatrixXd BMat(3, 12);
	BMat.row(0) = Jinv(0, 0) * Hx_xi + Jinv(0, 1) * Hx_eta;
	BMat.row(1) = Jinv(1, 0) * Hy_xi + Jinv(1, 1) * Hy_eta;
	BMat.row(2) = Jinv(0, 0) * Hy_xi + Jinv(0, 1) * Hy_eta + Jinv(1, 0) * Hx_xi + Jinv(1, 1) * Hx_eta;
	
	return BMat;
}

/*
Eigen::MatrixXd QuadPlateElement::BMatrix(double xi, double eta)
{
	Point p1 = plane.PointToCoord(Nodes[0]->Location);
	Point p2 = plane.PointToCoord(Nodes[1]->Location);
	Point p3 = plane.PointToCoord(Nodes[2]->Location);
	Point p4 = plane.PointToCoord(Nodes[3]->Location);

	Vector v12 = p1 - p2;
	Vector v23 = p2 - p3;
	Vector v34 = p3 - p4;
	Vector v41 = p4 - p1;

	double a5 = -v12.x / v12.squared_norm();
	double a6 = -v23.x / v23.squared_norm();
	double a7 = -v34.x / v34.squared_norm();
	double a8 = -v41.x / v41.squared_norm();

	double b5 = 0.75 * v12.x * v12.y / v12.squared_norm();
	double b6 = 0.75 * v23.x * v23.y / v23.squared_norm();
	double b7 = 0.75 * v34.x * v34.y / v34.squared_norm();
	double b8 = 0.75 * v41.x * v41.y / v41.squared_norm();

	double c5 = (0.25 * v12.x * v12.x - 0.5 * v12.y * v12.y) / v12.squared_norm();
	double c6 = (0.25 * v23.x * v23.x - 0.5 * v23.y * v23.y) / v23.squared_norm();
	double c7 = (0.25 * v34.x * v34.x - 0.5 * v34.y * v34.y) / v34.squared_norm();
	double c8 = (0.25 * v41.x * v41.x - 0.5 * v41.y * v41.y) / v41.squared_norm();

	double d5 = -v12.y / v12.squared_norm();
	double d6 = -v23.y / v23.squared_norm();
	double d7 = -v34.y / v34.squared_norm();
	double d8 = -v41.y / v41.squared_norm();

	double e5 = (-0.5 * v12.x * v12.x + 0.25 * v12.y * v12.y) / v12.squared_norm();
	double e6 = (-0.5 * v23.x * v23.x + 0.25 * v23.y * v23.y) / v23.squared_norm();
	double e7 = (-0.5 * v34.x * v34.x + 0.25 * v34.y * v34.y) / v34.squared_norm();
	double e8 = (-0.5 * v41.x * v41.x + 0.25 * v41.y * v41.y) / v41.squared_norm();

	Eigen::Vector<double, 8> dndxi, dndeta;
	dndxi << 0.25 * (1 - eta) * (2 * xi + eta),
		0.25 * (1 - eta) * (2 * xi - eta),
		0.25 * (1 + eta) * (2 * xi + eta),
		0.25 * (1 + eta) * (2 * xi - eta),
		-xi * (1 - eta),
		0.5 * (1 - eta * eta),
		-xi * (1 + eta),
		-0.5 * (1 - eta * eta);

	dndeta << 0.25 * (1 - xi) * (2 * eta + xi),
		0.25 * (1 + xi) * (2 * eta - xi),
		0.25 * (1 + xi) * (2 * eta + xi),
		0.25 * (1 - xi) * (2 * eta - xi),
		-0.5 * (1 - xi * xi),
		-eta * (1 + xi),
		0.5 * (1 - xi * xi),
		-eta * (1 - xi);

	Eigen::Vector<double, total_dof_local> Hx_xi, Hx_eta, Hy_xi, Hy_eta;
	Hx_xi << 1.5 * (a5 * dndxi(4) - a8 * dndxi(7)),
		b5* dndxi(4) + b8 * dndxi(7),
		dndxi(0) - c5 * dndxi(4) - c8 * dndxi(7),
		1.5 * (a6 * dndxi(5) - a5 * dndxi(4)),
		b6* dndxi(5) + b5 * dndxi(4),
		dndxi(1) - c6 * dndxi(5) - c5 * dndxi(4),
		1.5 * (a7 * dndxi(6) - a6 * dndxi(5)),
		b7* dndxi(6) + b6 * dndxi(5),
		dndxi(2) - c7 * dndxi(6) - c6 * dndxi(5),
		1.5 * (a8 * dndxi(7) - a7 * dndxi(6)),
		b8* dndxi(7) + b7 * dndxi(6),
		dndxi(3) - c8 * dndxi(7) - c7 * dndxi(6);

	Hx_eta << 1.5 * (a5 * dndeta(4) - a8 * dndeta(7)),
		b5* dndeta(4) + b8 * dndeta(7),
		dndeta(0) - c5 * dndeta(4) - c8 * dndeta(7),
		1.5 * (a6 * dndeta(5) - a5 * dndeta(4)),
		b6* dndeta(5) + b5 * dndeta(4),
		dndeta(1) - c6 * dndeta(5) - c5 * dndeta(4),
		1.5 * (a7 * dndeta(6) - a6 * dndeta(5)),
		b7* dndeta(6) + b6 * dndeta(5),
		dndeta(2) - c7 * dndeta(6) - c6 * dndeta(5),
		1.5 * (a8 * dndeta(7) - a7 * dndeta(6)),
		b8* dndeta(7) + b7 * dndeta(6),
		dndeta(3) - c8 * dndeta(7) - c7 * dndeta(6);

	Hy_xi << 1.5 * (d5 * dndxi(4) - d8 * dndxi(7)),
		-dndxi(0) + e5 * dndxi(4) + e8 * dndxi(7),
		-b5 * dndxi(4) - b8 * dndxi(7),
		1.5 * (d6 * dndxi(5) - d5 * dndxi(4)),
		-dndxi(1) + e6 * dndxi(5) + e5 * dndxi(4),
		-b6 * dndxi(5) - b5 * dndxi(4),
		1.5 * (d7 * dndxi(6) - d6 * dndxi(5)),
		-dndxi(2) + e7 * dndxi(6) + e6 * dndxi(5),
		-b7 * dndxi(6) - b6 * dndxi(5),
		1.5 * (d8 * dndxi(7) - d7 * dndxi(6)),
		-dndxi(3) + e8 * dndxi(7) + e7 * dndxi(6),
		-b8 * dndxi(7) - b7 * dndxi(6);

	Hy_eta << 1.5 * (d5 * dndeta(4) - d8 * dndeta(7)),
		-dndeta(0) + e5 * dndeta(4) + e8 * dndeta(7),
		-b5 * dndeta(4) - b8 * dndeta(7),
		1.5 * (d6 * dndeta(5) - d5 * dndeta(4)),
		-dndeta(1) + e6 * dndeta(5) + e5 * dndeta(4),
		-b6 * dndeta(5) - b5 * dndeta(4),
		1.5 * (d7 * dndeta(6) - d6 * dndeta(5)),
		-dndeta(2) + e7 * dndeta(6) + e6 * dndeta(5),
		-b7 * dndeta(6) - b6 * dndeta(5),
		1.5 * (d8 * dndeta(7) - d7 * dndeta(6)),
		-dndeta(3) + e8 * dndeta(7) + e7 * dndeta(6),
		-b8 * dndeta(7) - b7 * dndeta(6);

	Eigen::Matrix2d Jinv = JMatrix(xi, eta).inverse();
	Eigen::MatrixXd BMat(3, 12);
	BMat.row(0) = Jinv(0, 0) * Hx_xi + Jinv(0, 1) * Hx_eta;
	BMat.row(1) = Jinv(1, 0) * Hy_xi + Jinv(1, 1) * Hy_eta;
	BMat.row(2) = Jinv(0, 0) * Hy_xi + Jinv(0, 1) * Hy_eta + Jinv(1, 0) * Hx_xi + Jinv(1, 1) * Hx_eta;

	return BMat;
}
*/

Eigen::Matrix3d QuadPlateElement::DMatrix()
{
	Eigen::Matrix3d mat;
	mat << 1, Mat.Poisson, 0,
		Mat.Poisson, 1, 0,
		0, 0, (1.0 - Mat.Poisson) / 2.0;
	return Mat.Young * pow(thickness.plate_thick, 3) / 12 / (1 - Mat.Poisson * Mat.Poisson) * mat;
}

Eigen::MatrixXd QuadPlateElement::localStiffnessMatrix()
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
		K += BMat.transpose() * DMat * BMat * detJ;
	}
	return K;
}

Eigen::MatrixXd QuadPlateElement::geometric_local_stiffness_matrix(const std::vector<Displacement> &disp)
{
	const double intgp = 1.0 / sqrt(3.0);
	Eigen::Vector4d xi_list(-intgp, intgp, intgp, -intgp);
	Eigen::Vector4d eta_list(-intgp, -intgp, intgp, intgp);

	Eigen::MatrixXd Kg = Eigen::MatrixXd::Zero(12, 12);
	

	for (size_t i = 0; i < node_num; i++)
	{
		double xi = xi_list(i);
		double eta = eta_list(i);

		Eigen::MatrixXd hvecs = HVecs(shape_func(xi, eta));

		PlateStressData strs = stress(disp[0], disp[1], disp[2], disp[3], xi, eta);
		Eigen::Matrix2d SigMat;
		SigMat << strs.My, strs.Mxy, strs.Mxy, strs.Mx;

		double detJ = JMatrix(xi, eta).determinant();
		Kg += hvecs.transpose() * SigMat * hvecs * detJ;
		//Kg += hvecs.transpose() * SigMat * hvecs * intgp * detJ;
	}

	//Eigen::MatrixXd tr = trans_matrix();
	//return tr.transpose() * Kg * tr;
	return Kg;
}

QuadPlateElement::QuadPlateElement(Node* n0, Node* n1, Node* n2, Node* n3, double t, Material mat)
	: plane_element(n0, n1, n2, n3, t, mat)
{
	Nodes[0] = n0;
	Nodes[1] = n1;
	Nodes[2] = n2;
	Nodes[3] = n3;

	//thickness = t;
	thickness = Thickness(t);
	Mat = mat;

	Point p12 = (n1->Location + n2->Location) / 2;
	Point p23 = (n2->Location + n3->Location) / 2;
	Point p30 = (n3->Location + n0->Location) / 2;
	plane = Plane::CreateFromPoints(p30, p12, p23);
	//plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
}

QuadPlateElement::QuadPlateElement(Node* n0, Node* n1, Node* n2, Node* n3, Thickness t, Material mat)
	: plane_element(n0, n1, n2, n3, t, mat)
{
	Nodes[0] = n0;
	Nodes[1] = n1;
	Nodes[2] = n2;
	Nodes[3] = n3;

	thickness = t;
	Mat = mat;

	Point p12 = (n1->Location + n2->Location) / 2;
	Point p23 = (n2->Location + n3->Location) / 2;
	Point p30 = (n3->Location + n0->Location) / 2;
	plane = Plane::CreateFromPoints(p30, p12, p23);
}

QuadPlateElement::LocalMatrixd
QuadPlateElement::trans_matrix()
{
	Eigen::Matrix3d tr0 = trans_matrix3(plane);

	// ブロックの数を指定
	const int numBlocks = 8;

	// 大きな行列を作成し、対角ブロック行列として同じ行列を配置
	LocalMatrixd matrix;
	matrix.setZero();

	for (int i = 0; i < numBlocks; ++i)
		matrix.block(i * tr0.rows(), i * tr0.cols(), tr0.rows(), tr0.cols()) = tr0;

	// Debug
	// 250412 1bと2はテストで有効性が確認できたが、実際に導入すると微妙なので
	// もう少し確認が必要。特にメソッド2がかなり変な挙動をする。
	//matrix = WarpCorrectMatrix1a().transpose() * matrix;  // 不採用
	matrix = WarpCorrectMatrix1b().transpose() * matrix;
	matrix = WarpCorrectMatrix2().transpose() * matrix;

	return matrix;
}

Eigen::MatrixXd QuadPlateElement::WarpCorrectMatrix1a()
{
	double h1 = plane.DistanceTo(Nodes[0]->Location);
	Eigen::Vector4d sins, coss, lens;
	Eigen::Vector4d fx_fixs, fy_fixs;
	for (int i = 0; i < 4; i++)
	{
		double li = Nodes[i]->Location.distance_to(Nodes[(i + 1) % 4]->Location);
		lens(i) = li;
		sins(i) = (Nodes[(i + 1) % 4]->Location.y - Nodes[i]->Location.y) / li;
		coss(i) = (Nodes[(i + 1) % 4]->Location.x - Nodes[i]->Location.x) / li;
	}

	for (int i = 0; i < 4; i++)
	{
		double hi = h1 * pow(-1, i);
		double dsin = sins(i) * coss((i + 3) % 4) - coss(i) * sins((i + 3) % 4);
		
		double s = sins(i) / lens((i + 3) % 4) + sins((i + 3) % 4) / lens(i);
		double c = coss(i) / lens((i + 3) % 4) + coss((i + 3) % 4) / lens(i);
		fx_fixs(i) = s * hi / dsin;
		fy_fixs(i) = c * hi / dsin;
	}

	Eigen::MatrixXd mat = Eigen::MatrixXd::Identity(24, 24);
	for (int i = 0; i < 4; i++)
	{
		// 実際動かすと正負逆のような気がする。
		//mat(6 * i + 2, 6 * i) = fx_fixs(i);
		mat(6 * i + 2, 6 * i) = -fx_fixs(i);
		//mat(6 * i + 2, 6 * i + 1) = fy_fixs(i);
		mat(6 * i + 2, 6 * i + 1) = -fy_fixs(i);
	}

	return mat;
}

Eigen::MatrixXd QuadPlateElement::WarpCorrectMatrix1b()
{
	Point p0 = plane.PointToCoord(Nodes[0]->Location);
	Point p1 = plane.PointToCoord(Nodes[1]->Location);
	Point p2 = plane.PointToCoord(Nodes[2]->Location);
	Point p3 = plane.PointToCoord(Nodes[3]->Location);

	Eigen::Matrix2d matAinv, matA;
	matA << p2.x - p0.x, p3.x - p1.x,
		p2.y - p0.y, p3.y - p1.y;
	//matAinv << matA(1, 1), -matA(0, 1),
	//	-matA(1, 0), matA(0, 0);
	//matAinv /= (matA(0, 0) * matA(1, 1) - matA(0, 1) * matA(1, 0));

	 matAinv = matA.inverse();

	Eigen::MatrixXd mat = Eigen::MatrixXd::Identity(24, 24);
	mat(2, 0) = -matAinv(0, 0) * p0.z; mat(2, 1) = -matAinv(0, 1) * p0.z;
	mat(2, 6) = -matAinv(0, 0) * p1.z; mat(2, 7) = -matAinv(0, 1) * p1.z;
	mat(2, 12) = -matAinv(0, 0) * p2.z; mat(2, 13) = -matAinv(0, 1) * p2.z;
	mat(2, 18) = -matAinv(0, 0) * p3.z; mat(2, 19) = -matAinv(0, 1) * p3.z;

	mat(8, 0) = -matAinv(1, 0) * p0.z; mat(8, 1) = -matAinv(1, 1) * p0.z;
	mat(8, 6) = -matAinv(1, 0) * p1.z; mat(8, 7) = -matAinv(1, 1) * p1.z;
	mat(8, 12) = -matAinv(1, 0) * p2.z; mat(8, 13) = -matAinv(1, 1) * p2.z;
	mat(8, 18) = -matAinv(1, 0) * p3.z; mat(8, 19) = -matAinv(1, 1) * p3.z;

	mat(14, 0) = matAinv(0, 0) * p0.z; mat(14, 1) = matAinv(0, 1) * p0.z;
	mat(14, 6) = matAinv(0, 0) * p1.z; mat(14, 7) = matAinv(0, 1) * p1.z;
	mat(14, 12) = matAinv(0, 0) * p2.z; mat(14, 13) = matAinv(0, 1) * p2.z;
	mat(14, 18) = matAinv(0, 0) * p3.z; mat(14, 19) = matAinv(0, 1) * p3.z;

	mat(20, 0) = matAinv(1, 0) * p0.z; mat(20, 1) = matAinv(1, 1) * p0.z;
	mat(20, 6) = matAinv(1, 0) * p1.z; mat(20, 7) = matAinv(1, 1) * p1.z;
	mat(20, 12) = matAinv(1, 0) * p2.z; mat(20, 13) = matAinv(1, 1) * p2.z;
	mat(20, 18) = matAinv(1, 0) * p3.z; mat(20, 19) = matAinv(1, 1) * p3.z;

	return mat;
}

Eigen::MatrixXd QuadPlateElement::WarpCorrectMatrix2()
{
	Point p0 = plane.PointToCoord(Nodes[0]->Location);
	Point p1 = plane.PointToCoord(Nodes[1]->Location);
	Point p2 = plane.PointToCoord(Nodes[2]->Location);
	Point p3 = plane.PointToCoord(Nodes[3]->Location);
	Vector v01 = p1 - p0;
	Vector v12 = p2 - p1;
	Vector v23 = p3 - p2;
	Vector v30 = p0 - p3;

	Vector n0 = Vector::cross(v30, v01); n0 = n0 * (1 / n0.norm());
	Vector n1 = Vector::cross(v01, v12); n1 = n1 * (1 / n1.norm());
	Vector n2 = Vector::cross(v12, v23); n2 = n2 * (1 / n2.norm());
	Vector n3 = Vector::cross(v23, v30); n3 = n3 * (1 / n3.norm());

	double delta = p0.y + p0.x + p1.y - p1.x - p2.y - p2.x - p3.y + p3.x;

	Eigen::VectorXd fas(24);
	fas << 
		p0.y, -p0.x, 0, -n0.x / n0.z, -n0.y / n0.z, 0,
		p1.y, -p1.x, 0, -n1.x / n1.z, -n1.y / n1.z, 0,
		p2.y, -p2.x, 0, -n2.x / n2.z, -n2.y / n2.z, 0,
		p3.y, -p3.x, 0, -n3.x / n3.z, -n3.y / n3.z, 0;
	fas /= delta;

	Eigen::MatrixXd mat = Eigen::MatrixXd::Identity(24, 24);
	mat.row(0) += fas;
	mat.row(1) -= fas;
	mat.row(6) += fas;
	mat.row(7) += fas;
	mat.row(12) -= fas;
	mat.row(13) += fas;
	mat.row(18) -= fas;
	mat.row(19) -= fas;

	//std::cout << "Warp Correct 2: " << std::endl;
	//std::cout << mat << std::endl;

	return mat;
}

//Eigen::MatrixXd QuadPlateElement::NodeConsistentMass2()
//{
//	const double intgp = 1.0 / sqrt(3.0);
//	Eigen::Vector4d xi_list(-intgp, intgp, intgp, -intgp);
//	Eigen::Vector4d eta_list(-intgp, -intgp, intgp, intgp);
//
//	double t3 = thickness.plane_thick * thickness.plane_thick * thickness.plane_thick;
//
//	Eigen::MatrixXd M(total_dof, total_dof);
//	M.setZero();
//
//	//Eigen::Matrix3d DMat = DMatrix();
//	for (size_t i = 0; i < node_num; i++)
//	{
//		std::cout << "xi: " << xi_list(i) << ", eta: " << eta_list(i) << std::endl;
//		double detJ = JMatrix(xi_list(i), eta_list(i)).determinant();
//
//		// ui, vi, wi...
//		//Eigen::MatrixXd BMat = BMatrix(xi_list(i), eta_list(i));	
//		Eigen::Vector4d spf = ShapeFunction4(xi_list(i), eta_list(i));
//		Eigen::Vector<int, 12> shape_indices;
//		shape_indices << 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20;
//		Eigen::MatrixXd shape_mat = Eigen::MatrixXd::Zero(3, 12);
//		shape_mat << spf(0), 0, 0, spf(1), 0, 0, spf(2), 0, 0, spf(3), 0, 0,
//			0, spf(0), 0, 0, spf(1), 0, 0, spf(2), 0, 0, spf(3), 0,
//			0, 0, spf(0), 0, 0, spf(1), 0, 0, spf(2), 0, 0, spf(3);
//
//		// wi, theta_xi, theta_yi...
//		Eigen::Vector<int, 12> HVecs_indices;
//		HVecs_indices << 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22;
//		Eigen::MatrixXd H_mat = HVecs(ShapeFunctionSerendipity8(xi_list(i), eta_list(i)));
//		
//		Eigen::Matrix<double, 12, 12> MassComp;
//		MassComp = shape_mat.transpose() * shape_mat;
//		for (size_t j = 0; j < 12; j++)
//			for (size_t k = 0; k < 12; k++)
//				M(shape_indices(j), shape_indices(k)) += MassComp(j, k) * detJ * Mat.dense * thickness.plane_thick;;
//
//
//		MassComp = H_mat.transpose() * H_mat;
//		for (size_t j = 0; j < 12; j++)
//			for (size_t k = 0; k < 12; k++)
//				M(HVecs_indices(j), HVecs_indices(k)) += MassComp(j, k) * detJ * Mat.dense * t3;
//
//		//Eigen::MatrixXd BMat = BMatrix(xi_list(i), eta_list(i));
//		//M += BMat.transpose() * DMat * BMat * detJ;
//	}
//
//	Eigen::MatrixXd trMat = trans_matrix();
//	//return M;
//	return trMat.transpose() * M * trMat;
//}


Eigen::MatrixXd QuadPlateElement::NodeConsistentMass()
{
	// 2-Point Gauss Quadrature
	//const Eigen::Vector2d intg_weights(1, 1);
	//const Eigen::Vector2d intg_params(-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0));

	// 3-Point Gauss Quadrature
	const Eigen::Vector3d intg_weights(5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0);
	const Eigen::Vector3d intg_params(-sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0));

	// 4-Point Gauss Quadrature
	//const Eigen::Vector4d intg_weights(0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451);
	//const Eigen::Vector4d intg_params(-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116);
	
	// 5-Point Gauss Quadrature
	//Eigen::VectorXd intg_weights(5), intg_params(5);
	//intg_weights << 0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851;
	//intg_params << -0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459;

	double t3 = thickness.plane_thick * thickness.plane_thick * thickness.plane_thick;

	Eigen::MatrixXd M(total_dof, total_dof);
	M.setZero();
	
	//Eigen::Vector<int, 12> shape_indices;
	const int shape_indices[12] = { 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20 };
	// shape_indices << 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20;

	//Eigen::Vector<int, 12> HVecs_indices;
	//HVecs_indices << 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22;
	const int HVecs_indices[12] = { 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22 };

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

			// wi, theta_xi, theta_yi...
			Eigen::MatrixXd H_mat = HVecs(ShapeFunctionSerendipity8(intg_params(ix), intg_params(iy)));
			Eigen::Matrix<double, 12, 12> MassComp;
			MassComp = shape_mat.transpose() * shape_mat;
			for (size_t j = 0; j < 12; j++)
				for (size_t k = 0; k < 12; k++)
					M(shape_indices[j], shape_indices[k]) += MassComp(j, k) * intg_weights(ix) * intg_weights(iy) * detJ * Mat.dense * thickness.plane_thick;


			MassComp = H_mat.transpose() * H_mat;
			for (size_t j = 0; j < 12; j++)
				for (size_t k = 0; k < 12; k++)
					M(HVecs_indices[j], HVecs_indices[k]) += MassComp(j, k) * intg_weights(ix) * intg_weights(iy) * detJ * Mat.dense * t3;
		}

	}

	Eigen::MatrixXd trMat = trans_matrix();
	return trMat.transpose() * M * trMat;
}


//Eigen::MatrixXd QuadPlateElement::NodeConsistentMass2()
//{
//	const Eigen::Vector3d intg_weights(5.0/9.0, 8.0/9.0, 5.0 / 9.0);
//	const Eigen::Vector3d intg_params(-sqrt(3.0/5.0), 0, sqrt(3.0 / 5.0));
//	//Eigen::VectorXd	 xi_list(-intgp, intgp, intgp, -intgp);
//	//Eigen::Vector4d eta_list(-intgp, -intgp, intgp, intgp);
//
//	Eigen::MatrixXd M(total_dof, total_dof);
//	M.setZero();
//
//	//Eigen::Matrix3d DMat = DMatrix();
//	for (size_t ix = 0; ix < intg_params.count(); ix++)
//	{
//		for (size_t iy = 0; iy < intg_params.count(); iy++)
//		{
//			double detJ = JMatrix(intg_params(ix), intg_params(iy)).determinant();
//
//			// ui, vi, wi...
//			//Eigen::MatrixXd BMat = BMatrix(xi_list(i), eta_list(i));	
//			Eigen::Vector4d spf = ShapeFunction4(intg_params(ix), intg_params(iy));
//			Eigen::Vector<int, 12> shape_indices;
//			shape_indices << 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20;
//			Eigen::MatrixXd shape_mat = Eigen::MatrixXd::Zero(3, 12);
//			shape_mat << spf(0), 0, 0, spf(1), 0, 0, spf(2), 0, 0, spf(3), 0, 0,
//				0, spf(0), 0, 0, spf(1), 0, 0, spf(2), 0, 0, spf(3), 0,
//				0, 0, spf(0), 0, 0, spf(1), 0, 0, spf(2), 0, 0, spf(3);
//
//			// wi, theta_xi, theta_yi...
//			Eigen::Vector<int, 12> HVecs_indices;
//			HVecs_indices << 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22;
//			Eigen::MatrixXd H_mat = HVecs(ShapeFunctionSerendipity8(intg_params(ix), intg_params(iy)));
//
//			Eigen::Matrix<double, 12, 12> MassComp;
//			MassComp = shape_mat.transpose() * shape_mat;
//			for (size_t j = 0; j < 12; j++)
//				for (size_t k = 0; k < 12; k++)
//					M(shape_indices(j), shape_indices(k)) += MassComp(j, k) * intg_weights(ix) * intg_weights(iy) * detJ;
//
//
//			MassComp = H_mat.transpose() * H_mat;
//			for (size_t j = 0; j < 12; j++)
//				for (size_t k = 0; k < 12; k++)
//					M(HVecs_indices(j), HVecs_indices(k)) += MassComp(j, k) * intg_weights(ix) * intg_weights(iy) * detJ;
//
//			//Eigen::MatrixXd BMat = BMatrix(xi_list(i), eta_list(i));
//			//M += BMat.transpose() * DMat * BMat * detJ;
//		}
//
//	}
//
//	Eigen::MatrixXd trMat = trans_matrix();
//	return trMat.transpose() * M * trMat;
//}

std::vector<NodeLoadData> QuadPlateElement::InertialForceToNodeLoadData(Eigen::Vector3d accel_vec)
{
	Eigen::Matrix<double, total_dof, total_dof> m;
	m = NodeConsistentMass();
	Eigen::Vector<double, total_dof> f; f.setZero();

	f << accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0,
		accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0,
		accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0,
		accel_vec.x(), accel_vec.y(), accel_vec.z(), 0, 0, 0;

	f = m * f;
	std::vector<NodeLoadData> loads;
	loads.push_back(NodeLoadData(Nodes[0]->id, f(0), f(1), f(2), f(3), f(4), f(5)));
	loads.push_back(NodeLoadData(Nodes[1]->id, f(6), f(7), f(8), f(9), f(10), f(11)));
	loads.push_back(NodeLoadData(Nodes[2]->id, f(12), f(13), f(14), f(15), f(16), f(17)));
	loads.push_back(NodeLoadData(Nodes[3]->id, f(18), f(19), f(20), f(21), f(22), f(23)));
	return loads;
}

std::vector<NodeLoadData> QuadPlateElement::AreaForceToNodeLoadData(std::vector<Vector> load_vecs)
{
	// 2-Point Gauss Quadrature
	//const Eigen::Vector2d intg_weights(1, 1);
	//const Eigen::Vector2d intg_params(-sqrt(1.0 / 3.0), sqrt(1.0 / 3.0));

	// 3-Point Gauss Quadrature
	const Eigen::Vector3d intg_weights(5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0);
	const Eigen::Vector3d intg_params(-sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0));

	// 4-Point Gauss Quadrature
	//const Eigen::Vector4d intg_weights(0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451);
	//const Eigen::Vector4d intg_params(-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116);

	// 5-Point Gauss Quadrature
	//Eigen::VectorXd intg_weights(5), intg_params(5);
	//intg_weights << 0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851;
	//intg_params << -0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459;

	double t3 = thickness.plane_thick * thickness.plane_thick * thickness.plane_thick;

	Eigen::MatrixXd M(12, 12);
	M.setZero();

	Eigen::Vector<double, 12> p_vec;
	p_vec << load_vecs[0].x, load_vecs[0].y, load_vecs[0].z,
		load_vecs[1].x, load_vecs[1].y, load_vecs[1].z,
		load_vecs[2].x, load_vecs[2].y, load_vecs[2].z,
		load_vecs[3].x, load_vecs[3].y, load_vecs[3].z;
	//Eigen::Vector<int, 12> shape_indices;
	//const int shape_indices[12] = { 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20 };
	// shape_indices << 0, 1, 2, 6, 7, 8, 12, 13, 14, 18, 19, 20;

	//Eigen::Vector<int, 12> HVecs_indices;
	//HVecs_indices << 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22;
	//const int HVecs_indices[12] = { 2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22 };

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

double QuadPlateElement::Area()
{
	Vector v01 = Nodes[1]->Location - Nodes[0]->Location;
	Vector v03 = Nodes[3]->Location - Nodes[0]->Location;
	Vector vn13 = Vector::cross(v01, v03);

	Vector v23 = Nodes[3]->Location - Nodes[2]->Location;
	Vector v21 = Nodes[1]->Location - Nodes[2]->Location;
	Vector vn31 = Vector::cross(v01, v03);
	return (vn13.norm() + vn31.norm()) / 2;
}

Eigen::MatrixXd QuadPlateElement::StiffnessMatrix()
{
	Eigen::MatrixXd Kpln = plane_element.localStiffnessMatrix();
	Eigen::MatrixXd Kplt = localStiffnessMatrix();
	//std::cout << "K plane: \n" << Kpln << std::endl;
	//std::cout << "K plate: \n" << Kplt << std::endl;
	double max_diag = std::max(Kpln.diagonal().maxCoeff(), Kplt.diagonal().maxCoeff());

	Eigen::Matrix4d Krotz = Eigen::Matrix4d::Identity() * max_diag / 10000.0;

	//Krotz << 1, -0.5, -0.5,
	//	-0.5, 1, -0.5,
	//	-0.5, -0.5, 1;
	//Krotz *=  0.03 * Mat->Young * thickness * Area() / 1000;

	//Krotz.fill(100);
	//std::cout << "K rotz: \n" << Krotz << std::endl;
	//std::cout << "-------------------" << std::endl;
	int indices_rotz[4]{ 5,11,17,23 };
	int indices_pln[8]{ 0,1,6,7,12,13,18,19 };
	int indices_plt[12]{ 2,3,4,8,9,10,14,15,16,20,21,22 };

	Eigen::MatrixXd mat(total_dof, total_dof);
	mat.setZero();
	for (size_t i = 0; i < 8; i++)
	{
		int ir = indices_pln[i];
		for (size_t j = 0; j < 8; j++)
		{
			int ic = indices_pln[j];
			mat(ir, ic) = Kpln(i, j);
		}
	}

	for (size_t i = 0; i < 12; i++)
	{
		int ir = indices_plt[i];
		for (size_t j = 0; j < 12; j++)
		{
			int ic = indices_plt[j];
			mat(ir, ic) = Kplt(i, j);
		}
	}

	for (size_t i = 0; i < 4; i++)
	{
		int ir = indices_rotz[i];
		for (size_t j = 0; j < 4; j++)
		{
			int ic = indices_rotz[j];
			mat(ir, ic) = Krotz(i, j);
		}
	}

	//std::cout << "K: \n" << mat << std::endl;
	//std::cout << "-------------------" << std::endl;
	//std::cout << "rows: " << mat.rows() << "cols: " << mat.cols() << std::endl;
	Eigen::MatrixXd trMat = trans_matrix();

	//LocalMatrixd warp_matrix;
	//warp_matrix = WarpCorrectMatrix1b();
	//mat = warp_matrix * mat * warp_matrix.transpose();
	////matrix *= WarpCorrectMatrix().transpose(); // 不採用
	//warp_matrix = WarpCorrectMatrix2();
	//mat = warp_matrix * mat * warp_matrix.transpose();
	return trMat.transpose() * mat * trMat;
}

Eigen::MatrixXd QuadPlateElement::GeometricStiffnessMatrix(const std::vector<Displacement>& disp) {
	Eigen::MatrixXd KGpln = plane_element.geometric_local_stiffness_matrix(disp);
	Eigen::MatrixXd KGplt = geometric_local_stiffness_matrix(disp);

	int indices_pln[12]{ 0,1,2,6,7,8,12,13,14,18,19,20 };
	int indices_plt[12]{ 2,3,4,8,9,10,14,15,16,20,21,22 };

	Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(total_dof, total_dof);

	for (size_t i = 0; i < 12; i++)
	{
		int ir = indices_pln[i];
		for (size_t j = 0; j < 12; j++)
		{
			int ic = indices_pln[j];
			mat(ir, ic) += KGpln(i, j);
		}
	}

	for (size_t i = 0; i < 12; i++)
	{
		int ir = indices_plt[i];
		for (size_t j = 0; j < 12; j++)
		{
			int ic = indices_plt[j];
			mat(ir, ic) += KGplt(i, j);
		}
	}

	Eigen::MatrixXd trMat = trans_matrix();
	return trMat.transpose() * mat * trMat;
}

void QuadPlateElement::AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K)
{
	int indices[total_dof];
	for (size_t i = 0; i < node_num; i++)
		for (size_t j = 0; j < node_dof; j++)
			indices[node_dof * i + j] = Nodes[i]->id * 6 + j;

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

void QuadPlateElement::AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat)
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

void QuadPlateElement::AssembleGeometricStiffMatrix(
	Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp)
{
	AssembleMatrix(mat, GeometricStiffnessMatrix(disp));
}

void QuadPlateElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat)
{
	Eigen::VectorXd mass = NodeLumpedMass();
	for (size_t ni = 0; ni < node_num; ni++)
	{
		for (size_t i = 0; i < 3; i++) {
			int idx = Nodes[ni]->id * 6 + i;
			mat.coeffRef(idx, idx) += mass[ni];
		}
	}
}

PlateStressData QuadPlateElement::stress(Displacement d0, Displacement d1, 
	Displacement d2, Displacement d3, double xi, double eta)
{
	Eigen::Matrix3d tr = trans_matrix3(plane);
	Eigen::Matrix3d DMat = DMatrix();

	Displacement d0t = d0.translate(tr);
	Displacement d1t = d1.translate(tr);
	Displacement d2t = d2.translate(tr);
	Displacement d3t = d3.translate(tr);

	Eigen::VectorXd wvec(total_dof_local);
	wvec << d0t.Dz(), d0t.Rx(), d0t.Ry(), d1t.Dz(), d1t.Rx(), d1t.Ry(), 
		d2t.Dz(), d2t.Rx(), d2t.Ry(), d3t.Dz(), d3t.Rx(), d3t.Ry();

	// Moment
	Eigen::Vector3d strs = DMat * BMatrix(xi, eta) * wvec;

	// Shear
	Eigen::Matrix2d DJinv_dxi = dJinv_dxi(xi, eta);
	Eigen::Matrix2d DJinv_deta = dJinv_deta(xi, eta);
	Eigen::Matrix2d Jinv = JMatrix(xi, eta).inverse();

	Eigen::MatrixXd dH_dxi = HVecs(dndxi(xi, eta));
	Eigen::MatrixXd dH_deta = HVecs(dndeta(xi, eta));

	Eigen::MatrixXd dH_d2xi = HVecs(dn2_dxi2(xi, eta));
	Eigen::MatrixXd dH_d2eta = HVecs(dn2_deta2(xi, eta));
	Eigen::MatrixXd dH_deta_dxi = HVecs(dn2_dxieta(xi, eta));

	Eigen::MatrixXd dB_dxi(3, 12), dB_deta(3, 12);
	dB_dxi.row(0) = DJinv_dxi(0, 0) * dH_dxi.row(0) + Jinv(0, 0) * dH_d2xi.row(0)
		+ DJinv_dxi(0, 1) * dH_deta.row(0) + Jinv(0, 1) * dH_deta_dxi.row(0);
	dB_dxi.row(1) = DJinv_dxi(1, 0) * dH_dxi.row(1) + Jinv(1, 0) * dH_d2xi.row(1)
		+ DJinv_dxi(1, 1) * dH_deta.row(1) + Jinv(1, 1) * dH_deta_dxi.row(1);
	dB_dxi.row(2) = DJinv_dxi(0, 0) * dH_dxi.row(1) + Jinv(0, 0) * dH_d2xi.row(1)
		+ DJinv_dxi(0, 1) * dH_deta.row(1) + Jinv(0, 1) * dH_deta_dxi.row(1)
		+ DJinv_dxi(1, 0) * dH_dxi.row(0) + Jinv(1, 0) * dH_d2xi.row(0)
		+ DJinv_dxi(1, 1) * dH_deta.row(0) + Jinv(1, 1) * dH_deta_dxi.row(0);

	dB_deta.row(0) = DJinv_deta(0, 0) * dH_dxi.row(0) + Jinv(0, 0) * dH_deta_dxi.row(0)
		+ DJinv_deta(0, 1) * dH_deta.row(0) + Jinv(0, 1) * dH_d2eta.row(0);
	dB_deta.row(1) = DJinv_deta(1, 0) * dH_dxi.row(1) + Jinv(1, 0) * dH_deta_dxi.row(1)
		+ DJinv_deta(1, 1) * dH_deta.row(1) + Jinv(1, 1) * dH_d2eta.row(1);
	dB_deta.row(2) = DJinv_deta(0, 0) * dH_dxi.row(1) + Jinv(0, 0) * dH_deta_dxi.row(1)
		+ DJinv_deta(0, 1) * dH_deta.row(1) + Jinv(0, 1) * dH_d2eta.row(1)
		+ DJinv_deta(1, 0) * dH_dxi.row(0) + Jinv(1, 0) * dH_deta_dxi.row(0)
		+ DJinv_deta(1, 1) * dH_deta.row(0) + Jinv(1, 1) * dH_d2eta.row(0);

	Eigen::MatrixXd dB_dx(3, 12), dB_dy(3, 12);
	dB_dx = Jinv(0, 0) * dB_dxi + Jinv(0, 1) * dB_deta;
	dB_dy = Jinv(1, 0) * dB_dxi + Jinv(1, 1) * dB_deta;

	Eigen::Vector3d dM_dx = DMat * dB_dx * wvec;
	Eigen::Vector3d dM_dy = DMat * dB_dy * wvec;

	double qx = dM_dx(0) + dM_dy(2);
	double qy = dM_dy(1) + dM_dx(2);

	// Plane Composition
	MembraneStressData mstr = plane_element.stress(d0, d1, d2, d3, xi, eta);

	return PlateStressData(strs(0), strs(1), strs(2), qx, qy,
		mstr.sigx * thickness.plane_thick, mstr.sigy * thickness.plane_thick, mstr.sigxy * thickness.plane_thick);
}

//void QuadPlateElement::shearstress(Displacement d0, Displacement d1, 
//	Displacement d2, Displacement d3, double xi, double eta)
//{
//	Eigen::Matrix2d DJinv_dxi = dJinv_dxi(xi, eta);
//	Eigen::Matrix2d DJinv_deta = dJinv_deta(xi, eta);
//	Eigen::Matrix2d Jinv = JMatrix(xi, eta).inverse();
//
//	Eigen::MatrixXd dH_dxi = HVecs(dndxi(xi, eta));
//	Eigen::MatrixXd dH_deta = HVecs(dndeta(xi, eta));
//
//	Eigen::MatrixXd dH_d2xi = HVecs(dn2_dxi2(xi, eta));
//	Eigen::MatrixXd dH_d2eta = HVecs(dn2_deta2(xi, eta));
//	Eigen::MatrixXd dH_deta_dxi = HVecs(dn2_dxieta(xi, eta));
//
//	Eigen::MatrixXd dB_dxi(3, 12), dB_deta(3, 12);
//	dB_dxi.row(0) = DJinv_dxi(0, 0) * dH_dxi.row(0) + Jinv(0, 0) * dH_d2xi.row(0)
//		+ DJinv_dxi(0, 1) * dH_deta.row(0) + Jinv(0, 1) * dH_deta_dxi.row(0);
//	dB_dxi.row(1) = DJinv_dxi(1, 0) * dH_dxi.row(1) + Jinv(1, 0) * dH_d2xi.row(1) 
//		+ DJinv_dxi(1, 1) * dH_deta.row(1) + Jinv(1, 1) * dH_deta_dxi.row(1);
//	dB_dxi.row(2) = DJinv_dxi(0, 0) * dH_dxi.row(1) + Jinv(0, 0) * dH_d2xi.row(1)
//		+ DJinv_dxi(0, 1) * dH_deta.row(1) + Jinv(0, 1) * dH_deta_dxi.row(1)
//		+ DJinv_dxi(1, 0) * dH_dxi.row(0) + Jinv(1, 0) * dH_d2xi.row(0)
//		+ DJinv_dxi(1, 1) * dH_deta.row(0) + Jinv(1, 1) * dH_deta_dxi.row(0);
//
//	dB_deta.row(0) = DJinv_deta(0, 0) * dH_dxi.row(0) + Jinv(0, 0) * dH_deta_dxi.row(0)
//		+ DJinv_deta(0, 1) * dH_deta.row(0) + Jinv(0, 1) * dH_d2eta.row(0);
//	dB_deta.row(1) = DJinv_deta(1, 0) * dH_dxi.row(1) + Jinv(1, 0) * dH_deta_dxi.row(1)
//		+ DJinv_deta(1, 1) * dH_deta.row(1) + Jinv(1, 1) * dH_d2eta.row(1);
//	dB_deta.row(2) = DJinv_deta(0, 0) * dH_dxi.row(1) + Jinv(0, 0) * dH_deta_dxi.row(1)
//		+ DJinv_deta(0, 1) * dH_deta.row(1) + Jinv(0, 1) * dH_d2eta.row(1)
//		+ DJinv_deta(1, 0) * dH_dxi.row(0) + Jinv(1, 0) * dH_deta_dxi.row(0)
//		+ DJinv_deta(1, 1) * dH_deta.row(0) + Jinv(1, 1) * dH_d2eta.row(0);
//
//	Eigen::MatrixXd dB_dx(3, 12), dB_dy(3, 12);
//	dB_dx = Jinv(0, 0) * dB_dxi + Jinv(0, 1) * dB_deta;
//	dB_dy = Jinv(1, 0) * dB_dxi + Jinv(1, 1) * dB_deta;
//
//	Eigen::Matrix3d DMat = DMatrix();
//	Eigen::Matrix3d tr = trans_matrix3(plane);
//
//	Displacement d0t = d0.translate(tr);
//	Displacement d1t = d1.translate(tr);
//	Displacement d2t = d2.translate(tr);
//	Displacement d3t = d3.translate(tr);
//
//	Eigen::VectorXd wvec(total_dof_local);
//	wvec << d0t.Dz(), d0t.Rx(), d0t.Ry(), d1t.Dz(), d1t.Rx(), d1t.Ry(),
//		d2t.Dz(), d2t.Rx(), d2t.Ry(), d3t.Dz(), d3t.Rx(), d3t.Ry();
//
//	Eigen::Vector3d dM_dx = DMatrix() * dB_dx * wvec;
//	Eigen::Vector3d dM_dy = DMatrix() * dB_dy * wvec;
//
//	double qx = dM_dx(0) + dM_dy(2);
//	double qy = dM_dy(1) + dM_dx(2);
//
//	// std::cout << "M_x: " << dM_dx << ", M_y: " << dM_dy << std::endl;
//
//	std::cout << "qx: " << qx << ", qy: " << qy << std::endl;
//
//	// return PlateStressData(strs(0), strs(1), strs(2));
//}


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

Eigen::MatrixXd ComplexBeamElement::stiffness_matrix_local()
{
	double l = length();
	double lz_ = l - lzi - lzj;
	double ly_ = l - lyi - lyj;

	Eigen::Matrix4d Kbz = stiffness_matrix_beam_rot_z(Mat.Young, Sec->Iz, lz_);
	Eigen::Matrix4d Tbz = Eigen::Matrix4d::Identity();
	Tbz(0, 1) = lzi; Tbz(2, 3) = -lzj;

	Eigen::Matrix4d Kby = stiffness_matrix_beam_rot_y(Mat.Young, Sec->Iy, ly_);
	Eigen::Matrix4d Tby = Eigen::Matrix4d::Identity();
	Tby(0, 1) = -lzi; Tby(2, 3) = lzj;

	Eigen::Matrix2d Kx = stiffness_matrix_truss(Mat.Young, Sec->A, l);

	Eigen::Matrix2d Kt = stiffness_matrix_beam_rot_x(Mat.G(), Sec->K, l);

	// 端部バネの考慮
	Kbz = compute_Kprime_partial(Kbz, Lambda_sy, Lambda_bz, Lambda_sy_, Lambda_bz_);
	Kby = compute_Kprime_partial(Kby, Lambda_sz, Lambda_by, Lambda_sz_, Lambda_by_);

	// 剛域の考慮
	Kbz = Tbz.transpose() * Kbz * Tbz;
	Kby = Tby.transpose() * Kby * Tby;

	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(12, 12);
	K(0, 0) = Kx(0, 0); K(0, 6) = Kx(0, 1);
	K(6, 0) = Kx(1, 0); K(6, 6) = Kx(1, 1);
	
	K(1, 1) = Kbz(0, 0); K(1, 5) = Kbz(0, 1); K(1, 7) = Kbz(0, 2); K(1, 11) = Kbz(0, 3);
	K(5, 1) = Kbz(1, 0); K(5, 5) = Kbz(1, 1); K(5, 7) = Kbz(1, 2); K(5, 11) = Kbz(1, 3);
	K(7, 1) = Kbz(2, 0); K(7, 5) = Kbz(2, 1); K(7, 7) = Kbz(2, 2); K(7, 11) = Kbz(2, 3);
	K(11, 1) = Kbz(3, 0); K(11, 5) = Kbz(3, 1); K(11, 7) = Kbz(3, 2); K(11, 11) = Kbz(3, 3);
	
	K(2, 2) = Kby(0, 0); K(2, 4) = Kby(0, 1); K(2, 8) = Kby(0, 2); K(1, 10) = Kby(0, 3);
	K(4, 2) = Kby(1, 0); K(4, 4) = Kby(1, 1); K(4, 8) = Kby(1, 2); K(4, 10) = Kby(1, 3);
	K(8, 2) = Kby(2, 0); K(8, 4) = Kby(2, 1); K(8, 8) = Kby(2, 2); K(8, 10) = Kby(2, 3);
	K(10, 2) = Kby(3, 0); K(10, 4) = Kby(3, 1); K(10, 8) = Kby(3, 2); K(10, 10) = Kby(3, 3);

	K(3, 3) = Kt(0, 0); K(3, 9) = Kt(0, 1);
	K(9, 3) = Kt(1, 0); K(9, 9) = Kt(1, 1);
	
	return K;
}

ComplexBeamElement::ComplexBeamElement()
{
	Lambda_bz = 1; Lambda_bz_ = 1;
	Lambda_sy = 1; Lambda_sy_ = 1;
	Lambda_by = 1; Lambda_by_ = 1;
	Lambda_sz = 1; Lambda_sz_ = 1;
	
	lzi = 0; lzj = 0;
	lyi = 0; lyj = 0;
}

ComplexBeamElement::ComplexBeamElement(Node* n0, Node* n1, Section* sec, Material mat, double beta)
	:BeamElement(n0, n1, sec, mat, beta)
{
	Lambda_bz = 1; Lambda_bz_ = 1;
	Lambda_sy = 1; Lambda_sy_ = 1;
	Lambda_by = 1; Lambda_by_ = 1;
	Lambda_sz = 1; Lambda_sz_ = 1;

	lzi = 0; lzj = 0;
	lyi = 0; lyj = 0;
}

ComplexBeamElement::ComplexBeamElement(int _id, Node* n0, Node* n1, Section* sec, Material mat, double beta)
	:ComplexBeamElement(n0, n1, sec, mat, beta)
{
	id = _id;

	Lambda_bz = 1; Lambda_bz_ = 1;
	Lambda_sy = 1; Lambda_sy_ = 1;
	Lambda_by = 1; Lambda_by_ = 1;
	Lambda_sz = 1; Lambda_sz_ = 1;

	lzi = 0; lzj = 0;
	lyi = 0; lyj = 0;
}

Eigen::MatrixXd ComplexBeamElement::StiffnessMatrix() {
	Eigen::MatrixXd tr = trans_matrix();
	Eigen::MatrixXd k = stiffness_matrix_local();
	return tr.transpose() * k * tr;
}

//Eigen::MatrixXd ComplexBeamElement::NodeConsistentMass()
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
//}

double BarElementBase::length()
{
	return Nodes[0]->Location.distance_to(Nodes[1]->Location);
}

void BarElementBase::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat)
{
	int indices[6];
	Eigen::VectorXd mass = NodeLumpedMass();
	for (size_t i = 0; i < 3; i++) {
		int idx = Nodes[0]->id * 6 + i;
		mat.coeffRef(idx, idx) += mass[0];
	}
	for (size_t i = 0; i < 3; i++) {
		int idx = Nodes[1]->id * 6 + i;
		mat.coeffRef(idx, idx) += mass[1];
	}
}

PlatePrincipalStressData::PlatePrincipalStressData(PlateStressData psd)
{
	double dm = sqrt(pow(((psd.Mx - psd.My) * 0.5), 2) + psd.Mxy * psd.Mxy);
	double m1 = (psd.Mx + psd.My) * 0.5 + dm;
	double m2 = (psd.Mx + psd.My) * 0.5 - dm;
	
	double theta_m = 0;
	if (abs(psd.Mx - psd.My) > 0.00001) {
		double t2m = 2 * psd.Mxy / (psd.Mx - psd.My);
		theta_m = atan2(2 * psd.Mxy, (psd.Mx - psd.My)) / 2;
		//theta_m = atan(t2m) / 2;
	}
	
	double dn = sqrt(pow(((psd.Nx - psd.Ny) * 0.5), 2) + psd.Qxy * psd.Qxy);
	double n1 = (psd.Nx + psd.Ny) * 0.5 + dn;
	double n2 = (psd.Nx + psd.Ny) * 0.5 - dn;

	double theta_n = 0;
	if (abs(psd.Nx - psd.Ny) > 0.00001) {
		double t2n = 2 * psd.Qxy / (psd.Nx - psd.Ny);
		theta_n = atan2(2 * psd.Qxy, (psd.Nx - psd.Ny)) / 2;
		//theta_n = atan(t2n) / 2;
	}

	M1 = m1;
	M2 = m2;
	theta_m1 = theta_m;

	N1 = n1;
	N2 = n2;
	theta_n1 = theta_n;
}

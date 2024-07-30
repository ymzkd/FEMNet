#include <Eigen/Sparse>

#include "Components.h"
#include "Element.h"


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
	double EAl = Mat->Young * Sec->A / l;
	//double alpha = vij * Vector::XAxis() / l;
	//double beta = vij * Vector::YAxis() / l;
	//double gamma = vij * Vector::ZAxis() / l;

	Eigen::Matrix2d local_stiffMat;
	local_stiffMat << EAl, -EAl,
		-EAl, EAl;
	return local_stiffMat;
}

TrussElement::TrussElement(Node* n0, Node* n1, Section* sec, Material* mat)
{
	Nodes[0] = n0;
	Nodes[1] = n1;

	Sec = sec;
	Mat = mat;
}

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

void TrussElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat)
{
	int indices[6];
	for (size_t i = 0; i < 3; i++)
		indices[i] = Nodes[0]->id * 6 + i;
	for (size_t i = 0; i < 3; i++)
		indices[i + 3] = Nodes[1]->id * 6 + i;

	Eigen::MatrixXd smat = StiffnessMatrix();
	for (size_t i = 0; i < 6; i++)
	{

		for (size_t j = 0; j < i + 1; j++)
		{
			int ci, rj;
			if (indices[j] <= indices[i]) {
				ci = indices[i];
				rj = indices[j];
			}
			else {
				ci = indices[j];
				rj = indices[i];
			}
			mat.coeffRef(rj, ci) += smat(i, j);
		}
	}
}

BeamStress TrussElement::stress(Displacement d0, Displacement d1)
{
	Eigen::VectorXd disp(6);
	disp << d0.Dx(), d0.Dy(), d0.Dz(), d1.Dx(), d1.Dy(), d1.Dz();

	Eigen::Vector2d force = stiffness_matrix_local() * trans_matrix()* disp;
	
	BeamStressData str0(force(0), 0, 0, 0, 0, 0);
	BeamStressData str1(force(1), 0, 0, 0, 0, 0);
	return BeamStress(str0, str1);
}


BeamElement::BeamElement(Node* n0, Node* n1, Section* sec, Material* mat, double beta):
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

void BeamElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat)
{
	int indices[12];
	for (size_t i = 0; i < 6; i++)
		indices[i] = Nodes[0]->id * 6 + i;
	for (size_t i = 0; i < 6; i++)
		indices[i + 6] = Nodes[1]->id * 6 + i;
	
	Eigen::MatrixXd smat = StiffnessMatrix();
	for (size_t i = 0; i < 12; i++)
	{
		
		for (size_t j = 0; j < i + 1; j++)
		{
			int ci, rj;
			if (indices[j] <= indices[i]) {
				ci = indices[i];
				rj = indices[j];
			}
			else {
				ci = indices[j];
				rj = indices[i];
			}
			mat.coeffRef(rj,ci) += smat(i, j);
		}
	}
}

BeamStress BeamElement::stress(Displacement d0, Displacement d1)
{
	Eigen::VectorXd wvec(12);
	wvec.segment(0, 6) = Eigen::Map<Eigen::VectorXd>(d0.displace, 6);
	wvec.segment(6, 6) = Eigen::Map<Eigen::VectorXd>(d1.displace, 6);
	
	wvec = trans_matrix() * wvec;
	wvec = stiffness_matrix_local() * wvec;
	//wvec = StiffnessMatrix() * wvec;

	BeamStressData s0(wvec(0), wvec(1), wvec(2), wvec(3), wvec(4), wvec(5));
	BeamStressData s1(wvec(6), wvec(7), wvec(8), wvec(9), wvec(10), wvec(11));
	return BeamStress(s0, s1);
}

double BeamElement::element_length()
{
	Point p0 = Nodes[0]->Location;
	Point p1 = Nodes[1]->Location;
	double l2 = (p1.x - p0.x)* (p1.x - p0.x) + (p1.y - p0.y) * (p1.y - p0.y) + (p1.z - p0.z) * (p1.z - p0.z);
	return sqrt(l2);
}

Eigen::MatrixXd BeamElement::stiffness_matrix_local() {
	double E = Mat->Young;
	double G = Mat->G();
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
	double t, Material* mat)
{
	Nodes[0] = n0;
	Nodes[1] = n1;
	Nodes[2] = n2;
	
	thickness = t;
	Mat = mat;
	plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
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

void TriPlaneElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat)
{
	int indices[total_dof];
	for (size_t i = 0; i < node_dof; i++)
		indices[i] = Nodes[0]->id * 6 + i;
	for (size_t i = 0; i < node_dof; i++)
		indices[i + node_dof] = Nodes[1]->id * 6 + i;
	for (size_t i = 0; i < node_dof; i++)
		indices[i + node_dof * 2] = Nodes[2]->id * 6 + i;

	Eigen::MatrixXd smat = StiffnessMatrix();
	for (size_t i = 0; i < total_dof; i++)
	{

		for (size_t j = 0; j < i + 1; j++)
		{
			int ci, rj;
			if (indices[j] <= indices[i]) {
				ci = indices[i];
				rj = indices[j];
			}
			else {
				ci = indices[j];
				rj = indices[i];
			}
			mat.coeffRef(rj, ci) += smat(i, j);
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
	mat << 1, Mat->Poisson, 0,
		Mat->Poisson, 1, 0,
		0, 0, (1 - Mat->Poisson) / 2;
	return Mat->Young / (1 - Mat->Poisson * Mat->Poisson) * mat;
}

Eigen::MatrixXd TriPlaneElement::localStiffnessMatrix()
{
	Eigen::MatrixXd BMat = BMatrix();
	Eigen::MatrixXd DMat = DMatrix();

	Eigen::MatrixXd K = thickness * Area() * BMat.transpose() * DMat * BMat;
	return K;
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

Eigen::MatrixXd TriPlateElement::BMatrix(double xi, double eta)
{
	Point p1 = plane.PointToCoord(Nodes[0]->Location);
	Point p2 = plane.PointToCoord(Nodes[1]->Location);
	Point p3 = plane.PointToCoord(Nodes[2]->Location);

	Vector v23 = p2 - p3;
	Vector v31 = p3 - p1;
	Vector v12 = p1 - p2;

	double p4 = -6 * v23.x / v23.squared_norm();
	double p5 = -6 * v31.x / v31.squared_norm();
	double p6 = -6 * v12.x / v12.squared_norm();

	double q4 = 3 * v23.x * v23.y / v23.squared_norm();
	double q5 = 3 * v31.x * v31.y / v31.squared_norm();
	double q6 = 3 * v12.x * v12.y / v12.squared_norm();
	
	double r4 = 3 * v23.y * v23.y / v23.squared_norm();
	double r5 = 3 * v31.y * v31.y / v31.squared_norm();
	double r6 = 3 * v12.y * v12.y / v12.squared_norm();

	double t4 = -6 * v23.y / v23.squared_norm();
	double t5 = -6 * v31.y / v31.squared_norm();
	double t6 = -6 * v12.y / v12.squared_norm();

	Eigen::Vector<double, 9> Hx_xi, Hy_xi, Hx_eta, Hy_eta;
	Hx_xi << 
		p6 * (1 - 2 * xi) + (p5 - p6) * eta,
		q6 * (1 - 2 * xi) - (q5 + q6) * eta,
		-4 + 6 * (xi + eta) + r6 * (1 - 2 * xi) - eta * (r5 + r6),
		-p6 * (1 - 2 * xi) + eta * (p4 + p6),
		q6 * (1 - 2 * xi) - eta * (q6 - q4),
		-2 + 6 * xi + r6 * (1 - 2 * xi) + eta * (r4 - r6),
		-eta * (p5 + p4),
		eta * (q4 - q5),
		-eta * (r5 - r4);
	
	// 元論文
	Hy_xi <<
		t6 * (1 - 2 * xi) + (t5 - t6) * eta,
		1 + r6 * (1 - 2 * xi) - (r5 + r6) * eta,
		-q6 * (1 - 2 * xi) + eta * (q5 + q6),
		-t6 * (1 - 2 * xi) + eta * (t4 + t6),
		-1 + r6 * (1 - 2 * xi) + eta * (r4 - r6),
		-q6 * (1 - 2 * xi) - eta * (q4 - q6),
		-eta * (t4 + t5),
		eta * (r4 - r5),
		-eta * (q4 - q5);

	Hx_eta <<
		-p5 * (1 - 2 * eta) - xi * (p6 - p5),
		q5 * (1 - 2 * eta) - xi * (q5 + q6),
		-4 + 6 * (xi + eta) + r5 * (1 - 2 * eta) - xi * (r6 + r5),
		xi * (p4 + p6),
		xi * (q4 - q6),
		-xi * (r6 - r4),
		p5 * (1 - 2 * eta) - xi * (p4 + p5),
		q5 * (1 - 2 * eta) + xi * (q4 - q5),
		-2 + 6 * eta + r5 * (1 - 2 * eta) + xi * (r4 - r5);

	// 元論文
	Hy_eta <<
		-t5 * (1 - 2 * eta) - xi * (t6 - t5),
		1 + r5 * (1 - 2 * eta) - xi * (r5 + r6),
		-q5 * (1 - 2 * eta) + xi * (q5 + q6),
		xi * (t4 + t6),
		xi * (r4 - r6),
		-xi * (q4 - q6),
		t5 * (1 - 2 * eta) - xi * (t4 + t5),
		-1 + r5 * (1 - 2 * eta) + xi * (r4 - r5),
		-q5 * (1 - 2 * eta) - xi * (q4 - q5);

	Eigen::MatrixXd mat(3, 9);
	mat.row(0) = v31.y * Hx_xi + v12.y * Hx_eta;
	mat.row(1) = -v31.x * Hy_xi - v12.x * Hy_eta;
	mat.row(2) = -v31.x * Hx_xi - v12.x * Hx_eta + v31.y * Hy_xi + v12.y * Hy_eta;
	mat /= (2.0 * Area());
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
	mat << 1, Mat->Poisson, 0,
		Mat->Poisson, 1, 0,
		0, 0, (1 - Mat->Poisson) / 2;
	return Mat->Young * pow(thickness, 3) / 12 / (1 - Mat->Poisson * Mat->Poisson) * mat;
}

Eigen::MatrixXd TriPlateElement::localStiffnessMatrix()
{
	double weight = 1.0 / 3.0;
	double xi_list[3]{ 0.5 , 0.5 , 0 };
	double eta_list[3]{ 0, 0.5, 0.5 };
	double a2 = Area();
	//double a2 = 2.0 * Area();

	Eigen::Matrix3d DMat = DMatrix();
	Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(9, 9);
	for (size_t i = 0; i < 3; i++)
	{
		Eigen::MatrixXd BMat = BMatrix(xi_list[i], eta_list[i]);
		mat += BMat.transpose() * DMat * BMat;
	}
	return a2 * weight * mat;
}

TriPlateElement::TriPlateElement(Node* n0, Node* n1, Node* n2, double t, Material* mat) 
	: plane_element(n0, n1, n2, t, mat)
{
	Nodes[0] = n0;
	Nodes[1] = n1;
	Nodes[2] = n2;

	thickness = t;
	Mat = mat;
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

void TriPlateElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat)
{
	int indices[total_dof];
	for (size_t i = 0; i < node_dof; i++)
		indices[i] = Nodes[0]->id * 6 + i;
	for (size_t i = 0; i < node_dof; i++)
		indices[i + node_dof] = Nodes[1]->id * 6 + i;
	for (size_t i = 0; i < node_dof; i++)
		indices[i + node_dof * 2] = Nodes[2]->id * 6 + i;

	Eigen::MatrixXd smat = StiffnessMatrix();
	for (size_t i = 0; i < total_dof; i++)
	{

		for (size_t j = 0; j < i + 1; j++)
		{
			int ci, rj;
			if (indices[j] <= indices[i]) {
				ci = indices[i];
				rj = indices[j];
			}
			else {
				ci = indices[j];
				rj = indices[i];
			}
			mat.coeffRef(rj, ci) += smat(i, j);
		}
	}
}

PlateStressData TriPlateElement::stress(
	Displacement d0, Displacement d1, Displacement d2, double xi, double eta)
{
	Eigen::Matrix3d tr = trans_matrix3(plane);

	Displacement d0t = d0.translate(tr);
	Displacement d1t = d1.translate(tr);
	Displacement d2t = d2.translate(tr);
	
	// Moments
	Eigen::Matrix3d Dmat = DMatrix();
	Eigen::VectorXd wvec(9);
	wvec << d0t.Dz(), d0t.Rx(), d0t.Ry(), d1t.Dz(), d1t.Rx(), d1t.Ry(), d2t.Dz(), d2t.Rx(), d2t.Ry();
	Eigen::Vector3d strs = Dmat * BMatrix(xi, eta) * wvec;
	
	// Shear
	Point p1 = plane.PointToCoord(Nodes[0]->Location);
	Point p2 = plane.PointToCoord(Nodes[1]->Location);
	Point p3 = plane.PointToCoord(Nodes[2]->Location);

	Vector v23 = p2 - p3;
	Vector v31 = p3 - p1;
	Vector v12 = p1 - p2;

	double p4 = -6 * v23.x / v23.squared_norm();
	double p5 = -6 * v31.x / v31.squared_norm();
	double p6 = -6 * v12.x / v12.squared_norm();

	double q4 = 3 * v23.x * v23.y / v23.squared_norm();
	double q5 = 3 * v31.x * v31.y / v31.squared_norm();
	double q6 = 3 * v12.x * v12.y / v12.squared_norm();

	double r4 = 3 * v23.y * v23.y / v23.squared_norm();
	double r5 = 3 * v31.y * v31.y / v31.squared_norm();
	double r6 = 3 * v12.y * v12.y / v12.squared_norm();

	double t4 = -6 * v23.y / v23.squared_norm();
	double t5 = -6 * v31.y / v31.squared_norm();
	double t6 = -6 * v12.y / v12.squared_norm();

	double Jl = v12.y * v31.x - v12.x * v31.y;
	Eigen::Matrix2d Jinv;
	Jinv << v31.y / Jl, v12.y / Jl, -v31.x / Jl, -v12.x / Jl;

	Eigen::Vector<double, 9> Hy_xi_xi, Hy_xi_eta, Hy_eta_xi, Hy_eta_eta;
	Hy_xi_xi << -2.0 * t6, -2.0 * r6, 2.0 * q6, 2.0 * t6, -2.0 * r6, 2.0 * q6, 0, 0, 0;
	Hy_xi_eta << t5 - t6, -(r5 + r6), q5 + q6, t4 + t6, r4 - r6, -(q4 - q6), -(t4 + t5), r4 - r5, q5 - q4;
	Hy_eta_xi << t5 - t6, -(r5 + r6), q5 + q6, t4 + t6, r4 - r6, -(q4 - q6), -(t4 + t5), r4 - r5, q5 - q4;
	Hy_eta_eta << 2.0 * t5, -2.0 * r5, 2.0 * q5, 0, 0, 0, -2.0 * t5, -2.0 * r5, 2.0 * q5;

	Eigen::Vector<double, 9> Hx_xi_xi, Hx_xi_eta, Hx_eta_xi, Hx_eta_eta;
	Hx_xi_xi << -2.0 * p6, -2.0 * q6, 6.0 - 2.0 * r6, 2.0 * p6, -2.0 * q6, 6.0 - 2.0 * r6, 0, 0, 0;
	Hx_xi_eta << p5 - p6, -(q5 + q6), 6.0 - (r5 + r6), p4 + p6, q4 - q6, r4 - r6, -(p5 + p4), (q4 - q5), -(r5 - r4);
	Hx_eta_xi << p5 - p6, -(q5 + q6), 6.0 - (r5 + r6), p4 + p6, q4 - q6, r4 - r6, -(p4 + p5), q4 - q5, r4 - r5;
	Hx_eta_eta << 2.0 * p5, -2.0 * q5, 6.0 - 2.0 * r5, 0, 0, 0, -2.0 * p5, -2.0 * q5, 6.0 - 2.0 * r5;

	Eigen::Vector<double, 9> Hy_xi_x, Hy_xi_y, Hy_eta_x, Hy_eta_y;
	Hy_xi_x = Jinv(0, 0) * Hy_xi_xi + Jinv(0, 1) * Hy_xi_eta;
	Hy_xi_y = Jinv(1, 0) * Hy_xi_xi + Jinv(1, 1) * Hy_xi_eta;
	Hy_eta_x = Jinv(0, 0) * Hy_eta_xi + Jinv(0, 1) * Hy_eta_eta;
	Hy_eta_y = Jinv(1, 0) * Hy_eta_xi + Jinv(1, 1) * Hy_eta_eta;

	Eigen::Vector<double, 9> Hx_xi_x, Hx_xi_y, Hx_eta_x, Hx_eta_y;
	Hx_xi_x = Jinv(0, 0) * Hx_xi_xi + Jinv(0, 1) * Hx_xi_eta;
	Hx_xi_y = Jinv(1, 0) * Hx_xi_xi + Jinv(1, 1) * Hx_xi_eta;
	Hx_eta_x = Jinv(0, 0) * Hx_eta_xi + Jinv(0, 1) * Hx_eta_eta;
	Hx_eta_y = Jinv(1, 0) * Hx_eta_xi + Jinv(1, 1) * Hx_eta_eta;

	Eigen::MatrixXd Bmat_y(3, 9);
	Bmat_y.row(0) = v31.y * Hx_xi_y + v12.y * Hx_eta_y;
	Bmat_y.row(1) = -v31.x * Hy_xi_y - v12.x * Hy_eta_y;
	Bmat_y.row(2) = -v31.x * Hx_xi_y - v12.x * Hx_eta_y + v31.y * Hy_xi_y + v12.y * Hy_eta_y;
	Bmat_y /= (2.0 * Area());

	Eigen::MatrixXd Bmat_x(3, 9);
	Bmat_x.row(0) = v31.y * Hx_xi_x + v12.y * Hx_eta_x;
	Bmat_x.row(1) = -v31.x * Hy_xi_x - v12.x * Hy_eta_x;
	Bmat_x.row(2) = -v31.x * Hx_xi_x - v12.x * Hx_eta_x + v31.y * Hy_xi_x + v12.y * Hy_eta_x;
	Bmat_x /= (2.0 * Area());

	Eigen::Vector3d dMdx = Dmat * (Bmat_x * wvec);
	Eigen::Vector3d dMdy = Dmat * (Bmat_y * wvec);

	double qx = dMdx(0) + dMdy(2);
	double qy = dMdy(1) + dMdx(2);

	return PlateStressData(strs[0], strs[1], strs[2], qx, qy);
}

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
	mat << 1, Mat->Poisson, 0,
		Mat->Poisson, 1, 0,
		0, 0, (1 - Mat->Poisson) / 2.0;
	return Mat->Young / (1.0 - Mat->Poisson * Mat->Poisson) * mat;
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
		K += thickness * BMat.transpose() * DMat * BMat * detJ;
	}
	return K;
}

QuadPlaneElement::QuadPlaneElement(Node* n0, Node* n1, Node* n2, Node* n3, double t, Material* mat)
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
	// plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
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

Eigen::MatrixXd QuadPlaneElement::StiffnessMatrix()
{
	Eigen::MatrixXd K = localStiffnessMatrix();
	Eigen::MatrixXd TrMat = trans_matrix();
	return TrMat.transpose() * K * TrMat;
}

void QuadPlaneElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat)
{
	int indices[total_dof];
	for (size_t i = 0; i < node_num; i++)
		for (size_t j = 0; j < node_dof; j++)
			indices[node_dof * i + j] = Nodes[i]->id * 6 + j;
	
	Eigen::MatrixXd smat = StiffnessMatrix();
	for (size_t i = 0; i < total_dof; i++)
	{
		for (size_t j = 0; j < i + 1; j++)
		{
			int ci, rj;
			if (indices[j] <= indices[i]) {
				ci = indices[i];
				rj = indices[j];
			}
			else {
				ci = indices[j];
				rj = indices[i];
			}
			mat.coeffRef(rj, ci) += smat(i, j);
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

	Eigen::Matrix2d dJinv_dxi;
	dJinv_dxi = 0.25 / detJ * Matdfdxi - ddetJ_dxi / detJ * JMat.inverse();
	
	return dJinv_dxi;
}

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

	Eigen::Matrix2d dJinv_deta;
	dJinv_deta = 0.25 / detJ * Matdfdeta - ddetJ_deta / detJ * JMat.inverse();

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

	// Eigen::Vector<double, 8> shape_funcs;

	Eigen::Matrix<double, 2, total_dof_local> HVecs;
	// Eigen::Vector<double, total_dof_local> Hx_x, Hx_eta, Hy_xi, Hy_eta;
	
	HVecs.row(0) << 1.5 * (a5 * shape_funcs(4) - a8 * shape_funcs(7)),
		b5* shape_funcs(4) + b8 * shape_funcs(7),
		shape_funcs(0) - c5 * shape_funcs(4) - c8 * shape_funcs(7),
		1.5 * (a6 * shape_funcs(5) - a5 * shape_funcs(4)),
		b6* shape_funcs(5) + b5 * shape_funcs(4),
		shape_funcs(1) - c6 * shape_funcs(5) - c5 * shape_funcs(4),
		1.5 * (a7 * shape_funcs(6) - a6 * shape_funcs(5)),
		b7* shape_funcs(6) + b6 * shape_funcs(5),
		shape_funcs(2) - c7 * shape_funcs(6) - c6 * shape_funcs(5),
		1.5 * (a8 * shape_funcs(7) - a7 * shape_funcs(6)),
		b8* shape_funcs(7) + b7 * shape_funcs(6),
		shape_funcs(3) - c8 * shape_funcs(7) - c7 * shape_funcs(6);

	HVecs.row(1) << 1.5 * (d5 * shape_funcs(4) - d8 * shape_funcs(7)),
		-shape_funcs(0) + e5 * shape_funcs(4) + e8 * shape_funcs(7),
		-b5 * shape_funcs(4) - b8 * shape_funcs(7),
		1.5 * (d6 * shape_funcs(5) - d5 * shape_funcs(4)),
		-shape_funcs(1) + e6 * shape_funcs(5) + e5 * shape_funcs(4),
		-b6 * shape_funcs(5) - b5 * shape_funcs(4),
		1.5 * (d7 * shape_funcs(6) - d6 * shape_funcs(5)),
		-shape_funcs(2) + e7 * shape_funcs(6) + e6 * shape_funcs(5),
		-b7 * shape_funcs(6) - b6 * shape_funcs(5),
		1.5 * (d8 * shape_funcs(7) - d7 * shape_funcs(6)),
		-shape_funcs(3) + e8 * shape_funcs(7) + e7 * shape_funcs(6),
		-b8 * shape_funcs(7) - b7 * shape_funcs(6);

	return HVecs;
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
	mat << 1, Mat->Poisson, 0,
		Mat->Poisson, 1, 0,
		0, 0, (1.0 - Mat->Poisson) / 2.0;
	return Mat->Young * pow(thickness, 3) / 12 / (1 - Mat->Poisson * Mat->Poisson) * mat;
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

QuadPlateElement::QuadPlateElement(Node* n0, Node* n1, Node* n2, Node* n3, double t, Material* mat)
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
	//plane = Plane::CreateFromPoints(n0->Location, n1->Location, n2->Location);
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

	return matrix;
}

Eigen::MatrixXd QuadPlateElement::StiffnessMatrix()
{
	Eigen::MatrixXd Kpln = plane_element.localStiffnessMatrix();
	Eigen::MatrixXd Kplt = localStiffnessMatrix();
	//std::cout << "K plane: \n" << Kpln << std::endl;
	//std::cout << "K plate: \n" << Kplt << std::endl;
	double max_diag = std::max(Kpln.diagonal().maxCoeff(), Kplt.diagonal().maxCoeff());

	Eigen::Matrix4d Krotz = Eigen::Matrix4d::Identity() * max_diag / 1000;

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
	return trMat.transpose() * mat * trMat;
}

void QuadPlateElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat)
{
	int indices[total_dof];
	for (size_t i = 0; i < node_num; i++)
		for (size_t j = 0; j < node_dof; j++)
			indices[node_dof * i + j] = Nodes[i]->id * 6 + j;

	Eigen::MatrixXd smat = StiffnessMatrix();
	for (size_t i = 0; i < total_dof; i++)
	{
		for (size_t j = 0; j < i + 1; j++)
		{
			int ci, rj;
			if (indices[j] <= indices[i]) {
				ci = indices[i];
				rj = indices[j];
			}
			else {
				ci = indices[j];
				rj = indices[i];
			}
			mat.coeffRef(rj, ci) += smat(i, j);
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

	return PlateStressData(strs(0), strs(1), strs(2), qx, qy);
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

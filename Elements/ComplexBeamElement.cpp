#include "ComplexBeamElement.h"

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

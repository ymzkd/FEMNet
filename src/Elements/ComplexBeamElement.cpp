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
//   Kb        : 4x4 の対称行列 (Eigen::Matrix4d)
//   lambda_s  : float (double)
//   lambda_z  : float (double)
//   lambda_s_ : float (double)  // λ_s'
//   lambda_z_ : float (double)  // λ_z'
//
// 戻り値
//   Kp        : 4x4 の対称行列 (Eigen::Matrix4d)
//
Eigen::Matrix4d compute_Kprime_partial(const Eigen::Matrix4d& Kb,
	double lambda_si, double lambda_zi,
	double lambda_sj, double lambda_zj)
{
	// --- 1) Kb の各要素を取り出す ---
	double k11 = Kb(0, 0);
	double k12 = Kb(0, 1);
	double k13 = Kb(0, 2);
	double k14 = Kb(0, 3);
	double k22 = Kb(1, 1);
	double k23 = Kb(1, 2);
	double k24 = Kb(1, 3);
	double k33 = Kb(2, 2);
	double k34 = Kb(2, 3);
	double k44 = Kb(3, 3);

	// --- 2) Lambda|D| の計算 ---
	double LambdaD =
		k11 * k22 * k33 * k44
		+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		- k11 * (k24 * k24) * k33 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		- k11 * (k23 * k23) * k44 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
		- k11 * k22 * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
		- (k12 * k12) * k33 * k44 * (1.0 - lambda_si) * (1.0 - lambda_zi)
		- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
		+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
		- (k13 * k13) * k22 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj)
		- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		- (k14 * k14) * k22 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zj)
		+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj);

	// --- 3) K' 用の 4x4 行列を 0 で初期化 ---
	Eigen::Matrix4d Kp = Eigen::Matrix4d::Zero();

	// ---------------------------------------------------------------------
	// 上三角 (i <= j) の要素を計算し，下三角へコピー (対称行列を構築)
	// ---------------------------------------------------------------------

	// === (A) K_{11}' ===
	Kp(0, 0) =
		((lambda_si * k11) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_zi)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_sj)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_zj)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			);

	// === (B) K_{12}' ===
	Kp(0, 1) =
		-(lambda_si * lambda_zi * k11 * k22 / LambdaD) * (
			-k12 * k33 * k44
			- k13 * k34 * k24 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- k14 * k23 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			+ k14 * k33 * k24 * (1.0 - lambda_zj)
			+ k13 * k23 * k44 * (1.0 - lambda_sj)
			+ k12 * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			);

	// === (C) K_{13}' ===
	Kp(0, 2) =
		-(lambda_si * lambda_sj * k11 * k33 / LambdaD) * (
			(1.0 - lambda_zi) * k12 * k23 * k44
			+ (1.0 - lambda_zi) * (1.0 - lambda_zj) * k13 * (k24 * k24)
			+ (1.0 - lambda_zj) * k14 * k22 * k34
			- (1.0 - lambda_zi) * (1.0 - lambda_zj) * k14 * k23 * k24
			- k13 * k22 * k44
			- (1.0 - lambda_zi) * (1.0 - lambda_zj) * k12 * k24 * k34
			);

	// === (D) K_{14}' ===
	Kp(0, 3) =
		-(lambda_si * lambda_zj * k11 * k44 / LambdaD) * (
			-(1.0 - lambda_sj) * (1.0 - lambda_zi) * k12 * k23 * k34
			- (1.0 - lambda_sj) * (1.0 - lambda_zi) * k13 * k24 * k23
			- k14 * k22 * k33
			+ (1.0 - lambda_sj) * (1.0 - lambda_zi) * k14 * (k23 * k23)
			+ (1.0 - lambda_sj) * k13 * k22 * k34
			+ (1.0 - lambda_zi) * k12 * k24 * k33
			);

	// === (E) K_{22}' ===
	Kp(1, 1) =
		((lambda_zi * k22) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_zj)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_sj)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_si)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zj)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			);

	// === (F) K_{23}' ===
	Kp(1, 2) =
		-(lambda_zi * lambda_sj * k22 * k33 / LambdaD) * (
			-k11 * k23 * k44
			- (1.0 - lambda_si) * (1.0 - lambda_zj) * k13 * k24 * k14
			- (1.0 - lambda_si) * (1.0 - lambda_zj) * k14 * k12 * k34
			+ (1.0 - lambda_si) * (1.0 - lambda_zj) * (k14 * k14) * k23
			+ (1.0 - lambda_si) * k13 * k12 * k44
			+ (1.0 - lambda_zj) * k11 * k24 * k34
			);

	// === (G) K_{24}' ===
	Kp(1, 3) =
		-(lambda_zi * lambda_zj * k22 * k44 / LambdaD) * (
			(1.0 - lambda_sj) * k11 * k23 * k34
			+ (1.0 - lambda_si) * (1.0 - lambda_sj) * (k13 * k13) * k24
			+ (1.0 - lambda_si) * k14 * k12 * k33
			- (1.0 - lambda_si) * (1.0 - lambda_sj) * k14 * k23 * k13
			- (1.0 - lambda_si) * (1.0 - lambda_sj) * k13 * k12 * k34
			- k11 * k24 * k33
			);

	// === (H) K_{33}' ===
	Kp(2, 2) =
		((lambda_sj * k33) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_zi)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_zj)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_si) * (1.0 - lambda_zi)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_si) * (1.0 - lambda_zi)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_si) * (1.0 - lambda_zj)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_si)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zj)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			);

	// === (I) K_{34}' ===
	Kp(2, 3) =
		-(lambda_sj * lambda_zj * k33 * k44 / LambdaD) * (
			-k11 * k22 * k34
			- (1.0 - lambda_si) * (1.0 - lambda_zi) * k12 * k24 * k13
			- (1.0 - lambda_si) * (1.0 - lambda_zi) * k14 * k12 * k23
			+ (1.0 - lambda_si) * k14 * k22 * k13
			+ (1.0 - lambda_si) * (1.0 - lambda_zi) * (k12 * k12) * k34
			+ (1.0 - lambda_zi) * k11 * k24 * k23
			);

	// === (J) K_{44}' ===
	Kp(3, 3) =
		((lambda_zj * k44) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_zi)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_sj)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_si) * (1.0 - lambda_zi)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zi)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_si)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
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

// compute_Displacement_TransformMatrix 関数
// バネ両端の変位を梁端部の変位に変換する行列を計算
//
// 引数
//   Kb        : 4x4 の対称行列 (Eigen::Matrix4d)
//   lambda_s  : float (double)
//   lambda_z  : float (double)
//   lambda_s_ : float (double)  // λ_s'
//   lambda_z_ : float (double)  // λ_z'
//
// 戻り値
//   Kp        : 4x4 の対称行列 (Eigen::Matrix4d)
//
Eigen::Matrix4d compute_Displacement_TransformMatrix(const Eigen::Matrix4d& Kb,
	double lambda_si, double lambda_zi,
	double lambda_sj, double lambda_zj)
{
	// --- 1) Kb の各要素を取り出す ---
	double k11 = Kb(0, 0);
	double k12 = Kb(0, 1);
	double k13 = Kb(0, 2);
	double k14 = Kb(0, 3);
	double k22 = Kb(1, 1);
	double k23 = Kb(1, 2);
	double k24 = Kb(1, 3);
	double k33 = Kb(2, 2);
	double k34 = Kb(2, 3);
	double k44 = Kb(3, 3);

	// --- 2) Lambda|D| の計算 ---
	double LambdaD =
		k11 * k22 * k33 * k44
		+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		- k11 * (k24 * k24) * k33 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		- k11 * (k23 * k23) * k44 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
		- k11 * k22 * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
		- (k12 * k12) * k33 * k44 * (1.0 - lambda_si) * (1.0 - lambda_zi)
		- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
		+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
		- (k13 * k13) * k22 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj)
		- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
		- (k14 * k14) * k22 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zj)
		+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj);

	// --- 3) K' 用の 4x4 行列を 0 で初期化 ---
	Eigen::Matrix4d Kp = Eigen::Matrix4d::Zero();

	// ---------------------------------------------------------------------
	// 上三角 (i <= j) の要素を計算し，下三角へコピー (対称行列を構築)
	// ---------------------------------------------------------------------

	// === K_{11}' ===
	Kp(0, 0) =
		1.0 - ((1.0 - lambda_si) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_zi)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_sj)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_zj)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_sj) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			);

	// === K_{12}' ===
	Kp(0, 1) =
		((1.0 - lambda_si) * lambda_zi * k22 / LambdaD) * (
			-k12 * k33 * k44
			- k13 * k34 * k24 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- k14 * k23 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			+ k14 * k33 * k24 * (1.0 - lambda_zj)
			+ k13 * k23 * k44 * (1.0 - lambda_sj)
			+ k12 * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			);

	// === K_{21}' ===
	Kp(1, 0) =
		(lambda_si * (1.0 - lambda_zi) * k11 / LambdaD) * (
			-k12 * k33 * k44
			- k13 * k34 * k24 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- k14 * k23 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			+ k14 * k33 * k24 * (1.0 - lambda_zj)
			+ k13 * k23 * k44 * (1.0 - lambda_sj)
			+ k12 * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			);

	// === K_{13}' ===
	Kp(0, 2) =
		((1.0 - lambda_si) * lambda_sj * k33 / LambdaD) * (
			(1.0 - lambda_zi) * k12 * k23 * k44
			+ (1.0 - lambda_zi) * (1.0 - lambda_zj) * k13 * (k24 * k24)
			+ (1.0 - lambda_zj) * k14 * k22 * k34
			- (1.0 - lambda_zi) * (1.0 - lambda_zj) * k14 * k23 * k24
			- k13 * k22 * k44
			- (1.0 - lambda_zi) * (1.0 - lambda_zj) * k12 * k24 * k34
			);
	
	// === K_{31}' ===
	Kp(2, 0) =
		(lambda_si * (1.0 - lambda_sj) * k11 / LambdaD) * (
			(1.0 - lambda_zi) * k12 * k23 * k44
			+ (1.0 - lambda_zi) * (1.0 - lambda_zj) * k13 * (k24 * k24)
			+ (1.0 - lambda_zj) * k14 * k22 * k34
			- (1.0 - lambda_zi) * (1.0 - lambda_zj) * k14 * k23 * k24
			- k13 * k22 * k44
			- (1.0 - lambda_zi) * (1.0 - lambda_zj) * k12 * k24 * k34
			);

	// === K_{14}' ===
	Kp(0, 3) =
		((1.0 - lambda_si) * lambda_zj * k44 / LambdaD) * (
			-(1.0 - lambda_sj) * (1.0 - lambda_zi) * k12 * k23 * k34
			- (1.0 - lambda_sj) * (1.0 - lambda_zi) * k13 * k24 * k23
			- k14 * k22 * k33
			+ (1.0 - lambda_sj) * (1.0 - lambda_zi) * k14 * (k23 * k23)
			+ (1.0 - lambda_sj) * k13 * k22 * k34
			+ (1.0 - lambda_zi) * k12 * k24 * k33
			);

	// === K_{41}' ===
	Kp(3, 0) =
		(lambda_si * (1.0 - lambda_zj) * k11 / LambdaD) * (
			-(1.0 - lambda_sj) * (1.0 - lambda_zi) * k12 * k23 * k34
			- (1.0 - lambda_sj) * (1.0 - lambda_zi) * k13 * k24 * k23
			- k14 * k22 * k33
			+ (1.0 - lambda_sj) * (1.0 - lambda_zi) * k14 * (k23 * k23)
			+ (1.0 - lambda_sj) * k13 * k22 * k34
			+ (1.0 - lambda_zi) * k12 * k24 * k33
			);

	// === K_{22}' ===
	Kp(1, 1) =
		1.0 - ((1.0 - lambda_zi) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_zj)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_sj)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_si)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zj)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zj)
			);

	// === K_{23}' ===
	Kp(1, 2) =
		((1.0 - lambda_zi) * lambda_sj * k33 / LambdaD) * (
			-k11 * k23 * k44
			- (1.0 - lambda_si) * (1.0 - lambda_zj) * k13 * k24 * k14
			- (1.0 - lambda_si) * (1.0 - lambda_zj) * k14 * k12 * k34
			+ (1.0 - lambda_si) * (1.0 - lambda_zj) * (k14 * k14) * k23
			+ (1.0 - lambda_si) * k13 * k12 * k44
			+ (1.0 - lambda_zj) * k11 * k24 * k34
			);

	// === K_{32}' ===
	Kp(2, 1) =
		(lambda_zi * (1.0 - lambda_sj) * k22 / LambdaD) * (
			-k11 * k23 * k44
			- (1.0 - lambda_si) * (1.0 - lambda_zj) * k13 * k24 * k14
			- (1.0 - lambda_si) * (1.0 - lambda_zj) * k14 * k12 * k34
			+ (1.0 - lambda_si) * (1.0 - lambda_zj) * (k14 * k14) * k23
			+ (1.0 - lambda_si) * k13 * k12 * k44
			+ (1.0 - lambda_zj) * k11 * k24 * k34
			);

	// === K_{24}' ===
	Kp(1, 3) =
		((1.0 - lambda_zi) * lambda_zj * k44 / LambdaD) * (
			(1.0 - lambda_sj) * k11 * k23 * k34
			+ (1.0 - lambda_si) * (1.0 - lambda_sj) * (k13 * k13) * k24
			+ (1.0 - lambda_si) * k14 * k12 * k33
			- (1.0 - lambda_si) * (1.0 - lambda_sj) * k14 * k23 * k13
			- (1.0 - lambda_si) * (1.0 - lambda_sj) * k13 * k12 * k34
			- k11 * k24 * k33
			);

	// === K_{42}' ===
	Kp(3, 1) =
		(lambda_zi * (1.0 - lambda_zj) * k22 / LambdaD) * (
			(1.0 - lambda_sj) * k11 * k23 * k34
			+ (1.0 - lambda_si) * (1.0 - lambda_sj) * (k13 * k13) * k24
			+ (1.0 - lambda_si) * k14 * k12 * k33
			- (1.0 - lambda_si) * (1.0 - lambda_sj) * k14 * k23 * k13
			- (1.0 - lambda_si) * (1.0 - lambda_sj) * k13 * k12 * k34
			- k11 * k24 * k33
			);

	// === K_{33}' ===
	Kp(2, 2) =
		1.0 - ((1.0 - lambda_sj) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_zi)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_zj)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_si) * (1.0 - lambda_zi)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_si) * (1.0 - lambda_zi)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_si) * (1.0 - lambda_zj)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_si)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zj)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_si) * (1.0 - lambda_zi) * (1.0 - lambda_zj)
			);

	// === K_{34}' ===
	Kp(2, 3) =
		((1.0 - lambda_sj) * lambda_zj * k44 / LambdaD) * (
			-k11 * k22 * k34
			- (1.0 - lambda_si) * (1.0 - lambda_zi) * k12 * k24 * k13
			- (1.0 - lambda_si) * (1.0 - lambda_zi) * k14 * k12 * k23
			+ (1.0 - lambda_si) * k14 * k22 * k13
			+ (1.0 - lambda_si) * (1.0 - lambda_zi) * (k12 * k12) * k34
			+ (1.0 - lambda_zi) * k11 * k24 * k23
			);

	// === K_{43}' ===
	Kp(3, 2) =
		(lambda_sj * (1.0 - lambda_zj) * k33 / LambdaD) * (
			-k11 * k22 * k34
			- (1.0 - lambda_si) * (1.0 - lambda_zi) * k12 * k24 * k13
			- (1.0 - lambda_si) * (1.0 - lambda_zi) * k14 * k12 * k23
			+ (1.0 - lambda_si) * k14 * k22 * k13
			+ (1.0 - lambda_si) * (1.0 - lambda_zi) * (k12 * k12) * k34
			+ (1.0 - lambda_zi) * k11 * k24 * k23
			);

	// === K_{44}' ===
	Kp(3, 3) =
		1.0 - ((1.0 - lambda_zj) / LambdaD) * (
			k11 * k22 * k33 * k44
			+ 2.0 * k11 * k23 * k24 * k34 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- k11 * (k24 * k24) * k33 * (1.0 - lambda_zi)
			- k11 * (k23 * k23) * k44 * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- k11 * k22 * (k34 * k34) * (1.0 - lambda_sj)
			- (k12 * k12) * k33 * k44 * (1.0 - lambda_si) * (1.0 - lambda_zi)
			- 2.0 * k12 * k13 * k24 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- 2.0 * k12 * k14 * k23 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ 2.0 * k12 * k14 * k24 * k33 * (1.0 - lambda_si) * (1.0 - lambda_zi)
			+ 2.0 * k12 * k13 * k23 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ (k12 * k12) * (k34 * k34) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ (k13 * k13) * (k24 * k24) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			+ 2.0 * k13 * k14 * k22 * k34 * (1.0 - lambda_si) * (1.0 - lambda_sj)
			- (k13 * k13) * k22 * k44 * (1.0 - lambda_si) * (1.0 - lambda_sj)
			- 2.0 * k13 * k14 * k23 * k24 * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			- (k14 * k14) * k22 * k33 * (1.0 - lambda_si)
			+ (k14 * k14) * (k23 * k23) * (1.0 - lambda_si) * (1.0 - lambda_sj) * (1.0 - lambda_zi)
			);

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
    Kbz = compute_Kprime_partial(Kbz, Lambda_syi, Lambda_bzi, Lambda_syj, Lambda_bzj);
    Kby = compute_Kprime_partial(Kby, Lambda_szi, Lambda_byi, Lambda_szj, Lambda_byj);

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
    Lambda_bzi = 1;
    Lambda_bzj = 1;
    Lambda_syi = 1;
    Lambda_syj = 1;
    Lambda_byi = 1;
    Lambda_byj = 1;
    Lambda_szi = 1;
    Lambda_szj = 1;

    lzi = 0;
    lzj = 0;
    lyi = 0;
    lyj = 0;
}

ComplexBeamElement::ComplexBeamElement(Node *n0, Node *n1, Section *sec, Material mat, double beta)
    : BeamElement(n0, n1, sec, mat, beta)
{
    Lambda_bzi = 1;
    Lambda_bzj = 1;
    Lambda_syi = 1;
    Lambda_syj = 1;
    Lambda_byi = 1;
    Lambda_byj = 1;
    Lambda_szi = 1;
    Lambda_szj = 1;

    lzi = 0;
    lzj = 0;
    lyi = 0;
    lyj = 0;
}

ComplexBeamElement::ComplexBeamElement(int _id, Node *n0, Node *n1, Section *sec, Material mat, double beta)
    : ComplexBeamElement(n0, n1, sec, mat, beta)
{
    id = _id;

    Lambda_bzi = 1;
    Lambda_bzj = 1;
    Lambda_syi = 1;
    Lambda_syj = 1;
    Lambda_byi = 1;
    Lambda_byj = 1;
    Lambda_szi = 1;
    Lambda_szj = 1;

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

Displacement ComplexBeamElement::DisplaceAt(Displacement d0, Displacement d1, double p)
{
	Eigen::VectorXd wvec(12);
	wvec.segment(0, 6) = Eigen::Map<Eigen::VectorXd>(d0.displace, 6);
	wvec.segment(6, 6) = Eigen::Map<Eigen::VectorXd>(d1.displace, 6);
	
	double p2 = p * p;
	double p3 = p * p * p;
	double l = length();
	double lz = l - lzi - lzj;
	double ly = l - lyi - lyj;

	Eigen::MatrixXd T = trans_matrix();
	wvec = T * wvec;

	// z軸周り, y方向たわみ
	Eigen::Vector4d vvec_ij(wvec(1), wvec(5), wvec(7), wvec(11));
	Eigen::Matrix4d Tbz = Eigen::Matrix4d::Identity();
	Tbz(0, 1) = lzi; Tbz(2, 3) = -lzj;
	Eigen::Vector4d vvec_mq = Tbz * vvec_ij;
	Eigen::Matrix4d Kbz = stiffness_matrix_beam_rot_z(Mat.Young, Sec->Iz, lz);
	Kbz = compute_Displacement_TransformMatrix(Kbz, Lambda_syi, Lambda_bzi, Lambda_syj, Lambda_bzj);
	Eigen::Vector4d vvec_np = Kbz * vvec_mq;
	Eigen::Vector4d Nrz_v(1.0 - 3.0 * p2 + 2.0 * p3, p * lz - 2.0 * p2 * lz + p3 * lz, 3.0 * p2 - 2.0 * p3, -p2 * lz + p3 * lz);
	Eigen::Vector4d Nrz_sz(-6.0 * p / lz + 6.0 * p2 / lz, 1.0 - 4.0 * p + 3.0 * p2, 6.0 * p / lz - 6.0 * p2 / lz, -2.0 * p + 3.0 * p2);
	double v = Nrz_v.dot(vvec_np);
	double sz = Nrz_sz.dot(vvec_np);

	// y軸周り, z方向たわみ
	Eigen::Vector4d wvec_ij(wvec(2), wvec(4), wvec(8), wvec(10));
	Eigen::Matrix4d Tby = Eigen::Matrix4d::Identity();
	Tby(0, 1) = -lyi; Tby(2, 3) = lyj;
	Eigen::Vector4d wvec_mq = Tby * wvec_ij;
	Eigen::Matrix4d Kby = stiffness_matrix_beam_rot_y(Mat.Young, Sec->Iy, ly);
	Kby = compute_Displacement_TransformMatrix(Kby, Lambda_szi, Lambda_byi, Lambda_szj, Lambda_byj);
	Eigen::Vector4d wvec_np = Kby * wvec_mq;
	Eigen::Vector4d Nry_w(1.0 - 3.0 * p2 + 2.0 * p3, -p * ly + 2.0 * p2 * ly - p3 * ly, 3.0 * p2 - 2.0 * p3, p2 * ly - p3 * ly);
	Eigen::Vector4d Nry_sy(6.0 * p / ly - 6.0 * p2 / ly, 1.0 - 4.0 * p + 3.0 * p2, -6.0 * p / ly + 6.0 * p2 / ly, -2.0 * p + 3.0 * p2);

	double w = Nry_w.dot(wvec_np);
	double sy = Nry_sy.dot(wvec_np);

	// 軸方向, ねじり
	double u = wvec(0) * (1.0 - p) + wvec(6) * p;
	double sx = wvec(3) * (1.0 - p) + wvec(9) * p;

	Eigen::Vector<double, 6> ldisp;
	ldisp << u, v, w, sx, sy, sz;
	Eigen::VectorXd gdisp = T.block(0, 0, 6, 6).transpose() * ldisp;
	
	return Displacement(gdisp(0), gdisp(1), gdisp(2), gdisp(3), gdisp(4), gdisp(5));
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

double ComplexBeamElement::get_ksyi()
{
	double lz = length() - -lzi - lzj;
	return Lambda_syi / (1.0 - Lambda_syi) * Sec->Iz * Mat.Young * 12.0 / (lz * lz * lz);
}

double ComplexBeamElement::get_kbzi()
{
	double lz = length() - -lzi - lzj;
	return Lambda_bzi / (1.0 - Lambda_bzi) * Sec->Iz * Mat.Young * 4.0 / lz;
}

double ComplexBeamElement::get_ksyj()
{
	double lz = length() - -lzi - lzj;
	return Lambda_syj / (1.0 - Lambda_syj) * Sec->Iz * Mat.Young * 12.0 / (lz * lz * lz);
}

double ComplexBeamElement::get_kbzj()
{
	double lz = length() - -lzi - lzj;
	return Lambda_bzj / (1.0 - Lambda_bzj) * Sec->Iz * Mat.Young * 4.0 / lz;
}




double ComplexBeamElement::get_kszi()
{
	double ly = length() - -lyi - lyj;
	return Lambda_szi / (1.0 - Lambda_szi) * Sec->Iy * Mat.Young * 12.0 / (ly * ly * ly);
}

double ComplexBeamElement::get_kbyi()
{
	double ly = length() - -lyi - lyj;
	return Lambda_byi / (1.0 - Lambda_byi) * Sec->Iy * Mat.Young * 4.0 / ly;
}

double ComplexBeamElement::get_kszj()
{
	double ly = length() - -lyi - lyj;
	return Lambda_szj / (1.0 - Lambda_szj) * Sec->Iy * Mat.Young * 12.0 / (ly * ly * ly);
}

double ComplexBeamElement::get_kbyj()
{
	double ly = length() - -lyi - lyj;
	return Lambda_byj / (1.0 - Lambda_byj) * Sec->Iy * Mat.Young * 4.0 / ly;
}



void ComplexBeamElement::set_ksyi(double ksyi)
{
	double lz = length() - -lzi - lzj;
	Lambda_syi = ksyi / (Sec->Iz * Mat.Young * 12.0 / (lz * lz * lz) + ksyi);
}

void ComplexBeamElement::set_kbzi(double kbzi)
{
	double lz = length() - -lzi - lzj;
	Lambda_bzi = kbzi / (Sec->Iz * Mat.Young * 4.0 / lz + kbzi);
}

void ComplexBeamElement::set_ksyj(double ksyj)
{
	double lz = length() - -lzi - lzj;
	Lambda_syj = ksyj / (Sec->Iz * Mat.Young * 12.0 / (lz * lz * lz) + ksyj);
}

void ComplexBeamElement::set_kbzj(double kbzj)
{
	double lz = length() - -lzi - lzj;
	Lambda_bzj = kbzj / (Sec->Iz * Mat.Young * 4.0 / lz + kbzj);
}


void ComplexBeamElement::set_kszi(double kszi)
{
	double ly = length() - -lyi - lyj;
	Lambda_szi = kszi / (Sec->Iy * Mat.Young * 12.0 / (ly * ly * ly) + kszi);
}

void ComplexBeamElement::set_kbyi(double kbyi)
{
	double ly = length() - -lyi - lyj;
	Lambda_byi = kbyi / (Sec->Iy * Mat.Young * 4.0 / ly + kbyi);
}

void ComplexBeamElement::set_kszj(double kszj)
{
	double ly = length() - -lyi - lyj;
	Lambda_szj = kszj / (Sec->Iy * Mat.Young * 12.0 / (ly * ly * ly) + kszj);
}

void ComplexBeamElement::set_kbyj(double kbyj)
{
	double ly = length() - -lyi - lyj;
	Lambda_byj = kbyj / (Sec->Iy * Mat.Young * 4.0 / ly + kbyj);
}

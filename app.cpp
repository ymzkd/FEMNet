#include<iostream>

//#ifdef USE_MKL
//#define EIGEN_USE_MKL_ALL
//#endif

#include "Element.h"
#include "Model.h"
#include "LoadComponent.h"

void TestMethod1() {

	Point p0(0, 0, 0);
	Point p1(1, 0, 0);
	Point p2(0, 1, 0);
	Plane pl = Plane::CreateFromPoints(p0, p1, p2);

	Point pt(0.5, 0.5, 0.5);
	std::cout << "ex: " << pl.ex << std::endl;
	std::cout << "ey: " << pl.ey << std::endl;
	std::cout << "ez: " << pl.ez << std::endl;
	std::cout << "pt: " << pl.PointToCoord(pt) << std::endl;

	Material m0(5000, 0.2);
	Node n0(p0);
	Node n1(p1);
	Node n2(p2);
	TriPlaneElement pel(&n0, &n1, &n2, 12.0, m0);

	// std::cout << "TransMat: \n" << pel.trans_matrix() << std::endl;
	std::cout << "StiffMat: \n" << pel.StiffnessMatrix() << std::endl;
}

void TestMethod2() {

	std::cout << "\nTestMethod2 Start" << std::endl;

	double h_compare[] = { 0.0, 0.1, 1.0 };

	std::cout << "h\tNode2 dx\tNode2 dy\tNode2 dz" << std::endl;

	for each(double h in h_compare)
	{
		Point p0(0, 0, h);
		Point p1(10, 0, -h);
		Point p2(10, 10, h);
		Point p3(0, 10, -h);

		Material m0(1e6, 0.3);
		Node n0(p0); n0.id = 0;
		Node n1(p1); n1.id = 1;
		Node n2(p2); n2.id = 2;
		Node n3(p3); n3.id = 3;
		QuadPlateElement pel(&n0, &n1, &n2, &n3, 0.1, m0);

		FEModel model;
		model.Nodes.push_back(n0);
		model.Nodes.push_back(n1);
		model.Nodes.push_back(n2);
		model.Nodes.push_back(n3);

		model.Nodes[0].Fix.FixAll();
		model.Nodes[1].Fix.FixAll();

		model.add_element(pel);

		NodeLoad nl1(2, 0, 0.5, 0);
		NodeLoad nl2(3, 0, 0.5, 0);
		std::vector<std::shared_ptr<LoadBase>> loads;
		loads.push_back(std::make_shared<NodeLoad>(nl1));
		loads.push_back(std::make_shared<NodeLoad>(nl2));

		std::vector<Displacement> disp;
		std::vector<NodeLoad> react;
		model.SolveLinearStatic(loads, disp, react);

		std::cout << h << "\t" << disp[2].Dx() << "\t" << disp[2].Dy() << "\t" << disp[2].Dz() << std::endl;
	}
}

// #include <Eigen/Sparse>
// #include <Eigen/SparseCholesky>

// void CheckSparseSolver() {
	
// 	using namespace Eigen;

// 	// 型定義
// 	typedef SparseMatrix<double> SpMat;
// 	typedef Triplet<double> T;

// 	// 3x3の対称正定値疎行列を作成（上三角にデータを格納）
// 	std::vector<T> triplets;
// 	triplets.emplace_back(0, 0, 4.0);
// 	triplets.emplace_back(0, 1, 1.0); // 上三角だけ指定
// 	triplets.emplace_back(1, 1, 3.0);
// 	triplets.emplace_back(1, 2, 2.0);
// 	triplets.emplace_back(2, 2, 5.0);

// 	SpMat A(3, 3);
// 	A.setFromTriplets(triplets.begin(), triplets.end());

// 	// 右辺ベクトル b
// 	VectorXd b(3);
// 	b << 1.0, 2.0, 3.0;

// 	// LLT分解（上三角を使う）
// #ifdef EIGEN_USE_MKL_ALL
// 	Eigen::PardisoLLT<Eigen::SparseMatrix<double>> solver;
// #else
// 	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
// #endif
// 	solver.compute(A);

// 	if (solver.info() != Success) {
// 		std::cerr << "分解に失敗しました" << std::endl;
// 	}

// 	// 解 x を求める
// 	VectorXd x = solver.solve(b);

// 	if (solver.info() != Success) {
// 		std::cerr << "求解に失敗しました" << std::endl;
// 	}

// 	std::cout << "解 x:\n" << x << std::endl;


// 	// 解 x を求める
// 	VectorXd x2 = solver.solve(b);

// 	if (solver.info() != Success) {
// 		std::cerr << "求解に失敗しました" << std::endl;
// 	}

// 	std::cout << "解2 x:\n" << x2 << std::endl;
// }

void CheckQuadElement1() {

	Point p0(0, 0, 0);
	Point p1(1, 0, 0);
	Point p2(1, 1, 0);
	Point p3(0, 1, 0);
	// Plane pl = Plane::CreateFromPoints(p0, p1, p2);

	Material m0(5000, 0.2);
	m0.dense = 4.3;
	Node n0(p0);
	Node n1(p1);
	Node n2(p2);
	Node n3(p3);
	QuadPlateElement pel(&n0, &n1, &n2, &n3, 0.5, m0);

	std::cout << "Area: \n" << pel.Area() << std::endl;
	std::cout << "Mass: \n" << pel.NodeLumpedMass() << std::endl;

}

FEModel CantiBeamModel(double l, int n) {
	FEModel model;
	double dl = l / n;
	model.Nodes.push_back(Node(0, 0, 0, 0));
	model.Nodes[0].Fix.FixAll();
	for (size_t i = 0; i < n; i++)
		model.Nodes.push_back(Node(i+1, dl*(i+1), 0, 0));

	Material m0(5000.0, 0.2);
	m0.dense = 5.0 / 1000.0 / 1000.0;

	Section s0(100, 833.33, 833.33, 1406.25);

	model.Materials.push_back(m0);
	model.Sections.push_back(s0);

	for (size_t i = 0; i < n; i++) {
		std::shared_ptr<BeamElement> b1 = 
			std::make_shared<BeamElement>(i, &model.Nodes[i], &model.Nodes[i + 1], 
				&model.Sections[0], model.Materials[0]);
		model.Elements.push_back(b1);
		//model.add_element(BeamElement(i, &model.Nodes[i], &model.Nodes[i + 1], &s0, m0));
	}

	return model;
}

FEModel CantiColumnModel(double l, int n) {
	FEModel model;
	double dl = l / (double)n;
	model.Nodes.push_back(Node(0, 0, 0, 0));
	model.Nodes[0].Fix.FixAll();
	for (size_t i = 0; i < n; i++)
		model.Nodes.push_back(Node(i + 1, 0, 0, dl * (i + 1)));

	Material m0(5000.0, 0.2);
	m0.dense = 5.0 / 1000.0 / 1000.0;

	Section s0(100, 833.33, 833.33, 1406.25);

	model.Materials.push_back(m0);
	model.Sections.push_back(s0);

	for (size_t i = 0; i < n; i++) {
		std::shared_ptr<BeamElement> b1 =
			std::make_shared<BeamElement>(i, &model.Nodes[i], &model.Nodes[i + 1],
				&model.Sections[0], model.Materials[0]);
		model.Elements.push_back(b1);
		//model.add_element(BeamElement(i, &model.Nodes[i], &model.Nodes[i + 1], &s0, m0));
	}

	return model;
}


void CheckCantiBeamVibration() {

	std::cout << "CheckCantiBeamVibration Start" << std::endl;

	Material m0(5000, 0.2);
	m0.dense = 4.3;
	Node n0(0, 0, 0); n0.id = 0;
	Node n1(100, 0, 0); n1.id = 1;
	Node n2(200, 0, 0); n2.id = 2;
	n0.Fix.FixAll();

	Section s0(100, 833.33, 833.33, 1406.25);

	// NodeLoad nl = NodeLoad(2, 0, 0, -10, 0, 0, 0);
	// std::vector<std::shared_ptr<LoadBase>> loads;
	// loads.push_back(std::make_shared<NodeLoad>(nl));

	// std::vector<Displacement> disp;
	// std::vector<NodeLoad> react;
	

	FEModel model;

	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);

	model.Materials.push_back(m0);
	model.Sections.push_back(s0);

	//BeamElement b0(0, &n0, &n1, &s0, m0);
	//BeamElement b1(1, &n1, &n2, &s0, m0);
	//model.add_element(b0);
	//model.add_element(b1);
	std::shared_ptr<BeamElement> b1 = std::make_shared<BeamElement>(0, &n0, &n1, &s0, m0);
	std::shared_ptr<BeamElement> b2 = std::make_shared<BeamElement>(1, &n1, &n2, &s0, m0);
	model.Elements.push_back(b1);
	model.Elements.push_back(b2);



	// model.Solve(loads, disp, react);
	// std::cout << "Displacement: " << disp[2] << std::endl;
	//model.SolveVibrationTest();
	//model.SolveVibrationTest2();

}

void TestQuadPlateConsistMass() {

	std::cout << std::endl << "TestQuadPlateConsistMass Start" << std::endl;

	// ノードの作成
	//Node n0(0, 0, 0);
	//Node n1(1, 0, 0);
	//Node n2(1, 1, 0);
	//Node n3(0, 1, 0);
	Node n0(0.0, 0.0, 0.0);
	Node n1(1.05, -0.03, 0.0);
	Node n2(0.98, 1.02, 0.0);
	Node n3(-0.04, 0.97, 0.0);

	// 材料と厚さの設定
	Material mat(5000, 0.2, 0.5);
	Thickness thickness;
	thickness.plate_thick = 0.1;
	thickness.weight_thick = 0.1;
	thickness.plane_thick = 0.1;

	// QuadPlateElementのインスタンスを作成
	QuadPlateElement element(&n0, &n1, &n2, &n3, thickness, mat);

	// NodeConsistentMassメソッドを呼び出して結果を取得
	Eigen::MatrixXd massMatrix = element.NodeConsistentMass();

	// 結果を出力
	std::cout << "Node Consistent Mass Matrix:" << std::endl;
	std::cout << massMatrix.diagonal() << std::endl;

	//// NodeConsistentMassメソッドを呼び出して結果を取得
	//massMatrix = element.NodeConsistentMass2();

	//// 結果を出力
	//std::cout << "Node Consistent Mass Matrix2:" << std::endl;
	//std::cout << massMatrix.diagonal() << std::endl;

}

void TestTriPlateConsistentMass() {

	std::cout << std::endl << "TestTriPlateConsistentMass Start" << std::endl;
	// ノードの作成
	Node n0(0.0, 0.0, 0.0);
	Node n1(1.0, 0.0, 0.0);
	Node n2(0.5, 1.0, 0.0);

	// 材料と厚さの設定
	Material mat(5000, 0.3, 7.9); // ヤング率: 5000, ポアソン比: 0.3
	// mat.dense = 7.85;        // 材料密度 (例: 鉄の密度)

	Thickness thickness(0.1); // 板厚 0.1

	// TriPlateElement のインスタンスを作成
	TriPlateElement element(&n0, &n1, &n2, thickness, mat);

	// NodeConsistentMass メソッドを呼び出して結果を取得
	Eigen::MatrixXd massMatrix = element.NodeConsistentMass();

	// 結果を出力
	std::cout << "Node Consistent Mass Matrix:" << std::endl;
	std::cout << massMatrix.diagonal() << std::endl;
}

void CheckCantiBeamVibration2() {

	std::cout << "CheckCantiBeamVibration2 Start" << std::endl;
	FEModel model = CantiBeamModel(200, 10);

	// model.Solve(loads, disp, react);
	// std::cout << "Displacement: " << disp[2] << std::endl;
	std::vector<double> eigen_values;
	std::vector<std::vector<Displacement>> mode_vectors;
	model.SolveVibration(6, eigen_values, mode_vectors);
	//model.SolveVibrationTest2();
	double pi = 3.141592653589793238462643;
	std::cout << "Natural Periods" << std::endl;
	for each (double var in eigen_values)
	{
		std::cout << 2 * pi / var << "\n";
	}
}

// 座屈検討用の片持ち柱
void CheckCantiBeamBuckling() {

	std::cout << "CheckCantiBeamBuckling Start" << std::endl;
	FEModel model = CantiColumnModel(2000, 10);
	//FEModel model = CantiBeamModel(2000, 10);
	NodeLoad nl = NodeLoad(model.Nodes.size() - 1, 0, 0, -1.0);
	//NodeLoad nl = NodeLoad(model.Nodes.size() - 1, -1.0, 0, 0);
	std::vector<std::shared_ptr<LoadBase>> loads;
	loads.push_back(std::make_shared<NodeLoad>(nl));

	std::vector<Displacement> disp;
	std::vector<NodeLoad> react;
	model.Solve(loads, disp, react);

	FEStaticResult result = FEStaticResult(std::make_shared<FEModel>(model), loads, disp, react);

	std::cout << "Static Result: " << result.displace[model.Nodes.size() - 1] << std::endl;
	std::cout << "Beam Stress at Node 0 At 0: " << result.GetBeamStress(0, 0) << std::endl;

	FEBucklingAnalysis buckling_analysis(std::make_shared<FEStaticResult>(result));
	buckling_analysis.mode_num = 6;
	buckling_analysis.SolveBuckling();

	std::cout << "Buckling Analysis Results:" << std::endl;
	std::cout << "Mode Number: " << buckling_analysis.mode_num << std::endl;
	for (size_t i = 0; i < buckling_analysis.eigs.size(); i++)
	{
		std::cout << "Mode " << i << ": " << buckling_analysis.eigs[i] << std::endl;
	}
}

FEModel PyramidTrussModel(double D, int n, double h)
{
	FEModel model;
	double dl = D / (double)n;
	model.Nodes.push_back(Node(0, 0, 0, h));
	// pi / n
	double theta = 3.14159265358979323846 * 2.0 / n;
	for (size_t i = 0; i < n; i++)
	{
		double x = D * cos(theta * i);
		double y = D * sin(theta * i);
		Node ni(i + 1, x, y, 0);
		ni.Fix.FixAll();
		model.Nodes.push_back(ni);

	}

	Material m0(5000.0, 0.2);
	m0.dense = 5.0 / 1000.0 / 1000.0;

	Section s0(100, 833.33, 833.33, 1406.25);

	model.Materials.push_back(m0);
	model.Sections.push_back(s0);

	for (size_t i = 0; i < n; i++)
	{
		std::shared_ptr<TrussElement> t1 =
			std::make_shared<TrussElement>(i, &model.Nodes[i + 1], &model.Nodes[0],
				&model.Sections[0], model.Materials[0]);
		model.Elements.push_back(t1);
		// model.add_element(BeamElement(i, &model.Nodes[i], &model.Nodes[i + 1], &s0, m0));
	}

	return model;
}

// 座屈検討用のピラミッド型トラスサンプル
void CheckCantiPyramidTrussBuckling(double D, int n, double h)
{

	std::cout << "CheckCantiPyramidTrussBuckling Start" << std::endl;
	FEModel model = PyramidTrussModel(1000, 4, 100);
	// FEModel model = CantiBeamModel(2000, 10);
	NodeLoad nl = NodeLoad(0, 0, 0, -1.0);

	std::vector<std::shared_ptr<LoadBase>> loads;
	loads.push_back(std::make_shared<NodeLoad>(nl));

	std::vector<Displacement> disp;
	std::vector<NodeLoad> react;
	model.Solve(loads, disp, react);

	FEStaticResult result = FEStaticResult(std::make_shared<FEModel>(model), loads, disp, react);

	std::cout << "Static Result: " << result.displace[model.Nodes.size() - 1] << std::endl;
	std::cout << "Beam Stress at Node 0 At 0: " << result.GetBeamStress(0, 0) << std::endl;

	FEBucklingAnalysis buckling_analysis(std::make_shared<FEStaticResult>(result));
	buckling_analysis.mode_num = 1;
	buckling_analysis.SolveBuckling();

	std::cout << "Buckling Analysis Results:" << std::endl;
	std::cout << "Mode Number: " << buckling_analysis.mode_num << std::endl;
	for (size_t i = 0; i < buckling_analysis.eigs.size(); i++)
	{
		std::cout << "Mode " << i << ": " << buckling_analysis.eigs[i] << std::endl;
	}
}

void CheckQuadPlateBuckling() {
	std::cout << "CheckQuadPlateBuckling Start" << std::endl;
	// ノードの作成
	Node n0(0, 0.0, 0.0, 0.0);
	Node n1(1, 100.0, 0.0, 10.0);
	Node n2(2, 200.0, 0.0, 0.0);
	Node n3(3, 0.0, 100.0, 0.0);
	Node n4(4, 100.0, 100.0, 10.0);
	Node n5(5, 200.0, 100.0, 0.0);
	n0.Fix.PinFix();
	n2.Fix.PinFix();
	n3.Fix.PinFix();
	n5.Fix.PinFix();

	// 材料と厚さの設定
	Material mat(5000, 0.2);
	mat.dense = 4.3;
	Thickness thickness(2.0);
	// QuadPlateElementのインスタンスを作成
	QuadPlateElement element1(&n0, &n1, &n4, &n3, thickness, mat);
	QuadPlateElement element2(&n1, &n2, &n5, &n4, thickness, mat);
	// FEModelの作成
	FEModel model;
	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);
	model.Nodes.push_back(n3);
	model.Nodes.push_back(n4);
	model.Nodes.push_back(n5);
	model.Materials.push_back(mat);
	
	model.add_element(element1);
	model.add_element(element2);

	//std::shared_ptr<QuadPlateElement> qpe = std::make_shared<QuadPlateElement>(element1);
	//std::shared_ptr<QuadPlateElement> qpe2 = std::make_shared<QuadPlateElement>(element2);
	//model.Elements.push_back(qpe);
	//model.Elements.push_back(qpe2);

	NodeLoad nl1 = NodeLoad(1, 0, 0, -100.0);
	NodeLoad nl2 = NodeLoad(4, 0, 0, -100.0);
	std::vector<std::shared_ptr<LoadBase>> loads;
	loads.push_back(std::make_shared<NodeLoad>(nl1));
	loads.push_back(std::make_shared<NodeLoad>(nl2));
	std::vector<Displacement> disp;
	std::vector<NodeLoad> react;
	model.Solve(loads, disp, react);
	FEStaticResult result = FEStaticResult(std::make_shared<FEModel>(model), loads, disp, react);
	std::cout << "Static Result: " << result.displace[model.Nodes.size() - 1] << std::endl;
	FEBucklingAnalysis buckling_analysis(std::make_shared<FEStaticResult>(result));
	buckling_analysis.mode_num = 4;
	buckling_analysis.SolveBuckling();
	std::cout << "Buckling Analysis Results:" << std::endl;
	std::cout << "Mode Number: " << buckling_analysis.mode_num << std::endl;
	for (size_t i = 0; i < buckling_analysis.eigs.size(); i++)
	{
		std::cout << "Mode " << i << ": " << buckling_analysis.eigs[i] << std::endl;
	}
}

void TestBodyforceToNodeLoadData() {

	std::cout << "TestBodyforceToNodeLoadData Start" << std::endl;

	// サンプルの加速度ベクトルを定義
	Eigen::Vector3d accel_vec(0, 0.0, -9.81); // 重力加速度を x 軸方向に設定

	// サンプルのノードを作成
	Node node1(0, 0.0, 0.0, 0.0);
	Node node2(1, 500.0, 0.0, 0.0);

	// サンプルのセクションとマテリアルを作成
	Section section(100, 833.33, 833.33, 1406.25);
	Material material(5000, 0.2, 5.0);

	// サンプルの BeamElement を作成
	BeamElement beam(&node1, &node2, &section, material, 0.0);

	// BodyforceToNodeLoadData を呼び出し
	std::vector<NodeLoadData> loads = beam.InertialForceToNodeLoadData(accel_vec);

	// 結果を出力
	std::cout << "Node Load Data Results:" << std::endl;
	for (const auto& load : loads) {
		std::cout << "Node ID: " << load.id
			<< ", Px: " << load.Px()
			<< ", Py: " << load.Py()
			<< ", Pz: " << load.Pz()
			<< ", Mx: " << load.Mx()
			<< ", My: " << load.My()
			<< ", Mz: " << load.Mz()
			<< std::endl;
	}
}

void TestBeamInertialForceSolve() {

	std::cout << "TestBeamInertialForceSolve Start" << std::endl;

	// ノードの作成
	Node n0(0, 0, 0); n0.id = 0;
	Node n1(100, 0, 0); n1.id = 1;
	Node n2(200, 0, 0); n2.id = 2;
	n0.Fix.FixAll();

	// 材料と断面の設定
	Material mat(5000, 0.2, 7.85); // ヤング率: 5000, ポアソン比: 0.2, 密度: 7.85
	Section sec(100, 833.33, 833.33, 1406.25); // 断面積: 100, 慣性モーメント: 833.33, 1406.25

	// ビーム要素の作成
	BeamElement beam1(&n0, &n1, &sec, mat, 0.0);
	BeamElement beam2(&n1, &n2, &sec, mat, 0.0);

	// モデルの作成
	FEModel model;
	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);

	// 材料と断面をモデルに追加
	model.Materials.push_back(mat);
	model.Sections.push_back(sec);

	// ビーム要素をモデルに追加
	model.add_element(beam1);
	model.add_element(beam2);

	// 加速度ベクトル (重力加速度を z 軸方向に設定)
	Eigen::Vector3d accel_vec(0.0, 0.0, -9.81);

	// InertialForce を作成
	InertialForce inertial_force(accel_vec.x(), accel_vec.y(), accel_vec.z());

	//	std::vector<std::shared_ptr<LoadBase>> loads;
	//loads.push_back(std::make_shared<NodeLoad>(nl));

	// 荷重リストに追加
	std::vector<std::shared_ptr<LoadBase>> loads;
	loads.push_back(std::make_shared<InertialForce>(inertial_force));

	// 解析結果を格納する変数
	std::vector<Displacement> displacements;
	std::vector<NodeLoad> reactions;

	// 静的解析を実行
	model.SolveLinearStatic(loads, displacements, reactions);
	FEStaticResult result = FEStaticResult(std::make_shared<FEModel>(model), loads, displacements, reactions);
	
	// 結果を出力
	std::cout << "Displacements:" << std::endl;
	for (size_t i = 0; i < displacements.size(); ++i) {
		std::cout << "Node " << i << ": " << displacements[i] << std::endl;
	}

	std::cout << "Reactions:" << std::endl;
	for (const auto& reaction : reactions) {
		std::cout << reaction << std::endl;
	}

	std::cout << "Beam Stress at Node 0 At 0: " << result.GetBeamStress(0, 0) << std::endl;
	std::cout << "Beam Stress at Node 0 At 0.5: " << result.GetBeamStress(0, 0.5) << std::endl;
	std::cout << "Beam Stress at Node 0 At 1.0: " << result.GetBeamStress(0, 1.0) << std::endl;
	
}


// PlateLoadクラスのテスト
void TestPlatePressure() {
	
	std::cout << "TestPlatePressure Start" << std::endl;

	Node n0(0, 0, 0); n0.id = 0;
	Node n1(100, 0, 0); n1.id = 1;
	Node n2(200, 0, 0); n2.id = 2;
	Node n3(0, 100, 0); n3.id = 3;
	Node n4(100, 100, 0); n4.id = 4;
	Node n5(200, 100, 0); n5.id = 5;
	Node n6(0, 200, 0); n6.id = 6;
	Node n7(100, 200, 0); n7.id = 7;
	Node n8(200, 200, 0); n8.id = 8;
	n0.Fix.FixAll();
	n2.Fix.FixAll();
	n6.Fix.FixAll();
	n8.Fix.FixAll();

	Material m0(5000, 0.2);
	Thickness thickness(4.5);

	//QuadPlateElement element(&n0, &n1, &n2, &n3, thickness, mat);
	std::shared_ptr<QuadPlateElement> e0 = std::make_shared<QuadPlateElement>(0, &n0, &n1, &n4, &n3, thickness, m0);
	std::shared_ptr<QuadPlateElement> e1 = std::make_shared<QuadPlateElement>(1, &n1, &n2, &n5, &n4, thickness, m0);
	std::shared_ptr<QuadPlateElement> e2 = std::make_shared<QuadPlateElement>(2, &n3, &n4, &n7, &n6, thickness, m0);
	std::shared_ptr<QuadPlateElement> e3 = std::make_shared<QuadPlateElement>(3, &n4, &n5, &n8, &n7, thickness, m0);

	FEModel model;
	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);
	model.Nodes.push_back(n3);
	model.Nodes.push_back(n4);
	model.Nodes.push_back(n5);
	model.Nodes.push_back(n6);
	model.Nodes.push_back(n7);
	model.Nodes.push_back(n8);

	model.Materials.push_back(m0);

	model.Elements.push_back(e0);
	model.Elements.push_back(e1);
	model.Elements.push_back(e2);
	model.Elements.push_back(e3);

	double px = 0.0, py = 0.0, pz = -1.2;
	std::vector<std::shared_ptr<LoadBase>> loads;
	//PlateLoad pl0(e0.get(), px, py, pz);

	loads.push_back(std::make_shared<PlateLoad>(e0.get(), px, py, pz));
	loads.push_back(std::make_shared<PlateLoad>(e1.get(), px, py, pz));
	loads.push_back(std::make_shared<PlateLoad>(e2.get(), px, py, pz));
	loads.push_back(std::make_shared<PlateLoad>(e3.get(), px, py, pz));

	std::vector<Displacement> disp;// = model.Solve();
	std::vector<NodeLoad> react;// = model.Solve();
	model.Solve(loads, disp, react);
	FEStaticResult result = FEStaticResult(std::make_shared<FEModel>(model), loads, disp, react);

	std::cout << "Displacement at 4: " << result.displace[4] << std::endl;
	std::cout << "Displacement at 3: " << result.displace[3] << std::endl;
	std::cout << "Model Node Num: " << model.NodeNum() << std::endl;

}

void TestBeamTorsionMethod1() {

	std::cout << "TestBeamTorsionMethod1 Start" << std::endl;

	Node n0(0, 0, 0); n0.id = 0;
	Node n1(100, 0, 0); n1.id = 1;
	Node n2(100, 100, 0); n2.id = 2;
	n0.Fix.FixAll();

	Material m0(5000, 0.2);
	Section s0(100, 833.33, 833.33, 1406.25);

	std::shared_ptr<BeamElement> b1 = std::make_shared<BeamElement>(0, &n0, &n1, &s0, m0);
	std::shared_ptr<BeamElement> b2 = std::make_shared<BeamElement>(1, &n1, &n2, &s0, m0);

	FEModel model;
	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);

	model.Materials.push_back(m0);

	model.Sections.push_back(s0);

	model.Elements.push_back(b1);
	model.Elements.push_back(b2);

	NodeLoad nl = NodeLoad(2, 0, 0, -10, 0, 0, 0);
	std::vector<std::shared_ptr<LoadBase>> loads;
	loads.push_back(std::make_shared<NodeLoad>(nl));

	std::vector<Displacement> disp;// = model.Solve();
	std::vector<NodeLoad> react;// = model.Solve();
	model.Solve(loads, disp, react);
	FEStaticResult result = FEStaticResult(std::make_shared<FEModel>(model), loads, disp, react);

	std::cout << "Displacement: " << result.displace[2] << std::endl;

	std::cout << "Moment at 0 " << result.GetBeamStress(0, 0) << std::endl;
	std::cout << "Moment at 1 " << result.GetBeamStress(0, 1) << std::endl;

	std::cout << "Model Node Num: " << model.NodeNum() << std::endl;
	// std::cout << "TransMat: \n" << pel.trans_matrix() << std::endl;
}

// �Ў�����(������)�T���v��
void TestBeamSemiRigidMethod() {

	FEModel model;

	Node n0(0, 0, 0); n0.id = 0;
	Node n1(1000, 0, 0); n1.id = 1;
	n1.Fix.FixAll();

	Material m0(5000, 0.2);
	Section s0(600, 45000, 20000, 10000); // 20 x 30

	// ComplexBeamElement b0(&n0, &n1, &s0, m0);
	std::shared_ptr<ComplexBeamElement> b0 = std::make_shared<ComplexBeamElement>(&n0, &n1, &s0, m0);
	b0->Lambda_bz_ = 0.6;
	b0->Lambda_by_ = 0.6;
	//b0->Lambda_sz_ = 0.7;
	//b0->Lambda_sy_ = 0.4;

	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);

	model.Materials.push_back(m0);
	model.Sections.push_back(s0);

	model.Elements.push_back(b0);

	NodeLoad nl = NodeLoad(0, 0, -10, -10, 0, 0, 0);
	std::vector<std::shared_ptr<LoadBase>> loads;
	loads.push_back(std::make_shared<NodeLoad>(nl));
	//model.Loads.push_back(std::make_shared<NodeLoad>(nl));

	std::vector<Displacement> disp;// = model.Solve();
	std::vector<NodeLoad> react;// = model.Solve();
	model.Solve(loads, disp, react);
	FEStaticResult result = FEStaticResult(std::make_shared<FEModel>(model), loads, disp, react);

	std::cout << "React Forces" << std::endl;
	for each (NodeLoad nl in react)
	{
		std::cout << nl << std::endl;
	}
	
	for (size_t i = 0; i < disp.size(); i++)
	{
		std::cout << "id: " << i << std::endl;
		std::cout << disp[i] << std::endl;
	}

	std::cout << "OUTPUT BEAM STRESS" << std::endl;
	std::cout << result.GetBeamStress(0, 0) << std::endl;
	std::cout << result.GetBeamStress(0, 0.5) << std::endl;
	std::cout << result.GetBeamStress(0, 1) << std::endl;

}

void TestBeamDistLoad() {
	std::cout << "TestBeamDistLoad: " << std::endl;
	std::vector<double> w = { -1.0, -2.0, -0.3 };
	std::vector<double> params = { 0.1, 0.6, 0.85 };
	double length = 3000;


	Node n0(0, 0, 0, 0);
	Node n1(1, 3000, 0, 6);
	BeamElement b0(&n0, &n1, nullptr, Material(2.05 * 100000 * 1000, 0.2));
	BeamPolyLoad beam(w, params, &b0, BeamLoadAxis::ZAxis);

	std::cout << "R0: " << beam.R0() << std::endl;
	std::cout << "RA: " << beam.RA() << std::endl;
	std::cout << "M0: " << beam.M0() << std::endl;
	std::cout << "MA: " << beam.MA() << std::endl;

	double EI = 27333.33333;  // kN*m^2 (��)
	std::cout << "Shear Force at x=3: " << beam.shear_force(2800) << std::endl;
	std::cout << "Shear Force at x=3: " << beam.shear_force(2900) << std::endl;
	std::cout << "Bending Moment at x=3: " << beam.bending_moment(2500) << std::endl;
	std::cout << "Deflection at x=3: " << beam.deflection(2000, EI) << std::endl;
}

void TestAxialDistLoad() {
	std::cout << "TestAxialDistLoad: " << std::endl;
	std::vector<double> w = { 1.5, 0.5 };
	std::vector<double> params = { 0.5 / 3, 0.5 };
	double length = 3;


	Node n0(0, 0, 0, 0);
	Node n1(1, 3, 0, 0);
	BeamElement b0(&n0, &n1, nullptr, Material(2.05 * 100000 * 1000, 0.2));
	AxialPolyLoad beam(w, params, &b0);

	std::cout << "N2: " << beam.N0() << std::endl;
	std::cout << "N1: " << beam.N1() << std::endl;

	std::cout << "Axial Force at x=0: " << beam.axial_force(0) << std::endl;
	std::cout << "Axial Force at x=0.85: " << beam.axial_force(0.85) << std::endl;
	std::cout << "Axial Force at x=3.0: " << beam.axial_force(3) << std::endl;


	AxialTrapezoidalLoad beam_trap(1.5, 0.5, 0.5, 1.0, 1.5);
	std::cout << "N1: " << beam_trap.n1() << std::endl;
	std::cout << "N2: " << beam_trap.n2() << std::endl;
	std::cout << "Axial Force at x=0.85: " << beam_trap.axial_force(0.85) << std::endl;
}

void TestDynamicAnalysis() {

	std::cout << "TestDynamicAnalysis Start" << std::endl;

	int divnum = 4;
	FEModel model = CantiBeamModel(200, divnum);

	std::cout << "mode analysis" << std::endl;
	std::vector<double> eigen_values;
	std::vector<std::vector<Displacement>> mode_vectors;
	int computed_modenum = model.SolveVibration(6, eigen_values, mode_vectors);
	
	// eigen_valuesを出力
	std::cout << "Modenum: " << computed_modenum << std::endl;
	for (size_t i = 0; i < eigen_values.size(); i++)
	{
		std::cout << "mode: " << i << ", Eigen Value: " << eigen_values[i] 
			<< ", Natural Period: " << 2 * PI / eigen_values[i] << std::endl;
	}

	double delta_t = 0.01; // 時間刻み
	double total_time = 5.0; // 総時間
	int num_steps = static_cast<int>(total_time / delta_t); // ステップ数
	std::vector<double> gaccels(num_steps);

	std::cout << "num_steps: " << num_steps << std::endl;

	// 振幅
	double amplitude = 1.0;
	// 固有周期
	double natural_period = 0.5;
	
	// sin波を生成する。
	std::vector<double> time(num_steps);
	for (int i = 0; i < num_steps; ++i) {
		double time = i * delta_t;
		gaccels[i] = amplitude * sin(2 * PI * time / natural_period);
	}

	// DynamicAccelLoadを定義
	// ここでは、x軸方向に重力加速度を設定
	DynamicAccelLoad accel_load(delta_t, 0, 0, 1, gaccels);
	
	// gaccelsを出力
	//std::cout << "gaccels: " << std::endl;
	//for (size_t i = 0; i < gaccels.size(); i++)
	//{
	//	std::cout << gaccels[i] << std::endl;
	//}

	std::shared_ptr<FEModel> model_ptr = std::make_shared<FEModel>(model);
	DynamicAnalysis analysis(model_ptr, accel_load);
	analysis.Initialize();

	for (size_t i = 0; i < num_steps; i++)
	{
		analysis.ComputeStep();
		auto disp = analysis.GetDisplacements();
		std::cout << "Step: " << analysis.current_step << ", Time[s]: " << 
			analysis.current_step * analysis.accel_load.timestep << ", Z dir: " << disp[divnum].Dz() << std::endl;
	}
	//analysis.ComputeSteps(10);

	std::cout << "TestDynamicAnalysis End" << std::endl;

}

int main(void) {
	//std::cout << "TestMethod1 Start" << std::endl;
	//TestMethod1();
	//TestMethod2();

	//TestQuadPlateConsistMass();
	//TestTriPlateConsistentMass();
	
	//TestBodyforceToNodeLoadData();
	//TestBeamInertialForceSolve();
	
	//TestPlatePressure();

	//TestBeamTorsionMethod1();

	//CheckCantiBeamVibration();


	// 250511 Debug
	//TestDynamicAnalysis();
	//CheckCantiBeamVibration2();
	// CheckCantiBeamBuckling();

	// 座屈検討用のピラミッド型トラスサンプル
	CheckCantiPyramidTrussBuckling(1000, 4, 100);
	CheckQuadPlateBuckling();

	// CheckSparseSolver();




	//CheckQuadElement1();

	//TestBeamSemiRigidMethod();

}

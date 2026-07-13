#include<iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include<random>
#include<memory>

//#ifdef USE_MKL
//#define EIGEN_USE_MKL_ALL
//#endif

#include "Elements/Elements.h"
#include "Model.h"
#include "FELinearStaticOp.h"
#include "FEBucklingAnalysis.h"
#include "FEVibrateResult.h"
#include "FEAnalysis.h"
#include "FEDynamic.h"
#include "LoadComponent.h"
#include "SparseMatrixUtils.h"
#include "ResponseSpectrumMethod.h"

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

    for (double h : h_compare)
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


// 固有値解析の検証用: 固有値・周期・質量正規化(φ^T M φ)を出力
void PrintVibrationCheck(FEModel& model, int nev, const char* title) {
    std::vector<double> eigs;
    std::vector<std::vector<Displacement>> modes;
    int computed = model.SolveVibration(nev, eigs, modes);
    std::cout << "\n--- VibrationCheck: " << title
              << " (requested=" << nev << ", computed=" << computed << ") ---" << std::endl;
    for (size_t i = 0; i < eigs.size(); i++) {
        double mnorm = 0.0;
        for (size_t j = 0; j < model.Nodes.size(); j++) {
            double mj = model.Nodes[j].MassData.SumMass() / model.GraityAccel;
            const Displacement& d = modes[i][j];
            mnorm += mj * (d.Dx() * d.Dx() + d.Dy() * d.Dy() + d.Dz() * d.Dz());
        }
        std::cout << "  mode " << i + 1
                  << ": omega=" << eigs[i]
                  << "  T=" << 2 * PI / eigs[i]
                  << "  phi^T*M*phi=" << mnorm << std::endl;
    }
}

void TestVibrationCheck() {
    {
        FEModel model = CantiBeamModel(200, 10);
        model.ComputeElementNodeMass();
        PrintVibrationCheck(model, 6, "CantiBeam(200,10)");
    }
    {
        FEModel model = CantiColumnModel(300, 8);
        model.ComputeElementNodeMass();
        PrintVibrationCheck(model, 6, "CantiColumn(300,8)");
    }
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
    for (double var : eigen_values)
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

    FELinearStaticOp result(std::make_shared<FEModel>(model), loads);
    result.Compute();

    std::cout << "Static Result: " << result.displace[model.Nodes.size() - 1] << std::endl;
    std::cout << "Beam Stress at Node 0 At 0: " << result.GetBeamStress(0, 0) << std::endl;

    FEBucklingAnalysis buckling_analysis(std::make_shared<FELinearStaticOp>(result));
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
    }

    return model;
}

void TestVibrationCheckTruss() {
    FEModel model = PyramidTrussModel(1000, 4, 100);
    model.ComputeElementNodeMass();
    PrintVibrationCheck(model, 1, "PyramidTruss(1000,4,100)");
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

    FELinearStaticOp result = FELinearStaticOp(std::make_shared<FEModel>(model), loads);
    result.Compute();

    std::cout << "Static Result: " << result.displace[model.Nodes.size() - 1] << std::endl;
    std::cout << "Beam Stress at Node 0 At 0: " << result.GetBeamStress(0, 0) << std::endl;

    FEBucklingAnalysis buckling_analysis(std::make_shared<FELinearStaticOp>(result));
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
    QuadPlateElement element1(&n0, &n1, &n4, &n3, thickness, mat, 0.1);
    QuadPlateElement element2(&n1, &n2, &n5, &n4, thickness, mat, -0.1);
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

    NodeLoad nl1 = NodeLoad(1, 0, 0, -100.0);
    NodeLoad nl2 = NodeLoad(4, 0, 0, -100.0);
    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.push_back(std::make_shared<NodeLoad>(nl1));
    loads.push_back(std::make_shared<NodeLoad>(nl2));
    std::vector<Displacement> disp;
    std::vector<NodeLoad> react;
    model.SolveLinearStatic(loads, disp, react);
    FELinearStaticOp result = FELinearStaticOp(std::make_shared<FEModel>(model), loads);
    result.Compute();

    std::cout << "Static Result: " << result.displace[model.Nodes.size() - 1] << std::endl;
    FEBucklingAnalysis buckling_analysis(std::make_shared<FELinearStaticOp>(result));
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
    //std::vector<Displacement> displacements;
    //std::vector<NodeLoad> reactions;

    // 静的解析を実行
    //model.SolveLinearStatic(loads, displacements, reactions);
    FELinearStaticOp result = FELinearStaticOp(std::make_shared<FEModel>(model), loads);
    result.Compute();

    // 結果を出力
    std::cout << "Displacements:" << std::endl;
    for (size_t i = 0; i < result.GetDisplacements().size(); ++i) {
        std::cout << "Node " << i << ": " << result.GetDisplacements()[i] << std::endl;
    }

    std::cout << "Reactions:" << std::endl;
    for (const auto& reaction : result.GetReactForces()) {
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
    //model.SolveLinearStatic(loads, disp, react);
    FELinearStaticOp result = FELinearStaticOp(std::make_shared<FEModel>(model), loads);
    result.Compute();

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

    //std::vector<Displacement> disp;// = model.Solve();
    //std::vector<NodeLoad> react;// = model.Solve();
    //model.SolveLinearStatic(loads, disp, react);
    FELinearStaticOp result = FELinearStaticOp(std::make_shared<FEModel>(model), loads);
    result.Compute();

    std::cout << "Displacement: " << result.displace[2] << std::endl;

    std::cout << "Moment at 0 " << result.GetBeamStress(0, 0) << std::endl;
    std::cout << "Moment at 1 " << result.GetBeamStress(0, 1) << std::endl;

    std::cout << "Model Node Num: " << model.NodeNum() << std::endl;
    // std::cout << "TransMat: \n" << pel.trans_matrix() << std::endl;
}

void TestBeamSemiRigidMethod() {

    FEModel model;

    Node n0(0, 0, 0); n0.id = 0;
    Node n1(1000, 0, 0); n1.id = 1;
    n1.Fix.FixAll();

    Material m0(5000, 0.2);
    Section s0(600, 45000, 20000, 10000); // 20 x 30

    // ComplexBeamElement b0(&n0, &n1, &s0, m0);
    std::shared_ptr<ComplexBeamElement> b0 = std::make_shared<ComplexBeamElement>(&n0, &n1, &s0, m0);
    b0->Lambda_bzj = 0.6;
    b0->Lambda_byj = 0.6;
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

    //std::vector<Displacement> disp;// = model.Solve();
    //std::vector<NodeLoad> react;// = model.Solve();
    //model.SolveLinearStatic(loads, disp, react);
    FELinearStaticOp result = FELinearStaticOp(std::make_shared<FEModel>(model), loads);
    result.Compute();

    std::cout << "React Forces" << std::endl;
    for (const NodeLoad& nl : result.GetReactForces())
    {
        std::cout << nl << std::endl;
    }
    
    for (size_t i = 0; i < result.GetDisplacements().size(); i++)
    {
        std::cout << "id: " << i << std::endl;
        std::cout << result.GetDisplacements()[i] << std::endl;
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

// 1次共振加振→自由振動の時刻歴を計算し、先端変位履歴と最大変位を返す
double RunDampedResonance(std::shared_ptr<FEModel> model_ptr, FEDynamicDampInitializer* damp,
    double T1, int tip_node, int steps_per_cycle, int excite_cycles, int free_cycles,
    std::vector<double>& tip_history)
{
    double dt = T1 / steps_per_cycle;
    int excite_steps = steps_per_cycle * excite_cycles;
    int num_steps = excite_steps + steps_per_cycle * free_cycles;
    std::vector<double> gaccels(num_steps, 0.0);
    for (int i = 0; i < excite_steps; ++i)
        gaccels[i] = sin(2 * PI * (i * dt) / T1);

    DynamicAccelLoad accel_load(dt, 0, 0, 1, gaccels);
    DynamicAnalysis analysis(model_ptr, accel_load, damp);
    analysis.RecordEnabled = false;
    if (!analysis.Initialize()) {
        std::cout << "Failed to initialize dynamic analysis." << std::endl;
        return -1.0;
    }

    double peak = 0.0;
    tip_history.clear();
    for (int i = 0; i < num_steps; ++i) {
        analysis.ComputeStep();
        double dz = analysis.GetDisplacements()[tip_node].Dz();
        tip_history.push_back(dz);
        peak = std::max(peak, std::abs(dz));
    }
    return peak;
}

// 自由振動部分の正ピークの対数減衰率から減衰比を推定
double EstimateDampingFromDecay(const std::vector<double>& hist, size_t start)
{
    std::vector<double> peaks;
    for (size_t i = start + 1; i + 1 < hist.size(); i++) {
        if (hist[i] > hist[i - 1] && hist[i] > hist[i + 1] && hist[i] > 0)
            peaks.push_back(hist[i]);
    }
    if (peaks.size() < 3)
        return -1.0;
    double delta = std::log(peaks.front() / peaks.back()) / (double)(peaks.size() - 1);
    return delta / (2 * PI);
}

// 減衰初期化子(剛性比例・質量比例・レイリー)の比較検証
void TestDampInitializers() {

    std::cout << "TestDampInitializers Start" << std::endl;

    int divnum = 4;
    FEModel model = CantiBeamModel(200, divnum);
    std::shared_ptr<FEModel> model_ptr = std::make_shared<FEModel>(model);
    model_ptr->ComputeElementNodeMass();

    // 固有値解析(断面が対称でモードが縮退するため、1次と異なる振動数のモードを探す)
    std::vector<double> eigen_values;
    std::vector<std::vector<Displacement>> mode_vectors;
    int nconv = model_ptr->SolveVibration(4, eigen_values, mode_vectors);
    if (nconv < 2) {
        std::cout << "SolveVibration failed." << std::endl;
        return;
    }
    double w1 = eigen_values[0];
    double T1 = 2 * PI / w1;
    int mode_j = 2; // 1-based
    while (mode_j <= nconv && std::abs(eigen_values[mode_j - 1] - w1) / w1 < 1e-3)
        mode_j++;
    if (mode_j > nconv) {
        std::cout << "No distinct second mode found." << std::endl;
        return;
    }
    double wj = eigen_values[mode_j - 1];
    std::cout << "w1: " << w1 << " (T1: " << T1 << " s), mode_j: " << mode_j
        << ", wj: " << wj << std::endl;

    const double zeta = 0.05;
    const int steps_per_cycle = 40, excite_cycles = 20, free_cycles = 20;
    const size_t free_start = (size_t)steps_per_cycle * excite_cycles;
    std::vector<double> hist;

    // Case 1: 剛性比例(既存)
    FEDynamicStiffDampInitializer stiff_damp(zeta);
    double peak_stiff = RunDampedResonance(model_ptr, &stiff_damp, T1, divnum,
        steps_per_cycle, excite_cycles, free_cycles, hist);
    double zeta_stiff = EstimateDampingFromDecay(hist, free_start);
    std::cout << "[Stiffness] peak: " << peak_stiff << ", estimated zeta: " << zeta_stiff
        << ", w1(internal): " << stiff_damp.natural_angle_velocity << std::endl;

    // Case 2: 質量比例
    FEDynamicMassDampInitializer mass_damp(zeta);
    double peak_mass = RunDampedResonance(model_ptr, &mass_damp, T1, divnum,
        steps_per_cycle, excite_cycles, free_cycles, hist);
    double zeta_mass = EstimateDampingFromDecay(hist, free_start);
    std::cout << "[Mass]      peak: " << peak_mass << ", estimated zeta: " << zeta_mass
        << ", w1(internal): " << mass_damp.natural_angle_velocity << std::endl;

    // Case 3: レイリー(1次・mode_j次モードで zeta を指定)
    FEDynamicRayleighDampInitializer rayleigh_damp(zeta, zeta, 1, mode_j);
    double peak_ray = RunDampedResonance(model_ptr, &rayleigh_damp, T1, divnum,
        steps_per_cycle, excite_cycles, free_cycles, hist);
    double zeta_ray = EstimateDampingFromDecay(hist, free_start);
    std::cout << "[Rayleigh]  peak: " << peak_ray << ", estimated zeta: " << zeta_ray
        << ", alpha: " << rayleigh_damp.alpha << ", beta: " << rayleigh_damp.beta << std::endl;

    // 算出されたalpha, betaによる各モードの減衰比を逆算(= zetaになるはず)
    std::cout << "  zeta(w1): " << rayleigh_damp.alpha / (2 * w1) + rayleigh_damp.beta * w1 / 2
        << ", zeta(wj): " << rayleigh_damp.alpha / (2 * wj) + rayleigh_damp.beta * wj / 2 << std::endl;

    // Case 4: レイリー(alpha, beta 直接指定; Case 3 と同値になるはず)
    double alpha_direct = 2 * zeta * w1 * wj / (w1 + wj);
    double beta_direct = 2 * zeta / (w1 + wj);
    FEDynamicRayleighDampInitializer rayleigh_direct(alpha_direct, beta_direct);
    double peak_ray2 = RunDampedResonance(model_ptr, &rayleigh_direct, T1, divnum,
        steps_per_cycle, excite_cycles, free_cycles, hist);
    std::cout << "[RayDirect] peak: " << peak_ray2
        << ", alpha: " << alpha_direct << ", beta: " << beta_direct
        << ", diff vs Case3: " << std::abs(peak_ray2 - peak_ray) << std::endl;

    // DampRateAtPeriod の検証
    double Tj = 2 * PI / wj;
    std::cout << "DampRateAtPeriod checks:" << std::endl;
    std::cout << "  [Stiffness] at T1: " << stiff_damp.DampRateAtPeriod(T1)
        << " (expected " << zeta << "), at Tj: " << stiff_damp.DampRateAtPeriod(Tj)
        << " (expected " << zeta * wj / w1 << ")" << std::endl;
    std::cout << "  [Mass]      at T1: " << mass_damp.DampRateAtPeriod(T1)
        << " (expected " << zeta << "), at Tj: " << mass_damp.DampRateAtPeriod(Tj)
        << " (expected " << zeta * w1 / wj << ")" << std::endl;
    std::cout << "  [Rayleigh]  at T1: " << rayleigh_damp.DampRateAtPeriod(T1)
        << ", at Tj: " << rayleigh_damp.DampRateAtPeriod(Tj)
        << " (both expected " << zeta << ")" << std::endl;

    // alpha, beta 直接指定はInitialize前でも算定可能
    FEDynamicRayleighDampInitializer ray_pre(alpha_direct, beta_direct);
    std::cout << "  [RayDirect pre-Init] at T1: " << ray_pre.DampRateAtPeriod(T1)
        << " (expected " << zeta << ")" << std::endl;

    // Initialize前(w1未確定)は -1
    FEDynamicMassDampInitializer mass_pre(zeta);
    std::cout << "  [Mass pre-Init] returns: " << mass_pre.DampRateAtPeriod(T1)
        << " (expected -1)" << std::endl;

    std::cout << "TestDampInitializers End" << std::endl;
}

void TestSimplaFrame() {

    std::cout << "TestSimplaFrame Start" << std::endl;

    FEModel model;

    Node n0(0, 0, 0); n0.id = 0;
    Node n1(1000, 0, 0); n1.id = 1;
    Node n2(0, 0, 1000); n2.id = 2;
    Node n3(1000, 0, 1000); n3.id = 3;
    n0.Fix.FixAll();
    n1.Fix.FixAll();

    Material m0(5000, 0.2);
    Section s0(600, 45000, 20000, 10000); // 20 x 30

    // ComplexBeamElement b0(&n0, &n1, &s0, m0);
    std::shared_ptr<ComplexBeamElement> b0 = std::make_shared<ComplexBeamElement>(&n2, &n3, &s0, m0);
    b0->Lambda_bzj = 0.01;
    //b0->Lambda_byj = 0.6;
    //b0->Lambda_sz_ = 0.7;
    //b0->Lambda_sy_ = 0.4;
    std::shared_ptr<BeamElement> b1 = std::make_shared<BeamElement>(&n0, &n2, &s0, m0);
    std::shared_ptr<BeamElement> b2 = std::make_shared<BeamElement>(&n1, &n3, &s0, m0);


    model.Nodes.push_back(n0);
    model.Nodes.push_back(n1);
    model.Nodes.push_back(n2);
    model.Nodes.push_back(n3);

    model.Materials.push_back(m0);
    model.Sections.push_back(s0);

    model.Elements.push_back(b0);
    model.Elements.push_back(b1);
    model.Elements.push_back(b2);

    NodeLoad nl = NodeLoad(2, 10, 0, 0, 0, 0, 0);
    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.push_back(std::make_shared<NodeLoad>(nl));
    //model.Loads.push_back(std::make_shared<NodeLoad>(nl));

    //std::vector<Displacement> disp;// = model.Solve();
    //std::vector<NodeLoad> react;// = model.Solve();
    //model.SolveLinearStatic(loads, disp, react);
    FELinearStaticOp result = FELinearStaticOp(std::make_shared<FEModel>(model), loads);
    result.Compute();

    std::cout << "React Forces" << std::endl;
    for (const NodeLoad& nl : result.GetReactForces())
    {
        std::cout << nl << std::endl;
    }

    for (size_t i = 0; i < result.GetDisplacements().size(); i++)
    {
        std::cout << "id: " << i << std::endl;
        std::cout << result.GetDisplacements()[i] << std::endl;
    }

    std::cout << "OUTPUT BEAM STRESS" << std::endl;
    std::cout << result.GetBeamStress(0, 0) << std::endl;
    std::cout << result.GetBeamStress(0, 0.5) << std::endl;
    std::cout << result.GetBeamStress(0, 1) << std::endl;

    std::cout << "OUTPUT BEAM Displace" << std::endl;
    std::cout << result.GetBeamDisplace(0, 0) << std::endl;
    std::cout << result.GetBeamDisplace(0, 0.5) << std::endl;
    std::cout << result.GetBeamDisplace(0, 1) << std::endl;

    std::cout << "TestSimplaFrame End" << std::endl;

}

// ヘルパー関数: 疎行列をビジュアル表示
void printSparseMatrix(const char* name, const Eigen::SparseMatrix<double>& mat) {
    std::cout << "\n" << name << " (" << mat.rows() << "x" << mat.cols() << "):" << std::endl;

    // 密行列に変換して表示
    Eigen::MatrixXd dense = Eigen::MatrixXd(mat);

    for (int i = 0; i < dense.rows(); ++i) {
        std::cout << "  [ ";
        for (int j = 0; j < dense.cols(); ++j) {
            if (std::abs(dense(i, j)) < 1e-10) {
                std::cout << "  .  ";  // ゼロ要素は . で表示
            } else {
                std::printf("%5.1f", dense(i, j));
                std::cout << " ";
            }
        }
        std::cout << "]" << std::endl;
    }
}

/// <summary>
/// 節点IDから節点ポインタを検索するヘルパー関数
/// </summary>
Node* findNodeById(std::vector<Node>& nodes, int id) {
    for (auto& node : nodes) {
        if (node.id == id) {
            return &node;
        }
    }
    return nullptr;
}

/// <summary>
/// 階ごとの並進変位をサマリー表示
/// </summary>
void PrintFloorDisplacements(
    const Eigen::VectorXd& displacement,
    const RigidLinks& rigidLinks,
    int numStories,
    const std::string& caseLabel)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "Floor Displacements Summary: " << caseLabel << std::endl;
    std::cout << "========================================" << std::endl;

    for (int floor = 1; floor <= numStories; floor++) {
        std::cout << "\n--- Floor " << floor << " ---" << std::endl;

        if (floor - 1 < rigidLinks.links.size()) {
            const RigidLink& link = rigidLinks.links[floor - 1];

            // マスター節点情報（位置のみ）
            int masterNodeId = link.Master.id;
            std::cout << "Master Node " << masterNodeId << " @ ("
                      << link.Master.Location.x << ", "
                      << link.Master.Location.y << ")" << std::endl;

            // 各スレーブ節点の変位から逆算したマスター変位
            std::cout << "\nSlave nodes (calculated master from each):" << std::endl;

            std::vector<double> ux_m_calcs, uy_m_calcs;
            for (size_t i = 0; i < link.Slaves.size(); i++) {
                const Node& slave = link.Slaves[i];
                int slaveNodeId = slave.id;

                double ux_s = displacement(slaveNodeId * 6 + 0);
                double uy_s = displacement(slaveNodeId * 6 + 1);
                double rz_s = displacement(slaveNodeId * 6 + 5);

                double dx = slave.Location.x - link.Master.Location.x;
                double dy = slave.Location.y - link.Master.Location.y;

                double ux_m_calc = ux_s + rz_s * dy;
                double uy_m_calc = uy_s - rz_s * dx;

                ux_m_calcs.push_back(ux_m_calc);
                uy_m_calcs.push_back(uy_m_calc);

                std::cout << "  Node " << slaveNodeId << ": UX=" << ux_s << ", UY=" << uy_s << ", RZ=" << rz_s
                          << " -> Master: UX=" << ux_m_calc << ", UY=" << uy_m_calc << std::endl;
            }

            // スレーブ間での計算されたマスター変位の整合性チェック
            if (ux_m_calcs.size() > 1) {
                double max_err_ux = 0.0, max_err_uy = 0.0;
                for (size_t i = 1; i < ux_m_calcs.size(); i++) {
                    max_err_ux = std::max(max_err_ux, std::abs(ux_m_calcs[i] - ux_m_calcs[0]));
                    max_err_uy = std::max(max_err_uy, std::abs(uy_m_calcs[i] - uy_m_calcs[0]));
                }

                std::cout << "\nSlave Consistency: Max Error UX=" << max_err_ux
                          << " mm, UY=" << max_err_uy << " mm";
                if (max_err_ux < 1e-4 && max_err_uy < 1e-4) {
                    std::cout << "  OK" << std::endl;
                } else {
                    std::cout << "  FAIL" << std::endl;
                }
            }
        }
    }

    std::cout << "\n========================================" << std::endl;
}

/// <summary>
/// 剛体連結を用いた剛床サンプル
/// 節点順序をシャッフルしてインデクス順依存を検証する。
/// </summary>
void TestRigidFloorWithCenterMaster_Shuffled() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test: Rigid Floor with CENTER Master Nodes (SHUFFLED)" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // グリッド寸法
    double span = 5000.0;
    double height = 3000.0;
    int numStories = 3;

    Material steel(205000.0, 0.3, 79000.0);
    Section columnSection(3468.32, 14693530.62, 14693530.62, 6830000.0);
    Section beamSection(5558.94, 153247657.08, 10674961.13, 3520000.0);

    std::cout << "Building model with CENTER master nodes (SHUFFLED ORDER)..." << std::endl;
    std::cout << "  Master nodes at: (" << span/2 << ", " << span/2 << ", z)" << std::endl;
    std::cout << "  All frame nodes are slaves\n" << std::endl;

    FEModel model;
    model.Materials.push_back(steel);
    model.Sections.push_back(columnSection);
    model.Sections.push_back(beamSection);

    std::vector<std::array<double, 2>> gridPositions = {
        {0, 0}, {span, 0}, {span, span}, {0, span}
    };

    // ============================================
    // 節点情報を構造体に格納
    // ============================================
    struct NodeInfo {
        int id;
        Point position;
        bool isFixed;
    };

    std::vector<NodeInfo> nodeInfos;

    // 乱数生成器の初期化（固定シードで再現性を確保）
    std::random_device rd;
    std::mt19937 gen(12345);

    // 基礎レベル (z=0) - 完全固定
    for (int i = 0; i < 4; i++) {
        NodeInfo info;
        info.id = i;
        info.position = Point(gridPositions[i][0], gridPositions[i][1], 0);
        info.isFixed = true;
        nodeInfos.push_back(info);
    }

    // 節点IDをシャッフルするための準備（基礎を除く全節点ID）
    std::vector<int> shuffledNodeIds;
    for (int floor = 1; floor <= numStories; floor++) {
        for (int i = 0; i < 4; i++) {
            shuffledNodeIds.push_back(floor * 4 + i);
        }
    }

    // 節点IDをシャッフル（高さ方向にもばらつき）
    std::shuffle(shuffledNodeIds.begin(), shuffledNodeIds.end(), gen);

    // シャッフルされたIDを各位置に割り当て
    int idIndex = 0;
    for (int floor = 1; floor <= numStories; floor++) {
        double z = height * floor;
        for (int i = 0; i < 4; i++) {
            NodeInfo info;
            info.id = shuffledNodeIds[idIndex++];  // シャッフルされたID
            info.position = Point(gridPositions[i][0], gridPositions[i][1], z);
            info.isFixed = false;
            nodeInfos.push_back(info);
        }
    }

    // ============================================
    // 節点をID順にソートしてmodel.Nodesに追加
    // （Nodesには必ずID順で登録するルール）
    // ============================================
    std::sort(nodeInfos.begin(), nodeInfos.end(),
        [](const NodeInfo& a, const NodeInfo& b) { return a.id < b.id; });

    for (const auto& info : nodeInfos) {
        Node node(info.position);
        node.id = info.id;
        if (info.isFixed) {
            node.Fix.FixAll();
        }
        model.Nodes.push_back(node);
    }

    std::cout << "Created " << model.Nodes.size() << " frame nodes (ID順)" << std::endl;

    // 節点IDの分布を表示
    std::cout << "\nNode ID distribution by floor (shuffled):" << std::endl;
    for (int floor = 1; floor <= numStories; floor++) {
        double z = height * floor;
        std::cout << "  Floor " << floor << " (z=" << z << "): IDs = ";
        for (const auto& info : nodeInfos) {
            if (!info.isFixed && std::abs(info.position.z - z) < 1e-6) {
                std::cout << info.id << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // 仮想マスター節点を準備（モデルのNodesには追加しない！）
    std::vector<Node> virtualMasterNodes;
    int masterNodeStartId = 1000;

    for (int floor = 1; floor <= numStories; floor++) {
        double z = height * floor;
        Point p(span / 2.0, span / 2.0, z);
        Node node(p);
        node.id = masterNodeStartId + (floor - 1);
        virtualMasterNodes.push_back(node);
    }

    std::cout << "Created " << numStories << " virtual master nodes (NOT in model.Nodes)" << std::endl;
    std::cout << "  Virtual master positions: (" << span/2 << ", " << span/2 << ", z)" << std::endl;

    // ============================================
    // 要素作成（柱 + 梁）- IDベースで節点を検索
    // ============================================
    size_t elemId = 0;

    // 柱要素
    for (int floor = 0; floor < numStories; floor++) {
        for (int col = 0; col < 4; col++) {
            int nodeIId = floor * 4 + col;
            int nodeJId = (floor + 1) * 4 + col;

            Node* nodeI = findNodeById(model.Nodes, nodeIId);
            Node* nodeJ = findNodeById(model.Nodes, nodeJId);

            if (nodeI && nodeJ) {
                BeamElement elem(elemId++, nodeI, nodeJ,
                    &model.Sections[0], model.Materials[0]);
                model.add_element(elem);
            }
        }
    }

    size_t numColumns = elemId;

    // 梁要素
    for (int floor = 1; floor <= numStories; floor++) {
        int baseNodeId = floor * 4;
        std::vector<std::array<int, 2>> beamConnections = {
            {0, 1}, {1, 2}, {2, 3}, {3, 0}
        };

        for (auto& conn : beamConnections) {
            int nodeIId = baseNodeId + conn[0];
            int nodeJId = baseNodeId + conn[1];

            Node* nodeI = findNodeById(model.Nodes, nodeIId);
            Node* nodeJ = findNodeById(model.Nodes, nodeJId);

            if (nodeI && nodeJ) {
                BeamElement elem(elemId++, nodeI, nodeJ,
                    &model.Sections[1], model.Materials[0]);
                model.add_element(elem);
            }
        }
    }

    std::cout << "Created " << numColumns << " columns, " << (elemId - numColumns) << " beams" << std::endl;

    // 荷重条件: 幾何学的位置ベースで荷重点を決定
    int loadFloor = 2;
    int loadNodeLocal = 1;  // gridPositions[1] = (5000, 0)
    double loadZ = height * loadFloor;
    Point loadPos = Point(gridPositions[loadNodeLocal][0], gridPositions[loadNodeLocal][1], loadZ);
    double loadValue = 100000.0; // 100kN

    // この幾何学的位置にある節点IDを探す
    int loadNodeId = -1;
    for (const auto& info : nodeInfos) {
        if (!info.isFixed &&
            std::abs(info.position.x - loadPos.x) < 1e-6 &&
            std::abs(info.position.y - loadPos.y) < 1e-6 &&
            std::abs(info.position.z - loadPos.z) < 1e-6) {
            loadNodeId = info.id;
            break;
        }
    }

    if (loadNodeId == -1) {
        std::cerr << "Error: Could not find load node at position " << loadPos << std::endl;
        return;
    }

    NodeLoad nl(loadNodeId, loadValue, 0, 0);
    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.push_back(std::make_shared<NodeLoad>(nl));

    std::cout << "Applied load: " << loadValue << " N at Node " << loadNodeId
              << " (geometric position: " << loadPos << ")\n" << std::endl;

    // ============================================
    // 剛床設定：床中心にマスター節点
    // ============================================

    for (int floor = 1; floor <= numStories; floor++) {
        RigidLink link;

        link.flags[0] = true;   // UX
        link.flags[1] = true;   // UY
        link.flags[2] = false;  // UZ
        link.flags[3] = false;  // RX
        link.flags[4] = false;  // RY
        link.flags[5] = true;   // RZ

        // マスター節点: 床中心の仮想節点
        link.Master = virtualMasterNodes[floor - 1];

        // スレーブ節点: その階の幾何学的位置にある節点を収集
        double floorZ = height * floor;
        std::vector<Node> floorSlaves;

        // nodeInfosから該当する階の節点IDを取得
        for (const auto& info : nodeInfos) {
            if (!info.isFixed && std::abs(info.position.z - floorZ) < 1e-6) {
                Node* node = findNodeById(model.Nodes, info.id);
                if (node) {
                    floorSlaves.push_back(*node);
                }
            }
        }

        // スレーブリストをさらにシャッフル（階内での順序もランダム化）
        std::shuffle(floorSlaves.begin(), floorSlaves.end(), gen);
        for (size_t i = 0; i < floorSlaves.size(); i++) {
            link.Slaves.push_back(floorSlaves[i]);
        }

        model.RigidLinkData->links.push_back(link);
    }

    // ============================================
    // 解析実行
    // ============================================
    std::cout << "========================================" << std::endl;
    std::cout << "Solving with CENTER master nodes (SHUFFLED)..." << std::endl;
    std::cout << "========================================\n" << std::endl;

    std::vector<Displacement> disp_center;
    std::vector<NodeLoad> react_center;
    model.SolveLinearStatic(loads, disp_center, react_center);

    // 結果表示
    std::cout << "\n  Note: Virtual master nodes are not in Displacement vector" << std::endl;
    std::cout << "        Master DOF values are used only for transformation" << std::endl;

    // 反力チェック
    double sumRX = 0.0, sumRY = 0.0, sumRZ = 0.0;
    for (size_t i = 0; i < model.Nodes.size(); i++) {
        if (model.Nodes[i].Fix.IsAnyFix()) {
            sumRX += react_center[i].Px();
            sumRY += react_center[i].Py();
            sumRZ += react_center[i].Pz();
        }
    }

    std::cout << "\n--- Reaction Force Sum ---" << std::endl;
    std::cout << "  ΣRX = " << sumRX << " N (applied: " << loadValue << " N)" << std::endl;
    std::cout << "  ΣRY = " << sumRY << " N (should be ~0)" << std::endl;
    std::cout << "  ΣRZ = " << sumRZ << " N (should be ~0)" << std::endl;
    std::cout << "  Equilibrium: |ΣRX - Load| = " << std::abs(sumRX - loadValue) << " N" << std::endl;

    // ============================================
    // 剛体変形の検証
    // ============================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "Rigid Body Transformation Verification" << std::endl;
    std::cout << "========================================\n" << std::endl;

    Eigen::VectorXd dispVec = Eigen::VectorXd::Zero(disp_center.size() * 6);
    for (size_t i = 0; i < disp_center.size(); i++) {
        dispVec(i * 6 + 0) = disp_center[i].Dx();
        dispVec(i * 6 + 1) = disp_center[i].Dy();
        dispVec(i * 6 + 2) = disp_center[i].Dz();
        dispVec(i * 6 + 3) = disp_center[i].Rx();
        dispVec(i * 6 + 4) = disp_center[i].Ry();
        dispVec(i * 6 + 5) = disp_center[i].Rz();
    }

    // 階ごとの変位サマリー表示
    PrintFloorDisplacements(dispVec, *model.RigidLinkData, numStories, "Center Master (Shuffled)");

    // 剛体リンク付きモデルの固有値解析チェック
    model.ComputeElementNodeMass();
    PrintVibrationCheck(model, 8, "RigidFloor(CenterMaster,Shuffled)");

    std::cout << "\n========================================" << std::endl;
    std::cout << "Center Master Test (Shuffled) Complete" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

void TestSparseMatrixUtils() {

    std::cout << "\n========================================" << std::endl;
    std::cout << "TestSparseMatrixUtils Start" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // テスト用の6x6疎行列を作成
    // 対称行列の例:
    // [4  1  0  0  1  0]
    // [1  5  2  0  0  0]
    // [0  2  6  3  0  0]
    // [0  0  3  7  0  1]
    // [1  0  0  0  8  2]
    // [0  0  0  1  2  9]

    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;

    // 対称行列なので上三角のみ定義
    triplets.emplace_back(0, 0, 4.0);
    triplets.emplace_back(0, 1, 1.0);
    triplets.emplace_back(0, 4, 1.0);
    triplets.emplace_back(1, 1, 5.0);
    triplets.emplace_back(1, 2, 2.0);
    triplets.emplace_back(2, 2, 6.0);
    triplets.emplace_back(2, 3, 3.0);
    triplets.emplace_back(3, 3, 7.0);
    triplets.emplace_back(3, 5, 1.0);
    triplets.emplace_back(4, 4, 8.0);
    triplets.emplace_back(4, 5, 2.0);
    triplets.emplace_back(5, 5, 9.0);

    Eigen::SparseMatrix<double> A(6, 6);
    A.setFromTriplets(triplets.begin(), triplets.end());

    printSparseMatrix("Original matrix A", A);

    // =========================================
    // Test 1: splitMatrixWithResize (2x2 split, extract free_matrix only)
    // =========================================
    std::cout << "\n----------------------------------------" << std::endl;
    std::cout << "Test 1: splitMatrixWithResize (free_matrix only)" << std::endl;
    std::cout << "  Fixed indices: [1, 3]" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    std::vector<int> fixed_indices1 = { 1, 3 };
    Eigen::SparseMatrix<double> free_matrix;

    SparseMatrixUtils::splitMatrixWithResize(A, fixed_indices1, free_matrix);

    printSparseMatrix("Free DOF part (free_matrix)", free_matrix);

    // =========================================
    // Test 2: splitMatrixWithResize (2x2 split, extract all blocks)
    // =========================================
    std::cout << "\n----------------------------------------" << std::endl;
    std::cout << "Test 2: splitMatrixWithResize (all blocks)" << std::endl;
    std::cout << "  Fixed indices: [1, 3]" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    Eigen::SparseMatrix<double> free_matrix2, free_fixed_matrix, fixed_matrix;

    SparseMatrixUtils::splitMatrixWithResize(A, fixed_indices1,
        free_matrix2, free_fixed_matrix, fixed_matrix);

    printSparseMatrix("Free-Free block (free_matrix)", free_matrix2);
    printSparseMatrix("Free-Fixed block (free_fixed_matrix)", free_fixed_matrix);
    printSparseMatrix("Fixed-Fixed block (fixed_matrix)", fixed_matrix);

    // =========================================
    // Test 3: mergeMatrixWithResize (2x2 merge)
    // =========================================
    std::cout << "\n----------------------------------------" << std::endl;
    std::cout << "Test 3: mergeMatrixWithResize (2x2 merge)" << std::endl;
    std::cout << "  Reassemble split matrices" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    Eigen::SparseMatrix<double> A_merged;
    SparseMatrixUtils::mergeMatrixWithResize(free_matrix2, free_fixed_matrix,
        fixed_matrix, A_merged);

    printSparseMatrix("Merged matrix (A_merged)", A_merged);

    // =========================================
    // テスト4: splitMatrix3x3 (3x3分割)
    // =========================================
    std::cout << "\n----------------------------------------" << std::endl;
    std::cout << "Test 4: splitMatrix3x3 (3x3 split)" << std::endl;
    std::cout << "  Group1 indices: [1, 2]" << std::endl;
    std::cout << "  Group2 indices: [4]" << std::endl;
    std::cout << "  Remaining indices are free DOF group" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    std::vector<int> indices_group1 = { 1, 2 };
    std::vector<int> indices_group2 = { 4 };
    Eigen::SparseMatrix<double> mat_11, mat_12, mat_13, mat_22, mat_23, mat_33;

    SparseMatrixUtils::splitMatrix3x3(A, indices_group1, indices_group2,
        mat_11, mat_12, mat_13, mat_22, mat_23, mat_33);

    printSparseMatrix("Block(1,1) - Group1 x Group1", mat_11);
    printSparseMatrix("Block(1,2) - Group1 x Group2", mat_12);
    printSparseMatrix("Block(1,3) - Group1 x Free", mat_13);
    printSparseMatrix("Block(2,2) - Group2 x Group2", mat_22);
    printSparseMatrix("Block(2,3) - Group2 x Free", mat_23);
    printSparseMatrix("Block(3,3) - Free x Free", mat_33);

    // =========================================
    // Test 5: vstack (vertical stacking)
    // =========================================
    std::cout << "\n----------------------------------------" << std::endl;
    std::cout << "Test 5: vstack (vertical stacking)" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // Create two 2x3 matrices
    std::vector<T> triplets_A;
    triplets_A.emplace_back(0, 0, 1.0);
    triplets_A.emplace_back(0, 2, 2.0);
    triplets_A.emplace_back(1, 1, 3.0);

    std::vector<T> triplets_B;
    triplets_B.emplace_back(0, 0, 4.0);
    triplets_B.emplace_back(1, 1, 5.0);
    triplets_B.emplace_back(1, 2, 6.0);

    Eigen::SparseMatrix<double> A_vstack(2, 3);
    Eigen::SparseMatrix<double> B_vstack(2, 3);
    A_vstack.setFromTriplets(triplets_A.begin(), triplets_A.end());
    B_vstack.setFromTriplets(triplets_B.begin(), triplets_B.end());

    printSparseMatrix("Matrix A (2x3)", A_vstack);
    printSparseMatrix("Matrix B (2x3)", B_vstack);

    Eigen::SparseMatrix<double> C_vstack = SparseMatrixUtils::vstack(A_vstack, B_vstack);
    printSparseMatrix("Vertical stack result C = vstack(A, B) (4x3)", C_vstack);

    // =========================================
    // Test 6: hstack (horizontal stacking)
    // =========================================
    std::cout << "\n----------------------------------------" << std::endl;
    std::cout << "Test 6: hstack (horizontal stacking)" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // Create two 3x2 matrices
    std::vector<T> triplets_A2;
    triplets_A2.emplace_back(0, 0, 1.0);
    triplets_A2.emplace_back(1, 1, 2.0);
    triplets_A2.emplace_back(2, 0, 3.0);

    std::vector<T> triplets_B2;
    triplets_B2.emplace_back(0, 1, 4.0);
    triplets_B2.emplace_back(1, 0, 5.0);
    triplets_B2.emplace_back(2, 1, 6.0);

    Eigen::SparseMatrix<double> A_hstack(3, 2);
    Eigen::SparseMatrix<double> B_hstack(3, 2);
    A_hstack.setFromTriplets(triplets_A2.begin(), triplets_A2.end());
    B_hstack.setFromTriplets(triplets_B2.begin(), triplets_B2.end());

    printSparseMatrix("Matrix A (3x2)", A_hstack);
    printSparseMatrix("Matrix B (3x2)", B_hstack);

    Eigen::SparseMatrix<double> C_hstack = SparseMatrixUtils::hstack(A_hstack, B_hstack);
    printSparseMatrix("Horizontal stack result C = hstack(A, B) (3x4)", C_hstack);

    std::cout << "\n========================================" << std::endl;
    std::cout << "TestSparseMatrixUtils End" << std::endl;
    std::cout << "========================================\n" << std::endl;
}


// =============================================================
// シェル構造（浅い球状ドーム）モデル + CQC法ベンチマーク
// =============================================================

// 浅い球状ドーム: 21x21 = 441節点、周辺ピン
//  span: 平面投影スパン[mm]、rise: 中心ライズ[mm]、ndiv: 1辺の分割数
FEModel ShallowDomeModel(double span, double rise, int ndiv) {
    FEModel model;
    int n = ndiv + 1; // 1辺の節点数
    double L = span;
    double H = rise;
    double R = L * 0.5;

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            double x = -L * 0.5 + L * i / (double)(n - 1);
            double y = -L * 0.5 + L * j / (double)(n - 1);
            double r2 = x * x + y * y;
            double z = H * std::max(0.0, 1.0 - r2 / (R * R));
            int id = j * n + i;
            Node node(id, x, y, z);
            // 周辺ピン
            if (i == 0 || j == 0 || i == n - 1 || j == n - 1) {
                node.Fix.PinFix();
            }
            model.Nodes.push_back(node);
        }
    }

    // 鋼材（N, mm, t）相当
    Material mat(2.05e5, 0.3, 7.85e-9);
    model.Materials.push_back(mat);

    Thickness thickness(50.0); // 50mm厚

    // 各クアッドを2つの三角形に分割
    int eid = 0;
    for (int j = 0; j < n - 1; ++j) {
        for (int i = 0; i < n - 1; ++i) {
            int a = j * n + i;
            int b = j * n + (i + 1);
            int c = (j + 1) * n + (i + 1);
            int d = (j + 1) * n + i;
            auto el1 = std::make_shared<TriPlateElement>(
                eid++, &model.Nodes[a], &model.Nodes[b], &model.Nodes[c],
                thickness, model.Materials[0]);
            auto el2 = std::make_shared<TriPlateElement>(
                eid++, &model.Nodes[a], &model.Nodes[c], &model.Nodes[d],
                thickness, model.Materials[0]);
            model.Elements.push_back(el1);
            model.Elements.push_back(el2);
        }
    }

    return model;
}

// 告示風の応答スペクトル(概略形)
class NotifiedSpectrum : public IResponseSpectrum {
public:
    double Acceleration(double T) override {
        // mm/s^2
        double Sa_ms2;
        if (T < 0.16) {
            Sa_ms2 = 3.2 + 30.0 * T;
        } else if (T < 0.864) {
            Sa_ms2 = 8.0;
        } else {
            Sa_ms2 = 8.0 * 0.864 / T;
        }
        return Sa_ms2 * 1000.0;
    }
    double Velocity(double T) override {
        return Acceleration(T) * T / (2.0 * PI);
    }
    double Displacement(double T) override {
        double r = T / (2.0 * PI);
        return Acceleration(T) * r * r;
    }
};

void BenchResponseSpectrumCQC() {
    using clock = std::chrono::high_resolution_clock;
    auto sec = [](auto a, auto b) {
        return std::chrono::duration<double>(b - a).count();
    };

    std::cout << "\n=== BenchResponseSpectrumCQC ===" << std::endl;

    // (1) モデル構築
    auto t0 = clock::now();
    FEModel model = ShallowDomeModel(20000.0, 2000.0, 20); // 21x21=441節点
    model.ComputeElementNodeMass();
    auto t1 = clock::now();
    int N = model.NodeNum();
    int freeDof = model.FreeDOFNum();
    std::cout << "[1] Build model       : " << sec(t0, t1) << " s   "
              << "Nodes=" << N << " Elems=" << model.Elements.size()
              << " FreeDOF=" << freeDof << std::endl;

    // (2) モード解析
    int target_modes = 450;
    int nev = std::min(target_modes, freeDof - 2);
    std::vector<double> eigs;
    std::vector<std::vector<Displacement>> modes;

    auto t2 = clock::now();
    int computed = model.SolveVibration(nev, eigs, modes);
    auto t3 = clock::now();
    int M = (int)modes.size();
    std::cout << "[2] SolveVibration    : " << sec(t2, t3) << " s   "
              << "requested=" << nev << " computed=" << computed
              << " modes_size=" << M << std::endl;

    if (M == 0) {
        std::cout << "    !! No modes computed, abort." << std::endl;
        return;
    }

    // (3) CQC本体（Compute = 3 x calculate_responseCQC）
    NotifiedSpectrum spectrum;
    auto model_ptr = std::make_shared<FEModel>(model);
    FEVibrateResult vibresult(model_ptr, modes, eigs);

    auto t4 = clock::now();
    ResponseSpectrumMethod rsm(
        model_ptr, vibresult, Vector(1, 0, 0),
        &spectrum, ResponseSpectrumMethodType::CQC);
    auto t5 = clock::now();
    std::cout << "[3] CQC ctor+Compute  : " << sec(t4, t5)
              << " s   (= 3 x calculate_responseCQC)" << std::endl;
    std::cout << "    per response type : " << sec(t4, t5) / 3.0 << " s" << std::endl;

    // (4) ModeVectors() の値返しコストを単独計測
    {
        const int iter = 5;
        auto a = clock::now();
        volatile size_t sink = 0;
        for (int k = 0; k < iter; ++k) {
            auto mv = vibresult.ModeVectors(); // 値返しコピー
            sink += mv.size();
        }
        auto b = clock::now();
        double per_call_ms = sec(a, b) * 1000.0 / iter;
        long long calls_per_cqc = (long long)M * M + (long long)M;
        double est_copy_s_per_cqc = per_call_ms / 1000.0 * (double)calls_per_cqc;
        std::cout << "[4] ModeVectors() cost: " << per_call_ms << " ms/call   "
                  << "(値返しでvector-of-vectorを丸ごとコピー)" << std::endl;
        std::cout << "    CQC内での呼出回数 = M*M + M = " << calls_per_cqc << std::endl;
        std::cout << "    推定コピー時間/1CQC = " << est_copy_s_per_cqc << " s" << std::endl;
        std::cout << "    Compute()全体での推定コピー時間 = " << est_copy_s_per_cqc * 3.0 << " s" << std::endl;
    }

    // (4b) SRSS / ABS の Compute 時間も測る（修正前は M 回の値コピーで遅かった）
    {
        auto a1 = clock::now();
        ResponseSpectrumMethod rsm_srss(
            model_ptr, vibresult, Vector(1, 0, 0),
            &spectrum, ResponseSpectrumMethodType::SRSS);
        auto a2 = clock::now();
        std::cout << "[4b] SRSS ctor+Compute: " << sec(a1, a2) << " s (= 3 x SRSS)" << std::endl;

        auto b1 = clock::now();
        ResponseSpectrumMethod rsm_abs(
            model_ptr, vibresult, Vector(1, 0, 0),
            &spectrum, ResponseSpectrumMethodType::ABS);
        auto b2 = clock::now();
        std::cout << "[4c] ABS  ctor+Compute: " << sec(b1, b2) << " s (= 3 x ABS)" << std::endl;
    }

    // (5) 最適化版（コピーなし＋対称性活用＋スペクトル1回計算）
    auto t6 = clock::now();
    std::vector<double> periods(M);
    for (int i = 0; i < M; ++i) periods[i] = 2.0 * PI / eigs[i];
    std::vector<double> part_facs = vibresult.ParticipationFactors(Vector(1, 0, 0));
    std::vector<double> sa(M);
    for (int i = 0; i < M; ++i) sa[i] = spectrum.Acceleration(periods[i]);

    double damping = 0.05;
    std::vector<Displacement> resp(N);
    for (int j = 0; j < M; ++j) {
        const auto& uj = modes[j]; // ★参照 (コピー無し)
        double diag_fac = sa[j] * sa[j] * part_facs[j] * part_facs[j];
        for (int i = 0; i < N; ++i) {
            const auto& dj = uj[i];
            resp[i] += Displacement(
                diag_fac * dj.Dx() * dj.Dx(), diag_fac * dj.Dy() * dj.Dy(), diag_fac * dj.Dz() * dj.Dz(),
                diag_fac * dj.Rx() * dj.Rx(), diag_fac * dj.Ry() * dj.Ry(), diag_fac * dj.Rz() * dj.Rz());
        }
        for (int k = j + 1; k < M; ++k) {
            const auto& uk = modes[k];
            double rjk = periods[k] / periods[j];
            double corr = 8.0 * damping * damping * (1.0 + rjk) * pow(rjk, 1.5) /
                          (pow(1.0 - rjk * rjk, 2.0) + 4.0 * damping * damping * rjk * pow(1.0 + rjk, 2.0));
            double fac = 2.0 * sa[j] * sa[k] * part_facs[j] * part_facs[k] * corr;
            for (int i = 0; i < N; ++i) {
                const auto& dj = uj[i];
                const auto& dk = uk[i];
                resp[i] += Displacement(
                    fac * dj.Dx() * dk.Dx(), fac * dj.Dy() * dk.Dy(), fac * dj.Dz() * dk.Dz(),
                    fac * dj.Rx() * dk.Rx(), fac * dj.Ry() * dk.Ry(), fac * dj.Rz() * dk.Rz());
            }
        }
    }
    for (int i = 0; i < N; ++i) {
        resp[i] = Displacement(
            sqrt(std::max(0.0, resp[i].Dx())), sqrt(std::max(0.0, resp[i].Dy())), sqrt(std::max(0.0, resp[i].Dz())),
            sqrt(std::max(0.0, resp[i].Rx())), sqrt(std::max(0.0, resp[i].Ry())), sqrt(std::max(0.0, resp[i].Rz())));
    }
    auto t7 = clock::now();
    std::cout << "[5] Optimized 1xCQC   : " << sec(t6, t7) << " s   "
              << "(参照化 + 上三角のみ + spectrum事前計算)" << std::endl;

    // 結果の整合性チェック (Z方向最大値の比較)
    double max_rsm = 0, max_opt = 0;
    auto rsm_disp = rsm.GetAccelerations();
    for (int i = 0; i < N; ++i) {
        if (rsm_disp[i].Dx() > max_rsm) max_rsm = rsm_disp[i].Dx();
        if (resp[i].Dx() > max_opt) max_opt = resp[i].Dx();
    }
    std::cout << "    sanity: rsm.maxAx=" << max_rsm
              << " opt.maxAx=" << max_opt << std::endl;
}

// パフォーマンス計測用の立体ラーメンモデル（nx×ny柱 × nz層）
FEModel FrameBuildingModel(int nx, int ny, int nz, double bay, double hstory) {
    FEModel model;
    int id = 0;
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                Node n(id, i * bay, j * bay, k * hstory);
                if (k == 0) n.Fix.FixAll();
                model.Nodes.push_back(n);
                id++;
            }

    Material m0(5000.0, 0.2);
    m0.dense = 5.0 / 1000.0 / 1000.0;
    Section s0(100, 833.33, 833.33, 1406.25);
    model.Materials.push_back(m0);
    model.Sections.push_back(s0);

    auto idx = [&](int i, int j, int k) { return (k * ny + j) * nx + i; };
    int eid = 0;
    // 柱
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                model.Elements.push_back(std::make_shared<BeamElement>(eid++,
                    &model.Nodes[idx(i, j, k)], &model.Nodes[idx(i, j, k + 1)],
                    &model.Sections[0], model.Materials[0]));
    // X方向梁
    for (int k = 1; k <= nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i + 1 < nx; i++)
                model.Elements.push_back(std::make_shared<BeamElement>(eid++,
                    &model.Nodes[idx(i, j, k)], &model.Nodes[idx(i + 1, j, k)],
                    &model.Sections[0], model.Materials[0]));
    // Y方向梁
    for (int k = 1; k <= nz; k++)
        for (int j = 0; j + 1 < ny; j++)
            for (int i = 0; i < nx; i++)
                model.Elements.push_back(std::make_shared<BeamElement>(eid++,
                    &model.Nodes[idx(i, j, k)], &model.Nodes[idx(i, j + 1, k)],
                    &model.Sections[0], model.Materials[0]));
    return model;
}

// SolveBuckling の規模×モード数スケーリング計測
void BenchBucklingScaling() {
    using clock = std::chrono::high_resolution_clock;
    auto sec = [](auto a, auto b) {
        return std::chrono::duration<double>(b - a).count();
    };

    std::cout << "\n=== BenchBucklingScaling ===" << std::endl;

    auto run = [&](const char* name, FEModel& model, std::initializer_list<int> nevs) {
        // 全自由節点に鉛直荷重（柱に軸圧縮を入れる）
        std::vector<std::shared_ptr<LoadBase>> loads;
        for (size_t i = 0; i < model.Nodes.size(); i++)
            if (!model.Nodes[i].Fix.IsAnyFix())
                loads.push_back(std::make_shared<NodeLoad>(NodeLoad((int)i, 0, 0, -1.0)));

        auto model_ptr = std::make_shared<FEModel>(model);
        FELinearStaticOp st(model_ptr, loads);
        auto t0 = clock::now();
        st.Compute();
        auto t1 = clock::now();
        std::cout << "\n[" << name << "]  Nodes=" << model.NodeNum()
                  << "  FreeDOF=" << model.FreeDOFNum()
                  << "  static=" << sec(t0, t1) << " s" << std::endl;

        auto st_ptr = std::make_shared<FELinearStaticOp>(st);
        for (int nev : nevs) {
            FEBucklingAnalysis ba(st_ptr);
            ba.mode_num = nev;
            try {
                auto a = clock::now();
                int r = ba.SolveBuckling();
                auto b = clock::now();
                std::cout << "  nev=" << nev << "  ret=" << r
                          << "  time=" << sec(a, b) << " s"
                          << "  lambda1=" << (ba.eigs.empty() ? 0.0 : ba.eigs[0]) << std::endl;
            }
            catch (const std::exception& e) {
                std::cout << "  nev=" << nev << "  EXCEPTION: " << e.what() << std::endl;
            }
        }
    };

    { FEModel m = FrameBuildingModel(8, 8, 15, 6000, 4000); run("Frame 8x8x15", m, { 5, 20, 50 }); }
    { FEModel m = FrameBuildingModel(10, 10, 20, 6000, 4000); run("Frame 10x10x20", m, { 5, 20, 50 }); }
    { FEModel m = FrameBuildingModel(12, 12, 25, 6000, 4000); run("Frame 12x12x25", m, { 5, 20 }); }
    { FEModel m = ShallowDomeModel(20000.0, 2000.0, 32); run("Dome ndiv=32", m, { 5, 20 }); }
}

// 座屈解析の検証用: 固有値と モード形状の代表値(符号非依存の max|φ|)を出力
void PrintBucklingCheck(FEModel& model, std::vector<std::shared_ptr<LoadBase>>& loads,
    int nev, const char* title) {
    try {
        auto model_ptr = std::make_shared<FEModel>(model);
        FELinearStaticOp st(model_ptr, loads);
        st.Compute();

        FEBucklingAnalysis ba(std::make_shared<FELinearStaticOp>(st));
        ba.mode_num = nev;
        int r = ba.SolveBuckling();
        std::cout << "\n--- BucklingCheck: " << title
                  << " (requested=" << nev << ", ret=" << r << ") ---"
                  << std::setprecision(10) << std::endl;
        for (size_t i = 0; i < ba.eigs.size(); i++) {
            double vmax = 0.0;
            for (const Displacement& d : ba.mode_vectors[i]) {
                vmax = std::max(vmax, std::abs(d.Dx()));
                vmax = std::max(vmax, std::abs(d.Dy()));
                vmax = std::max(vmax, std::abs(d.Dz()));
            }
            std::cout << "  mode " << i + 1 << ": lambda=" << ba.eigs[i]
                      << "  max|phi|=" << vmax << std::endl;
        }
        // 縮退ペア（λ重複）は固有空間内の基底が任意なので、基底回転に
        // 不変な合成振幅 max sqrt(φi^2 + φj^2) で固有空間の一致を確認する
        for (size_t i = 0; i + 1 < ba.eigs.size(); i++) {
            if (std::abs(ba.eigs[i + 1] - ba.eigs[i]) > 1e-8 * std::abs(ba.eigs[i]))
                continue;
            double pmax = 0.0;
            for (size_t n = 0; n < ba.mode_vectors[i].size(); n++) {
                const Displacement& a = ba.mode_vectors[i][n];
                const Displacement& b = ba.mode_vectors[i + 1][n];
                pmax = std::max(pmax, std::sqrt(a.Dx() * a.Dx() + b.Dx() * b.Dx()));
                pmax = std::max(pmax, std::sqrt(a.Dy() * a.Dy() + b.Dy() * b.Dy()));
                pmax = std::max(pmax, std::sqrt(a.Dz() * a.Dz() + b.Dz() * b.Dz()));
            }
            std::cout << "  pair(" << i + 1 << "," << i + 2
                      << "): max|phi_pair|=" << pmax << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cout << "\n--- BucklingCheck: " << title
                  << "  EXCEPTION: " << e.what() << " ---" << std::endl;
    }
}

void TestBucklingCheck() {
    // 頂部集中荷重の片持ち柱（オイラー座屈）
    {
        FEModel model = CantiColumnModel(2000, 10);
        std::vector<std::shared_ptr<LoadBase>> loads;
        loads.push_back(std::make_shared<NodeLoad>(
            NodeLoad((int)model.Nodes.size() - 1, 0, 0, -1.0)));
        PrintBucklingCheck(model, loads, 6, "CantiColumn(2000,10)");
    }
    // ピラミッドトラス（頂点の回転は固定: 未固定だと静的解析段階で
    // 剛性ゼロ自由度により全体Kが特異になり解けないため）
    {
        FEModel model = PyramidTrussModel(1000, 4, 100);
        model.Nodes[0].Fix = Support(false, false, false, true, true, true);
        std::vector<std::shared_ptr<LoadBase>> loads;
        loads.push_back(std::make_shared<NodeLoad>(NodeLoad(0, 0, 0, -1.0)));
        PrintBucklingCheck(model, loads, 1, "PyramidTruss(1000,4,100)");
    }
    // 小型立体ラーメン（全節点鉛直荷重）
    {
        FEModel model = FrameBuildingModel(4, 4, 5, 6000, 4000);
        std::vector<std::shared_ptr<LoadBase>> loads;
        for (size_t i = 0; i < model.Nodes.size(); i++)
            if (!model.Nodes[i].Fix.IsAnyFix())
                loads.push_back(std::make_shared<NodeLoad>(NodeLoad((int)i, 0, 0, -1.0)));
        PrintBucklingCheck(model, loads, 12, "Frame 4x4x5");
    }
    // ドーム板（全節点鉛直荷重）
    {
        FEModel model = ShallowDomeModel(20000.0, 2000.0, 20);
        std::vector<std::shared_ptr<LoadBase>> loads;
        for (size_t i = 0; i < model.Nodes.size(); i++)
            if (!model.Nodes[i].Fix.IsAnyFix())
                loads.push_back(std::make_shared<NodeLoad>(NodeLoad((int)i, 0, 0, -1.0)));
        PrintBucklingCheck(model, loads, 5, "Dome ndiv=20");
    }
}

// SolveVibration の規模×モード数スケーリング計測
void BenchVibrationScaling() {
    using clock = std::chrono::high_resolution_clock;
    auto sec = [](auto a, auto b) {
        return std::chrono::duration<double>(b - a).count();
    };

    std::cout << "\n=== BenchVibrationScaling ===" << std::endl;

    auto run = [&](const char* name, FEModel& model, std::initializer_list<int> nevs) {
        model.ComputeElementNodeMass();
        std::cout << "\n[" << name << "]  Nodes=" << model.NodeNum()
                  << "  FreeDOF=" << model.FreeDOFNum() << std::endl;
        for (int nev : nevs) {
            std::vector<double> eigs;
            std::vector<std::vector<Displacement>> modes;
            auto t0 = clock::now();
            int computed = model.SolveVibration(nev, eigs, modes);
            auto t1 = clock::now();
            std::cout << "  nev=" << nev << "  computed=" << computed
                      << "  time=" << sec(t0, t1) << " s"
                      << "  T1=" << (eigs.empty() ? 0.0 : 2 * PI / eigs[0]) << std::endl;
        }
    };

    { FEModel m = ShallowDomeModel(20000.0, 2000.0, 20); run("Dome ndiv=20", m, { 10, 50, 150, 450 }); }
    { FEModel m = ShallowDomeModel(20000.0, 2000.0, 32); run("Dome ndiv=32", m, { 10, 50, 150, 300 }); }
    { FEModel m = ShallowDomeModel(20000.0, 2000.0, 40); run("Dome ndiv=40", m, { 10, 150 }); }
    { FEModel m = FrameBuildingModel(8, 8, 15, 6000, 4000); run("Frame 8x8x15", m, { 10, 50, 150 }); }
    { FEModel m = FrameBuildingModel(10, 10, 20, 6000, 4000); run("Frame 10x10x20", m, { 10, 150 }); }
}

void TestWallNodalForces() {
    std::cout << "\n========== TestWallNodalForces ==========" << std::endl;

    // 鉛直壁: XZ平面(法線=+Y), 幅W(X方向), 高さH(Z方向)
    // 節点並び: n0=脚部左, n1=脚部右, n2=頂部右, n3=頂部左
    const double W = 2.0, H = 3.0, t = 0.2;
    Point p0(0, 0, 0), p1(W, 0, 0), p2(W, 0, H), p3(0, 0, H);

    Material m0(2.05e5, 0.3);
    Node n0(p0); n0.id = 0;
    Node n1(p1); n1.id = 1;
    Node n2(p2); n2.id = 2;
    Node n3(p3); n3.id = 3;
    QuadPlateElement pel(&n0, &n1, &n2, &n3, t, m0);

    FEModel model;
    model.Nodes.push_back(n0);
    model.Nodes.push_back(n1);
    model.Nodes.push_back(n2);
    model.Nodes.push_back(n3);
    model.Nodes[0].Fix.FixAll(); // 脚部固定（片持ち壁）
    model.Nodes[1].Fix.FixAll();
    model.add_element(pel);

    const double Fx = 10.0;  // 面内水平(X) → 壁の面内せん断
    const double Fz = -20.0; // 鉛直(Z)     → 壁の軸力
    std::vector<std::shared_ptr<LoadBase>> loads;
    loads.push_back(std::make_shared<NodeLoad>(NodeLoad(2, Fx / 2, 0, Fz / 2)));
    loads.push_back(std::make_shared<NodeLoad>(NodeLoad(3, Fx / 2, 0, Fz / 2)));

    FELinearStaticOp op(std::make_shared<FEModel>(model), loads);
    op.Compute();

    std::cout << "Wall W=" << W << " H=" << H << " t=" << t << std::endl;
    std::cout << "Top load total: Fx(in-plane shear)=" << Fx << ", Fz(axial)=" << Fz << std::endl;
    std::cout << "Expected at base center: sumFx=" << Fx << ", sumFz=" << Fz
              << ", overturning(My about +Y)= +Fx*H=" << (Fx * H)
              << "  (reaction side = negatives)" << std::endl;

    auto printNL = [](const char *tag, NodeLoadData &d) {
        std::cout << "  " << tag << " id=" << d.id
                  << "  F=(" << d.Px() << ", " << d.Py() << ", " << d.Pz() << ")"
                  << "  M=(" << d.Mx() << ", " << d.My() << ", " << d.Mz() << ")" << std::endl;
    };

    std::cout << "\n-- Reactions (GetReactForces) --" << std::endl;
    std::vector<NodeLoad> react = op.GetReactForces();
    for (auto &r : react)
        std::cout << "  node " << r.id
                  << "  F=(" << r.Px() << ", " << r.Py() << ", " << r.Pz() << ")"
                  << "  M=(" << r.Mx() << ", " << r.My() << ", " << r.Mz() << ")" << std::endl;

    std::cout << "\n-- GetPlateNodalForces local=false (global) --" << std::endl;
    std::vector<NodeLoadData> nfg = op.GetPlateNodalForces(0, false);
    for (auto &d : nfg) printNL("global", d);

    std::cout << "\n-- GetPlateNodalForces local=true (element plane axes) --" << std::endl;
    std::vector<NodeLoadData> nfl = op.GetPlateNodalForces(0, true);
    for (auto &d : nfl) printNL("local ", d);

    // 脚部(node0,1)の節点力を基底中心 c=(W/2,0,0) へ剛体換算（global）
    const double cx = W / 2.0, cy = 0.0, cz = 0.0;
    const double px[4] = {0, W, W, 0};
    const double py[4] = {0, 0, 0, 0};
    const double pz[4] = {0, 0, H, H};
    double Fsum[3] = {0, 0, 0}, Msum[3] = {0, 0, 0};
    for (auto &d : nfg) {
        if (d.id != 0 && d.id != 1) continue;
        double fx = d.Px(), fy = d.Py(), fz = d.Pz();
        double rx = px[d.id] - cx, ry = py[d.id] - cy, rz = pz[d.id] - cz;
        Fsum[0] += fx; Fsum[1] += fy; Fsum[2] += fz;
        Msum[0] += d.Mx() + (ry * fz - rz * fy);
        Msum[1] += d.My() + (rz * fx - rx * fz);
        Msum[2] += d.Mz() + (rx * fy - ry * fx);
    }
    std::cout << "\n-- Base reduction (sum node0,1 nodal forces about base center) --" << std::endl;
    std::cout << "  F_a(global)=(" << Fsum[0] << ", " << Fsum[1] << ", " << Fsum[2] << ")" << std::endl;
    std::cout << "  M_a(global)=(" << Msum[0] << ", " << Msum[1] << ", " << Msum[2] << ")" << std::endl;
    std::cout << "  => in-plane shear ~ F_a.x=" << Fsum[0]
              << ", axial ~ F_a.z=" << Fsum[2]
              << ", overturning ~ M_a.y=" << Msum[1] << std::endl;

    // ---- 部材フレームへ射影（z=部材軸(鉛直), x=面外法線, y=面内水平）----
    auto dot = [](const double a[3], const double b[3]) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; };
    const double ez[3] = {0, 0, 1};  // 部材軸(鉛直)
    const double ex[3] = {0, 1, 0};  // 面外法線
    const double ey[3] = {-1, 0, 0}; // 面内水平 (ez×ex)

    // 脚部a: 反力側のため符号反転 → 連続断面力
    {
        double F[3] = {-Fsum[0], -Fsum[1], -Fsum[2]};
        double M[3] = {-Msum[0], -Msum[1], -Msum[2]};
        std::cout << "\n-- Section force @ base a (z=axis, x=out-of-plane, y=in-plane) --" << std::endl;
        std::cout << "  N=" << dot(F, ez) << "  Qx(out)=" << dot(F, ex) << "  Qy(in)=" << dot(F, ey)
                  << "  Mx(in-plane/overturn)=" << dot(M, ex) << "  My(out)=" << dot(M, ey)
                  << "  Mz(torsion)=" << dot(M, ez) << std::endl;
    }
    // 頂部b: 載荷側のためそのまま
    {
        const double bx = W / 2.0, by = 0.0, bz = H;
        double F[3] = {0, 0, 0}, M[3] = {0, 0, 0};
        for (auto &d : nfg) {
            if (d.id != 2 && d.id != 3) continue;
            double fx = d.Px(), fy = d.Py(), fz = d.Pz();
            double rx = px[d.id] - bx, ry = py[d.id] - by, rz = pz[d.id] - bz;
            F[0] += fx; F[1] += fy; F[2] += fz;
            M[0] += d.Mx() + (ry * fz - rz * fy);
            M[1] += d.My() + (rz * fx - rx * fz);
            M[2] += d.Mz() + (rx * fy - ry * fx);
        }
        std::cout << "-- Section force @ top b --" << std::endl;
        std::cout << "  N=" << dot(F, ez) << "  Qx(out)=" << dot(F, ex) << "  Qy(in)=" << dot(F, ey)
                  << "  Mx(in-plane/overturn)=" << dot(M, ex) << "  My(out)=" << dot(M, ey)
                  << "  Mz(torsion)=" << dot(M, ez) << std::endl;
    }
    std::cout << "==========================================" << std::endl;
}

int main(void) {
    // 壁→部材 節点力(手法1)の検証: TestWallNodalForces();
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
    //CheckCantiPyramidTrussBuckling(1000, 4, 100);
    //CheckQuadPlateBuckling();

    //BenchResponseSpectrumCQC();

    // Test mergeMatrixWithResize
    // TestMergeMatrixWithResize();

    // CheckSparseSolver();

    // 250929_Debug
    // TestSimplaFrame();

    //CheckQuadElement1();

    //TestBeamSemiRigidMethod();

    // 260103 Debug - SparseMatrixUtils のテスト
    // TestSparseMatrixUtils();

    // 固有値解析の前後比較検証
    TestVibrationCheck();
    TestVibrationCheckTruss();

    // 剛体連結を用いた剛床サンプル（末尾で固有値解析チェックも実行）
    TestRigidFloorWithCenterMaster_Shuffled();

    // 多数モード（450本）のベンチマーク
    BenchResponseSpectrumCQC();

    // 座屈解析の前後比較検証
    TestBucklingCheck();

    // 実行速度計測（規模×モード数）
    //BenchVibrationScaling();
    //BenchBucklingScaling();

    // 260713 Debug - 減衰初期化子(質量比例・レイリー)のテスト
    //TestDampInitializers();

}

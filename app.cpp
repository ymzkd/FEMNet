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

	Material m0(5000, 0.2);
	m0.dense = 5.0 / 1000 / 1000;

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

int main(void) {
	//std::cout << "TestMethod1 Start" << std::endl;
	//TestMethod1();
	
	//TestBeamTorsionMethod1();

	//CheckCantiBeamVibration();
	CheckCantiBeamVibration2();

	//CheckQuadElement1();

	//TestBeamSemiRigidMethod();

}

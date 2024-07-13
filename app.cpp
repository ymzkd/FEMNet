#include<iostream>

#include "Element.h"
#include "Model.h"

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
	TriPlaneElement pel(&n0, &n1, &n2, 12.0, &m0);

	// std::cout << "TransMat: \n" << pel.trans_matrix() << std::endl;
	std::cout << "StiffMat: \n" << pel.StiffnessMatrix() << std::endl;
}

// Tri Plane Element
void TestMethod2() {
	Node n0(0, 0, 0); n0.id = 0;
	Node n1(1000, 0, 0); n1.id = 1;
	Node n2(1000, 1000, 0); n2.id = 2;
	Node n3(0, 1000, 0); n3.id = 3;
	n0.Fix.FixAll();
	n1.Fix.FixAll();
	n2.Fix.FixAll();
	n3.Fix.FixAll();
	n2.Fix.Ux() = false;
	n2.Fix.Uy() = false;
	n3.Fix.Ux() = false;
	n3.Fix.Uy() = false;

	Material m0(5000, 0.2);
	
	TriPlaneElement tpe0(&n0, &n1, &n2, 12, &m0);
	TriPlaneElement tpe1(&n2, &n3, &n0, 12, &m0);
	std::cout << "area1: " << tpe0.Area() << std::endl;
	std::cout << "area2: " << tpe1.Area() << std::endl;

	SSModel model;
	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);
	model.Nodes.push_back(n3);

	model.add_element(&tpe0);
	model.add_element(&tpe1);

	std::list<Load> loads;
	loads.push_back(Load(2, 1000, 0, 0));
	loads.push_back(Load(3, 1000, 0, 0));

	std::vector<Displacement> disp = model.Solve(loads);
	for (size_t i = 0; i < disp.size(); i++)
	{
		std::cout << "id: " << i << std::endl;
		std::cout << disp[i] << std::endl;
	}

	MembraneStressData ms = tpe0.stress(disp[0], disp[1], disp[2]);

	std::cout << "OUTPUT Plane STRESS" << std::endl;
	std::cout << ms << std::endl;
}

void TestMethod3() {

	SSModel model;

	Node n0(0, 0, 0); n0.id = 0;
	//Node n1(0, 0, 1000); n1.id = 1;
	Node n1(1000, 0, 0); n1.id = 1;
	n0.Fix.FixAll();

	Material m0(5000, 0.2);
	Section s0(600, 45000, 20000, 10000);

	BeamElement b0(&n0, &n1, &s0, &m0);
	ElementBase* eb = &b0;

	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.add_element(&b0);

	std::list<Load> loads;
	loads.push_back(Load(1, 0, 10, -10, 0, 0, 0));

	std::vector<Displacement> disp = model.Solve(loads);
	for (size_t i = 0; i < disp.size(); i++)
	{
		std::cout << "id: " << i << std::endl;
		std::cout << disp[i] << std::endl;
	}

	BeamStress bs = b0.stress(disp[0], disp[1]);

	std::cout << "OUTPUT BEAM STRESS" << std::endl;
	std::cout << bs << std::endl;

}

// Tri Plate
void TestMethod4() {
	Node n0(0, 0, 0); n0.id = 0;
	Node n1(1000, 0, 0); n1.id = 1;
	Node n2(1000, 1000, 0); n2.id = 2;
	Node n3(0, 1000, 0); n3.id = 3;
	Node n4(500, 500, 0); n4.id = 4;
	n0.Fix.PinFix();
	n1.Fix.PinFix();
	n2.Fix.PinFix();
	n3.Fix.PinFix();

	//n0.Fix.Rz() = true;
	//n1.Fix.Rz() = true;
	//n2.Fix.Rz() = true;
	//n3.Fix.Rz() = true;
	//n4.Fix.Rz() = true;

	Material m0(5000, 0.2);

	TriPlateElement tpe0(&n0, &n1, &n4, 12, &m0);
	TriPlateElement tpe1(&n1, &n2, &n4, 12, &m0);
	TriPlateElement tpe2(&n2, &n3, &n4, 12, &m0);
	TriPlateElement tpe3(&n3, &n0, &n4, 12, &m0);
	//std::cout << "tpe0 Stiff Matrix: \n" << tpe3.StiffnessMatrix() << std::endl;
   // std::cout << "Sym?: " << tpe0.StiffnessMatrix().isApprox(tpe0.StiffnessMatrix().transpose(), 0.000001);
   // std::cout << "area2: " << tpe1.Area() << std::endl;

	SSModel model;
	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);
	model.Nodes.push_back(n3);
	model.Nodes.push_back(n4);

	model.add_element(&tpe0);
	model.add_element(&tpe1);
	model.add_element(&tpe2);
	model.add_element(&tpe3);

	std::list<Load> loads;
	loads.push_back(Load(4, 0, 0, -1000));

	std::vector<Displacement> disp = model.Solve(loads);
	for (size_t i = 0; i < disp.size(); i++)
	{
		std::cout << "id: " << i << std::endl;
		std::cout << disp[i] << std::endl;
	}

	PlateStressData ms = tpe0.stress(disp[0], disp[1], disp[4], 0, 0);

	std::cout << "OUTPUT Plate STRESS" << std::endl;
	std::cout << ms << std::endl;

	std::cout << "check shear stress" << std::endl;
	// tpe0.shearstress(disp[0], disp[1], disp[4], 0, 0);
	tpe0.shearstress(disp[0], disp[1], disp[4]);

}

// Quad Plane
void TestMethod5() {
	Node n0(0, 0, 0); n0.id = 0;
	Node n1(1000, 0, 0); n1.id = 1;
	Node n2(1000, 1000, 0); n2.id = 2;
	Node n3(0, 1000, 0); n3.id = 3;
	Node n4(500, 500, 0); n4.id = 4;
	n0.Fix.FixAll();
	n1.Fix.FixAll();
	n2.Fix.FixAll();
	n3.Fix.FixAll();
	n2.Fix.Ux() = false; n2.Fix.Uy() = false;
	n3.Fix.Ux() = false; n3.Fix.Uy() = false;
	//n3.Fix.Uz() = true; n3.Fix.Rz() = true;
	//n3.Fix.PinFix();

	//n0.Fix.Rz() = true;
	//n1.Fix.Rz() = true;
	//n2.Fix.Rz() = true;
	//n3.Fix.Rz() = true;
	//n4.Fix.Rz() = true;
	
	Material m0(5000, 0.2);

	QuadPlaneElement tpe0(&n0, &n1, &n2, &n3, 12, &m0);
	//TriPlateElement tpe0(&n0, &n1, &n4, 12, &m0);
	//TriPlateElement tpe1(&n1, &n2, &n4, 12, &m0);
	//TriPlateElement tpe2(&n2, &n3, &n4, 12, &m0);
	//TriPlateElement tpe3(&n3, &n0, &n4, 12, &m0);
	// std::cout << "tpe0 Stiff Matrix: \n" << tpe3.StiffnessMatrix() << std::endl;
	// std::cout << "Sym?: " << tpe0.StiffnessMatrix().isApprox(tpe0.StiffnessMatrix().transpose(), 0.000001);
	// std::cout << "area2: " << tpe1.Area() << std::endl;

	SSModel model;
	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);
	model.Nodes.push_back(n3);
	//model.Nodes.push_back(n4);

	model.add_element(&tpe0);
	//model.add_element(&tpe1);
	//model.add_element(&tpe2);
	//model.add_element(&tpe3);

	std::list<Load> loads;
	loads.push_back(Load(2, 1000, 0, 0));
	loads.push_back(Load(3, 1000, 0, 0));

	std::vector<Displacement> disp = model.Solve(loads);
	for (size_t i = 0; i < disp.size(); i++)
	{
		std::cout << "id: " << i << std::endl;
		std::cout << disp[i] << std::endl;
	}

	MembraneStressData ms = tpe0.stress(disp[0], disp[1], disp[2], disp[3], 0.5, 0.5);
	
	std::cout << "OUTPUT Plane STRESS" << std::endl;
	std::cout << ms << std::endl;

	//std::cout << "check shear stress" << std::endl;
	// tpe0.shearstress(disp[0], disp[1], disp[4], 0, 0);
	//tpe0.shearstress(disp[0], disp[1], disp[4]);
}


// Quad Plate
void TestMethodQuadPlate() {
	Node n0(0, 0, 0); n0.id = 0;
	Node n1(500, 0, 0); n1.id = 1;
	Node n2(1000, 0, 0); n2.id = 2;
	Node n3(0, 500, 0); n3.id = 3;
	Node n4(500, 500, 0); n4.id = 4;
	Node n5(1000, 500, 0); n5.id = 5;
	Node n6(0, 1000, 0); n6.id = 6;
	Node n7(500, 1000, 0); n7.id = 7;
	Node n8(1000, 1000, 0); n8.id = 8;
	
	n0.Fix.PinFix();
	n2.Fix.PinFix();
	n6.Fix.PinFix();
	n8.Fix.PinFix();
	
	//n2.Fix.Ux() = false; n2.Fix.Uy() = false;
	//n3.Fix.Ux() = false; n3.Fix.Uy() = false;
	//n3.Fix.Uz() = true; n3.Fix.Rz() = true;
	//n3.Fix.PinFix();

	//n0.Fix.Rz() = true;
	//n1.Fix.Rz() = true;
	//n2.Fix.Rz() = true;
	//n3.Fix.Rz() = true;
	//n4.Fix.Rz() = true;

	Material m0(5000, 0.2);

	QuadPlateElement tpe0(&n0, &n1, &n4, &n3, 12, &m0);
	QuadPlateElement tpe1(&n1, &n2, &n5, &n4, 12, &m0);
	QuadPlateElement tpe2(&n3, &n4, &n7, &n6, 12, &m0);
	QuadPlateElement tpe3(&n4, &n5, &n8, &n7, 12, &m0);

	// std::cout << "tpe0 Stiff Matrix: \n" << tpe3.StiffnessMatrix() << std::endl;
	// std::cout << "Sym?: " << tpe0.StiffnessMatrix().isApprox(tpe0.StiffnessMatrix().transpose(), 0.000001);
	// std::cout << "area2: " << tpe1.Area() << std::endl;

	SSModel model;
	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);
	model.Nodes.push_back(n3);
	model.Nodes.push_back(n4);
	model.Nodes.push_back(n5);
	model.Nodes.push_back(n6);
	model.Nodes.push_back(n7);
	model.Nodes.push_back(n8);

	model.add_element(&tpe0);
	model.add_element(&tpe1);
	model.add_element(&tpe2);
	model.add_element(&tpe3);

	std::list<Load> loads;
	loads.push_back(Load(4, 0, 0, -1000));

	std::vector<Displacement> disp = model.Solve(loads);
	for (size_t i = 0; i < disp.size(); i++)
	{
		std::cout << "id: " << i << std::endl;
		std::cout << disp[i] << std::endl;
	}

	PlateStressData ms = tpe0.stress(disp[0], disp[1], disp[4], disp[3], 1.0, 1.0);

	std::cout << "OUTPUT Plate STRESS" << std::endl;
	std::cout << ms << std::endl;


	std::cout << "Shear Stress" << std::endl;
	//tpe0.shearstress(disp[0], disp[1], disp[4], disp[3], 0.0, 0.0);
	tpe0.shearstress(disp[0], disp[1], disp[4], disp[3], 1.0, 1.0);
	//std::cout << "check shear stress" << std::endl;
	// tpe0.shearstress(disp[0], disp[1], disp[4], 0, 0);
	//tpe0.shearstress(disp[0], disp[1], disp[4]);
}

// Truss Test
void TestTrussMethod() {

	SSModel model;

	Node n0(0, 0, 0); n0.id = 0;
	Node n1(1000, 0, 0); n1.id = 1;
	Node n2(1000, 1000, 0); n2.id = 2;
	Node n3(0, 1000, 0); n3.id = 3;
	Node n4(500, 500, 100); n4.id = 4;
	n0.Fix.FixAll();
	n1.Fix.FixAll();
	n2.Fix.FixAll();
	n3.Fix.FixAll();
	n4.Fix.FixAll();
	n4.Fix.Ux() = false;
	n4.Fix.Uy() = false;
	n4.Fix.Uz() = false;

	Material m0(5000, 0.2);
	// 20 x 30
	Section s0(600, 45000, 20000, 10000);

	TrussElement b0(&n0, &n4, &s0, &m0);
	TrussElement b1(&n1, &n4, &s0, &m0);
	TrussElement b2(&n2, &n4, &s0, &m0);
	TrussElement b3(&n3, &n4, &s0, &m0);
	ElementBase* eb = &b0;

	model.Nodes.push_back(n0);
	model.Nodes.push_back(n1);
	model.Nodes.push_back(n2);
	model.Nodes.push_back(n3);
	model.Nodes.push_back(n4);
	
	model.add_element(&b0);
	model.add_element(&b1);
	model.add_element(&b2);
	model.add_element(&b3);

	std::list<Load> loads;
	loads.push_back(Load(4, 0, 0, -1000, 0, 0, 0));

	std::vector<Displacement> disp = model.Solve(loads);
	for (size_t i = 0; i < disp.size(); i++)
	{
		std::cout << "id: " << i << std::endl;
		std::cout << disp[i] << std::endl;
	}

	BeamStress bs = b0.stress(disp[0], disp[4]);

	std::cout << "OUTPUT Truss STRESS" << std::endl;
	std::cout << bs << std::endl;

}

int main(void) {
	//TestMethod1();
	
	//TestMethod4();
	
	//TestMethod5();
	
	// TestMethodQuadPlate();
	
	TestTrussMethod();
	
	//TestMethod2();
}

#include <numeric>

#include "LoadComponent.h"


double BeamPolyLoad::R0() {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.R0();
        });
}

double BeamPolyLoad::RA() {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.RA();
        });
}

double BeamPolyLoad::M0() {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.M0();
        });
}

double BeamPolyLoad::MA() {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.MA();
        });
}

double BeamPolyLoad::shear_force(double x) {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [&x](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.shear_force(x);
        });
}

double BeamPolyLoad::bending_moment(double x) {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [&x](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.bending_moment(x);
        });
}

double BeamPolyLoad::deflection(double x, double EI) {
    return std::accumulate(traps.begin(), traps.end(), 0.0, [&x, &EI](double sum, BeamTrapezoidalLoad trap) {
        return sum + trap.deflection(x, EI);
        });
}

NodeLoadData BeamPolyLoad::load_i()
{
    Eigen::Matrix3d tr = trans_matrix3(
        element->Nodes[0]->Location, element->Nodes[1]->Location, element->Beta).transpose();

    Eigen::Vector3d m, q;
    if (axis == BeamLoadAxis::ZAxis) {
        m = tr * Eigen::Vector3d(0, -M0(), 0);
		q = tr * Eigen::Vector3d(0, 0, R0());
	}
	else // YAxis
	{
		m = tr * Eigen::Vector3d(0, 0, M0());
		q = tr * Eigen::Vector3d(0, R0(), 0);
    }

    return NodeLoadData(element->Nodes[0]->id, q[0], q[1], q[2], m[0], m[1], m[2]);
}

NodeLoadData BeamPolyLoad::load_j()
{
    Eigen::Matrix3d tr = trans_matrix3(
        element->Nodes[0]->Location, element->Nodes[1]->Location, element->Beta).transpose();

    Eigen::Vector3d m, q;
    if (axis == BeamLoadAxis::ZAxis) {
        m = tr * Eigen::Vector3d(0, MA(), 0);
        q = tr * Eigen::Vector3d(0, 0, RA());
    }
    else // YAxis
    {
        m = tr * Eigen::Vector3d(0, 0, -MA());
        q = tr * Eigen::Vector3d(0, RA(), 0);
    }

    return NodeLoadData(element->Nodes[1]->id, q[0], q[1], q[2], m[0], m[1], m[2]);
}

BeamStressData BeamPolyLoad::GetBeamStress(double p)
{
    BeamStressData bsd;
    if (axis == BeamLoadAxis::ZAxis) {
	    bsd.My = bending_moment(p * element->length());
        bsd.Qz = shear_force(p * element->length());
	}
	else {
        bsd.Mz = bending_moment(p * element->length());
        bsd.Qy = shear_force(p * element->length());
	}
    return bsd;
}

Displacement BeamPolyLoad::GetDisplacement(double p)
{
    double I = (axis == BeamLoadAxis::ZAxis) ? element->Sec->Iy : element->Sec->Iz;
    double d = deflection(p * element->length(), element->Mat.Young * I);
    Eigen::Matrix3d tr = trans_matrix3(
        element->Nodes[0]->Location,
        element->Nodes[1]->Location, element->Beta);

    int vid = (axis == BeamLoadAxis::ZAxis) ? 2 : 1;
    Eigen::Vector3d dt = tr.row(vid) * d;
    return Displacement(dt[0], dt[1], dt[2], 0, 0, 0);
}

double AxialPolyLoad::axial_force(double x)
{
    return std::accumulate(traps.begin(), traps.end(), 0.0, [&x](double sum, AxialTrapezoidalLoad trap) {
        return sum + trap.axial_force(x);
        });
}
double AxialPolyLoad::N0()
{
    return std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, AxialTrapezoidalLoad trap) {
        return sum + trap.n1();
        });
}
double AxialPolyLoad::N1()
{
    return -std::accumulate(traps.begin(), traps.end(), 0.0, [](double sum, AxialTrapezoidalLoad trap) {
        return sum + trap.n2();
        });
}

NodeLoadData AxialPolyLoad::load_i()
{
    Vector dir = element->Nodes[1]->Location - element->Nodes[0]->Location;
    dir = dir * (N0() / dir.norm());
    return NodeLoadData(element->Nodes[0]->id, dir.x, dir.y, dir.z);
}

NodeLoadData AxialPolyLoad::load_j()
{
    Vector dir = element->Nodes[1]->Location - element->Nodes[0]->Location;
    dir = dir * (N1() / dir.norm());
    return NodeLoadData(element->Nodes[1]->id, dir.x, dir.y, dir.z);
}

BeamStressData AxialPolyLoad::GetBeamStress(double p)
{
    BeamStressData bsd;
    bsd.Nx = axial_force(p * element->length());
    return bsd;
}

Displacement AxialPolyLoad::GetDisplacement(double p)
{
    return Displacement();
}

NodeLoadData::NodeLoadData(int _id, double px, double py, double pz)
{
    id = _id;
    loads[0] = px;
    loads[1] = py;
    loads[2] = pz;
    loads[3] = 0;
    loads[4] = 0;
    loads[5] = 0;
}

NodeLoadData::NodeLoadData(int _id, double px, double py, double pz, double mx, double my, double mz)
{
    id = _id;
    loads[0] = px;
    loads[1] = py;
    loads[2] = pz;
    loads[3] = mx;
    loads[4] = my;
    loads[5] = mz;
}

NodeLoad::NodeLoad(int _id, double px, double py, double pz)
{
    id = _id;
    data = NodeLoadData(_id, px, py, pz);
}

NodeLoad::NodeLoad(int _id, double px, double py, double pz, double mx, double my, double mz)
{
    id = _id;
    data = NodeLoadData(_id, px, py, pz, mx, my, mz);
}

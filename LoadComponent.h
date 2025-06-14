#ifndef _LOAD_COMPONENT_
#define _LOAD_COMPONENT_

#ifndef SWIGCSHARP
#include<iostream>
#include <vector>
//#include <numeric>
#endif

#include "Components.h"
#include "Element.h"

class LoadBase {
public:
    LoadBase(){};
    
    virtual std::vector<NodeLoadData> NodeLoads() = 0;
};

class InertialForce: public LoadBase {
public:
	Vector accels;
	InertialForce() : InertialForce(0, 0, 0) {};
	InertialForce(double x, double y, double z) : accels(x, y, z) {};
	InertialForce(Vector v) : accels(v) {};
	
	std::vector<NodeLoadData> NodeLoads() override {
		return std::vector<NodeLoadData>();
	}
};

class NodeLoad : public LoadBase {
public:
    NodeLoadData data;
public:
    NodeLoad() : NodeLoad(-1, 0, 0, 0, 0, 0, 0) {};
    NodeLoad(int _id, double px, double py, double pz);
    NodeLoad(int _id, double px, double py, double pz, double mx, double my, double mz);

    int id = -1;
    double* loads() { return data.loads; };
    double& Px() { return data.loads[0]; }
    double& Py() { return data.loads[1]; }
    double& Pz() { return data.loads[2]; }
    double& Mx() { return data.loads[3]; }
    double& My() { return data.loads[4]; }
    double& Mz() { return data.loads[5]; }

    std::vector<NodeLoadData> NodeLoads() override { return { data }; }

    // operator<<
    friend std::ostream &operator<<(std::ostream &os, const NodeLoad &nodeLoad)
    {
        os << "NodeLoad ID: " << nodeLoad.id << "\n";
        os << "Px: " << nodeLoad.data.Px() << ", Py: " << nodeLoad.data.Py() << ", Pz: " << nodeLoad.data.Pz() << "\n";
        os << "Mx: " << nodeLoad.data.Mx() << ", My: " << nodeLoad.data.My() << ", Mz: " << nodeLoad.data.Mz() << "\n";
        return os;
    }
};


// 面荷重(WIP)
class PlateLoad : public LoadBase {
public:
	PlaneElementBase* element;
	std::vector<Vector> load_vecs;

    PlateLoad() {};
    
    PlateLoad(ElementBase* el, double px, double py, double pz) {

        element = dynamic_cast<PlaneElementBase*>(el);
        //if (!element) {
        //    throw std::runtime_error("Failed to cast ElementBase* to PlaneElementBase*");
        //}
        load_vecs.push_back(Vector(px, py, pz));
        load_vecs.push_back(Vector(px, py, pz));
        load_vecs.push_back(Vector(px, py, pz));
        if (element->NodeNum() == 4)
            load_vecs.push_back(Vector(px, py, pz));

    }

    PlateLoad(PlaneElementBase* el, double px, double py, double pz){ 
		element = el;
        load_vecs.push_back(Vector(px, py, pz));
        load_vecs.push_back(Vector(px, py, pz));
        load_vecs.push_back(Vector(px, py, pz));
        if (element->NodeNum() == 4)
            load_vecs.push_back(Vector(px, py, pz));
        
    }
	
    PlateLoad(PlaneElementBase* el, Vector v) {
		element = el;
		load_vecs.push_back(v);
		load_vecs.push_back(v);
		load_vecs.push_back(v);

        if (element->NodeNum() == 4)
            load_vecs.push_back(v);
	}
	
    PlateLoad(PlaneElementBase* el, Vector p1, Vector p2, Vector p3) {
		element = el;
		load_vecs.push_back(p1);
		load_vecs.push_back(p2);
		load_vecs.push_back(p3);
	}

    PlateLoad(PlaneElementBase* el, Vector p1, Vector p2, Vector p3, Vector p4) {
        element = el;
        load_vecs.push_back(p1);
        load_vecs.push_back(p2);
        load_vecs.push_back(p3);
        load_vecs.push_back(p4);
    }

    std::vector<NodeLoadData> NodeLoads() override;
};


enum BeamLoadAxis
{
    YAxis, ZAxis, XAxis
};

class BeamLoadBase : public LoadBase {
public:
    BeamElement* element;
    BeamLoadAxis axis;
    BeamLoadBase(BeamElement* element, BeamLoadAxis axis)
        : element(element), axis(axis) {}

    virtual NodeLoadData load_i() = 0;
    virtual NodeLoadData load_j() = 0;

    virtual BeamStressData GetBeamStress(double p) = 0;
    virtual Displacement GetDisplacement(double p) = 0;
};

// 梁の台形分布荷重
class BeamTrapezoidalLoad {
public:
    double w1, w2, L1, L2, L3, L;

    BeamTrapezoidalLoad(double w1, double w2, double L1, double L2, double L3)
        : w1(w1), w2(w2), L1(L1), L2(L2), L3(L3) {
        L = L1 + L2 + L3;
    }

    double R0() {
        double R0_EQ = (w1 * L2 / 2) * (2 * L3 / L + L2 / L - (L1 / L - L3 / L) * (2 * L1 * L3 / (L * L) + L2 * L3 / (L * L) + L1 * L2 / (L * L)));
        double R0_TR = ((w2 - w1) * L2 / 6) * (-3.0 / 5 * std::pow(L2, 3) / std::pow(L, 3) + 3.0 / 2 * std::pow(L2, 2) / std::pow(L, 2) * (1 - 2 * L3 / L) + 6 * L2 * L3 / std::pow(L, 2) * (1 - L3 / L) + 3 * std::pow(L3, 2) / std::pow(L, 2) * (3 - 2 * L3 / L));
        return R0_EQ + R0_TR;
    }

    double RA() {
        return L2 * (w1 + w2) / 2 - R0();
    }

    double M0() {
        double M0_EQ = (w1 * L2 * L / 8) * (std::pow((L2 / L + 2 * L3 / L), 2) * (2 * L1 / L + L2 / L) + 1.0 / 3 * std::pow((L2 / L), 2) * (2 - 6 * L3 / L - 3 * L2 / L));
        double M0_TR = ((w2 - w1) * L2 / 6) * ((std::pow((3 * L3 + L2), 2) / L / 3 + std::pow(L2, 2) / 6 / L - std::pow((3 * L3 + L2), 3) / (9 * std::pow(L, 2)) - 17.0 / 90 * std::pow(L2, 3) / std::pow(L, 2) - std::pow(L2, 2) * L3 / 2 / std::pow(L, 2)));
        return M0_EQ + M0_TR;
    }

    double MA() {
        double MA_EQ = (w1 * L * L2 / 8) * (std::pow((2 * L1 / L + L2 / L), 2) * (L2 / L + 2 * L3 / L) + 1.0 / 3 * std::pow((L2 / L), 2) * (2 - 6 * L1 / L - 3 * L2 / L));
        double MA_TR = ((w2 - w1) * L2 / 6) * (1.0 / 9 * std::pow((3 * L3 + L2), 3) / std::pow(L, 2) + 17.0 / 90 * std::pow(L2, 3) / std::pow(L, 2) + std::pow(L2, 2) * L3 / 2 / std::pow(L, 2) - 2 * std::pow((3 * L3 + L2), 2) / 3 / L - std::pow(L2, 2) / 3 / L + 3 * L3 + L2);
        return MA_EQ + MA_TR;
    }

    double shear_force(double x) {
        double R0 = this->R0();
        double RA = this->RA();
        double S;

        if (x < this->L1) {
            S = R0;
        }
        else if (this->L1 <= x && x <= this->L1 + this->L2) {
            S = R0 - (this->w1 * (x - this->L1) + ((this->w2 - this->w1) / 2 / this->L2) * pow((x - this->L1), 2));
        }
        else {
            S = -RA;
        }

        return S;
    }

    double bending_moment(double x) {
        double r0 = R0();
        double rA = RA();
        double m0 = M0();
        double mA = MA();
        double M;

        if (x < L1) {
            M = r0 * x - m0;
        }
        else if (L1 <= x && x <= L1 + L2) {
            M = r0 * x - m0 - (w1 / 2 * pow(x - L1, 2) + (w2 - w1) / 6 / L2 * pow(x - L1, 3));
        }
        else {
            M = rA * (L - x) - mA;
        }
        return M;
    }

    double deflection(double x, double EI) {
        double R0 = this->R0();
        double RA = this->RA();
        double M0 = this->M0();
        double MA = this->MA();
        double delta;

        if (x < L1) {
            delta = (1.0 / 6.0 / EI) * (3.0 * M0 * x * x - R0 * x * x * x);
        }
        else if (L1 <= x && x <= L1 + L2) {
            delta = (1.0 / 60.0 / EI) * (30.0 * M0 * x * x - 10.0 * R0 * x * x * x + ((w2 - w1) / 2.0 / L2) * pow(x - L1, 5) + 5.0 * w1 / 2.0 * pow(x - L1, 4));
        }
        else {
            delta = (1.0 / 6.0 / EI) * (3.0 * MA * pow(L - x, 2) - RA * pow(L - x, 3));
        }
        return delta;
    }

};

// 梁の多角形分布荷重
class BeamPolyLoad : public BeamLoadBase {
private:
    std::vector<BeamTrapezoidalLoad> traps;
public:
    std::vector<double> w;
    std::vector<double> params;
    
    // BeamElement* element;

    BeamPolyLoad(const std::vector<double> w, 
        const std::vector<double> params, 
        BeamElement* element, BeamLoadAxis axis)
            : w(w), params(params), BeamLoadBase(element, axis) {

        double length = element->length();

        // double length 
        for (int i = 0; i < w.size() - 1; i++) {
            double L1 = params[i] * length;
            double L2 = params[i + 1] * length - L1;
            double L3 = length - L1 - L2;
            traps.push_back(BeamTrapezoidalLoad(w[i], w[i + 1], L1, L2, L3));
        }
    }

    BeamPolyLoad(const std::vector<double> w,
        const std::vector<double> params, BeamLoadAxis axis)
        : w(w), params(params), BeamLoadBase(NULL, axis){}

    // i端反力
    double R0();
    // j端反力
    double RA();
    // i端曲げモーメント
    double M0();
    // j端曲げモーメント
    double MA();
    // せん断力
    double shear_force(double x);
    // 曲げモーメント
    double bending_moment(double x);
    // たわみ
    double deflection(double x, double EI);
    
    NodeLoadData load_i() override;
    NodeLoadData load_j() override;
    BeamStressData GetBeamStress(double p) override;
    Displacement GetDisplacement(double p) override;
    std::vector<NodeLoadData> NodeLoads() override {
        return { load_i(), load_j() };
    }

    static BeamPolyLoad CreateFromUniLoad(double w, BeamElement* element, BeamLoadAxis axis) {
        return BeamPolyLoad(std::vector<double>{w, w}, std::vector<double>{0, 1}, element, axis);
    }
};

class AxialTrapezoidalLoad {
public:
    double w1, w2, L1, L2, L3;
    AxialTrapezoidalLoad(double w1, double w2, double L1, double L2, double L3)
        : w1(w1), w2(w2), L1(L1), L2(L2), L3(L3) {
    }

    double n1() {
        double l = L1 + L2 + L3;
        double p = (w1 + w2) * L2 / 2.0;
        double dg = (w1 + 2 * w2) / 3.0 / (w1 + w2) * L2;
        return (L2 + L3 - dg) / l * p;
    }

    double n2() {
        double l = L1 + L2 + L3;
        double p = (w1 + w2) * L2 / 2.0;
        double dg = (w1 + 2 * w2) / 3.0 / (w1 + w2) * L2;
        return -(L1 + dg) / l * p;
    }

    double axial_force(double x) {
        double n1 = this->n1();
        double n2 = this->n2();

        if (x < L1)
            return n1;
        else if (x >= L1 && x <= L1 + L2)
            return (w1 - w2) / 2.0 / L2 * std::pow(x - L1, 2) - w1 * (x - L1) + n1;
        else
            return n2;
    }
};


// 軸方向の多角形分布荷重
class AxialPolyLoad : public BeamLoadBase {
private:
    std::vector<AxialTrapezoidalLoad> traps;
public:
    std::vector<double> params;
    std::vector<double> w;
    //BeamElement* element;

    AxialPolyLoad(const std::vector<double>& w, 
        const std::vector<double>& params, BeamElement* element)
        : w(w), params(params), BeamLoadBase(element, BeamLoadAxis::XAxis) {

        double length = element->length();

        for (int i = 0; i < w.size() - 1; i++) {
            double L1 = params[i] * length;
            double L2 = params[i + 1] * length - L1;
            double L3 = length - L1 - L2;
            traps.push_back(AxialTrapezoidalLoad(w[i], w[i + 1], L1, L2, L3));
        }
    }

    AxialPolyLoad(const std::vector<double>& w, 
		const std::vector<double>& params)
		: w(w), params(params), BeamLoadBase(NULL, BeamLoadAxis::XAxis) {}

    double axial_force(double x);
    double N0();
    double N1();

    NodeLoadData load_i() override;
    NodeLoadData load_j() override;

    std::vector<NodeLoadData> NodeLoads() override {
        return { load_i(), load_j() };
    }

    BeamStressData GetBeamStress(double p) override;
    Displacement GetDisplacement(double p) override;
};

class DynamicAccelLoad {
public:
    double timestep;
    Vector Direction;
	std::vector<double> Accels;

    size_t DataCount() {
		return Accels.size();
    }

	DynamicAccelLoad(double timestep, Vector Direction, std::vector<double> Accels)
		: timestep(timestep), Direction(Direction), Accels(Accels) {
	}
	DynamicAccelLoad(double timestep, double x, double y, double z, std::vector<double> Accels)
		: timestep(timestep), Direction(x, y, z), Accels(Accels) {
	}
	DynamicAccelLoad(double timestep, double x, double y, double z)
		: timestep(timestep), Direction(x, y, z) {
	}
};


#endif
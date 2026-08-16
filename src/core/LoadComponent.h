#ifndef _LOAD_COMPONENT_
#define _LOAD_COMPONENT_

#ifndef SWIG
#include<iostream>
#include <vector>
#endif

#include "Components.h"
#include "Elements/Elements.h"

enum class LoadType {
    None,
    BodyForce,
};

class LoadBase {
public:
    LoadBase(){};
    virtual ~LoadBase() = default;

    virtual LoadType Type() { return LoadType::None; }
    virtual std::vector<NodeLoadData> NodeLoads() = 0;

    // 荷重値を係数倍した複製を返す。
    // 荷重組み合わせを単一の静的荷重ケースへ合成する際に使用する。
    virtual std::shared_ptr<LoadBase> scaled(double factor) const = 0;
};

class InertialForce: public LoadBase {
public:

    Vector accels;
    InertialForce() : InertialForce(0, 0, 0) {};
    InertialForce(double x, double y, double z) : accels(x, y, z) {};
    InertialForce(Vector v) : accels(v) {};

    LoadType Type() override { return LoadType::BodyForce; }

    std::vector<NodeLoadData> NodeLoads() override {
        return std::vector<NodeLoadData>();
    }

    std::shared_ptr<LoadBase> scaled(double factor) const override {
        return std::make_shared<InertialForce>(
            accels.x * factor, accels.y * factor, accels.z * factor);
    }
};

enum class BodyForaceSelector : unsigned int {
    None = 0,
    NodeMass = 1 << 0,
    LoadMass = 1 << 1,
    ElementMass = 1 << 2,
    All = NodeMass | LoadMass | ElementMass
};

class NodeBodyForce : public LoadBase {
public:
    Node* node;
    Vector Accels;
    BodyForaceSelector selector = BodyForaceSelector::All;

    NodeBodyForce() : node(nullptr), Accels(0, 0, 0) {};
    NodeBodyForce(Node* node, double x, double y, double z) : node(node), Accels(x, y, z) {};

    LoadType Type() override { return LoadType::BodyForce; }

    std::vector<NodeLoadData> NodeLoads() override {
        Vector f = Accels * node->MassData.SumMass(); // ここで質量を掛けて加速度を計算
        return std::vector<NodeLoadData>{ NodeLoadData(node->id, f.x, f.y, f.z) };
    }

    std::shared_ptr<LoadBase> scaled(double factor) const override {
        auto copy = std::make_shared<NodeBodyForce>(*this);
        copy->Accels = copy->Accels * factor;
        return copy;
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

    std::shared_ptr<LoadBase> scaled(double factor) const override {
        auto copy = std::make_shared<NodeLoad>(*this);
        copy->data *= factor;
        return copy;
    }

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
    
    PlateLoad(ElementBase* el, double px, double py, double pz, bool local = false) {
        element = dynamic_cast<PlaneElementBase*>(el);
        Vector lv(px, py, pz);
        
        if (local) {
            Plane pl = element->plane;
            lv = lv.x * pl.ex + lv.y * pl.ey + lv.z * pl.ez;
        }
        
        load_vecs.push_back(Vector(lv.x, lv.y, lv.z));
        load_vecs.push_back(Vector(lv.x, lv.y, lv.z));
        load_vecs.push_back(Vector(lv.x, lv.y, lv.z));
        if (element->NodeNum() == 4)
            load_vecs.push_back(Vector(lv.x, lv.y, lv.z));
    }

    PlateLoad(PlaneElementBase* el, double px, double py, double pz, bool local = false) { 
        element = el;
        Vector lv(px, py, pz);
        
        if (local) {
            Plane pl = element->plane;
            lv = lv.x * pl.ex + lv.y * pl.ey + lv.z * pl.ez;
        }
        
        load_vecs.push_back(lv);
        load_vecs.push_back(lv);
        load_vecs.push_back(lv);
        if (element->NodeNum() == 4)
            load_vecs.push_back(lv);
    }
    
    PlateLoad(PlaneElementBase* el, Vector v, bool local = false) {
        element = el;
        Vector lv = v;
        if (local) {
            Plane pl = element->plane;
            lv = lv.x * pl.ex + lv.y * pl.ey + lv.z * pl.ez;
        }
        load_vecs.push_back(lv);
        load_vecs.push_back(lv);
        load_vecs.push_back(lv);
        if (element->NodeNum() == 4)
            load_vecs.push_back(lv);
    }
    
    PlateLoad(PlaneElementBase* el, Vector p1, Vector p2, Vector p3, bool local = false) {
        element = el;
        Vector lv1 = p1;
        Vector lv2 = p2;
        Vector lv3 = p3;

        if (local) {
            Plane pl = element->plane;
            lv1 = lv1.x * pl.ex + lv1.y * pl.ey + lv1.z * pl.ez;
            lv2 = lv2.x * pl.ex + lv2.y * pl.ey + lv2.z * pl.ez;
            lv3 = lv3.x * pl.ex + lv3.y * pl.ey + lv3.z * pl.ez;
        }

        load_vecs.push_back(lv1);
        load_vecs.push_back(lv2);
        load_vecs.push_back(lv3);
    }

    PlateLoad(PlaneElementBase* el, Vector p1, Vector p2, Vector p3, Vector p4, bool local = false) {
        element = el;
        Vector lv1 = p1;
        Vector lv2 = p2;
        Vector lv3 = p3;
        Vector lv4 = p4;

        if (local) {
            Plane pl = element->plane;
            lv1 = lv1.x * pl.ex + lv1.y * pl.ey + lv1.z * pl.ez;
            lv2 = lv2.x * pl.ex + lv2.y * pl.ey + lv2.z * pl.ez;
            lv3 = lv3.x * pl.ex + lv3.y * pl.ey + lv3.z * pl.ez;
            lv4 = lv4.x * pl.ex + lv4.y * pl.ey + lv4.z * pl.ez;
        }

        load_vecs.push_back(lv1);
        load_vecs.push_back(lv2);
        load_vecs.push_back(lv3);
        load_vecs.push_back(lv4);
    }

    std::vector<NodeLoadData> NodeLoads() override;

    std::shared_ptr<LoadBase> scaled(double factor) const override {
        auto copy = std::make_shared<PlateLoad>(*this);
        for (auto& v : copy->load_vecs)
            v = v * factor;
        return copy;
    }
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

    BeamPolyLoad() : BeamLoadBase(nullptr, BeamLoadAxis::YAxis) {};

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

    // コンストラクタ経由で再構築することで内部のtrapsにも係数を反映する
    std::shared_ptr<LoadBase> scaled(double factor) const override {
        std::vector<double> ws = w;
        for (auto& x : ws) x *= factor;
        if (element == nullptr)
            return std::make_shared<BeamPolyLoad>(ws, params, axis);
        return std::make_shared<BeamPolyLoad>(ws, params, element, axis);
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

    // コンストラクタ経由で再構築することで内部のtrapsにも係数を反映する
    std::shared_ptr<LoadBase> scaled(double factor) const override {
        std::vector<double> ws = w;
        for (auto& x : ws) x *= factor;
        if (element == nullptr)
            return std::make_shared<AxialPolyLoad>(ws, params);
        return std::make_shared<AxialPolyLoad>(ws, params, element);
    }

    BeamStressData GetBeamStress(double p) override;
    Displacement GetDisplacement(double p) override;
};


#endif
#ifndef _COMPONENTS_
#define _COMPONENTS_

#ifndef SWIGCSHARP
#include<iostream>
#include <cmath>
#include <vector>
#include <numeric>
#endif

struct Vector {
public:
    double x, y, z;


    Vector() : x(0), y(0), z(0) {};
    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    double norm() { return sqrt(multiply(*this, *this)); }
    double squared_norm() { return multiply(*this, *this); }

    static Vector XAxis() { return Vector(1, 0, 0); }
    static Vector YAxis() { return Vector(0, 1, 0); }
    static Vector ZAxis() { return Vector(0, 0, 1); }

    static Vector add(const Vector v0, const Vector v1);
    static Vector subtract(const Vector v0, const Vector v1);
    static double multiply(const Vector v0, const Vector v1);
    static Vector multiply(const Vector v0, const double v1);
    static Vector cross(const Vector v0, const Vector v1);
    Vector operator+(const Vector v1) { return add(*this, v1); }
    Vector operator-(const Vector v1) { return subtract(*this, v1); }
    double operator*(const Vector v1) { return multiply(*this, v1); }
    Vector operator*(const double v1) { return multiply(*this, v1); }

    friend std::ostream& operator<<(std::ostream& os, const Vector& m);
};

std::ostream& operator<<(std::ostream& os, const Vector& m);

struct Point
{
public:
    double x, y, z;

    Point() : x(0), y(0), z(0){};
    Point(Vector v) : x(v.x), y(v.y), z(v.z) {};
    Point(double x, double y, double z) : x(x), y(y), z(z){};

    double distance_to(Point p1);
    
    static Vector subtract(const Point p0, const Point p1);
    static Point add(const Point p0, const Point p1);
    static Point divide(const Point p0, const double t);
    Vector operator-(const Point p1) { return subtract(*this, p1); }
    Point operator+(const Point p1) { return add(*this, p1); }
    Point operator/(const double t) { return divide(*this, t); }
    friend std::ostream& operator<<(std::ostream& os, const Point& m);

};

std::ostream& operator<<(std::ostream& os, const Point& m);


struct Plane {
    Point origin;
    Vector ex, ey, ez;
    
    Plane() {}
    Plane(Point o, Vector ex, Vector ey, Vector ez) 
        : origin(o), ex(ex), ey(ey), ez(ez) {}

    Point PointToCoord(Point p);
    double DistanceTo(Point p);
    static Plane CreateFromPoints(Point p0, Point p1, Point p2);
};

struct Displacement {
    Displacement() :Displacement(0, 0, 0, 0, 0, 0) {};
    Displacement(double dx, double dy, double dz, double rx, double ry, double rz);
    double displace[6];
    double Dx() { return displace[0]; }
    double Dy() { return displace[1]; }
    double Dz() { return displace[2]; }
    double Rx() { return displace[3]; }
    double Ry() { return displace[4]; }
    double Rz() { return displace[5]; }

    friend std::ostream& operator<<(std::ostream& os, const Displacement& m);
    Displacement translate(Eigen::Matrix3d transmat);
};

std::ostream& operator<<(std::ostream& os, const Displacement& m);

struct Section
{
public:
    int id = -1;
    double A;
    double Iy;
    double Iz;
    double K;

    Section(double A, double Iy, double Iz, double K)
        : A(A), Iy(Iy), Iz(Iz), K(K){};
};

struct Load {
    Load() : Load(-1, 0, 0, 0, 0, 0, 0) {};
    Load(int _id, double px, double py, double pz);
    Load(int _id, double px, double py, double pz, double mx, double my, double mz);

    int id = -1;
    double loads[6];
    double& Px() { return loads[0]; }
    double& Py() { return loads[1]; }
    double& Pz() { return loads[2]; }
    double& Mx() { return loads[3]; }
    double& My() { return loads[4]; }
    double& Mz() { return loads[5]; }
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
class BeamPolyLoad {
private:
    std::vector<BeamTrapezoidalLoad> traps;
public:
    std::vector<double> w;
    std::vector<double> params;
    double length;

    BeamPolyLoad(const std::vector<double>& w, const std::vector<double>& params, double length)
        : w(w), params(params), length(length) {
        for (int i = 0; i < w.size() - 1; i++) {
            double L1 = params[i] * length;
            double L2 = params[i + 1] * length - L1;
            double L3 = length - L1 - L2;
            traps.push_back(BeamTrapezoidalLoad(w[i], w[i + 1], L1, L2, L3));
        }
    }

    double R0();
    double RA();
    double M0();
    double MA();
    double shear_force(double x);
    double bending_moment(double x);
    double deflection(double x, double EI);
};

struct Support {
    bool flags[6];
    Support() : Support(false, false, false, false, false, false) {};
    Support(bool ux, bool uy, bool uz, bool rx, bool ry, bool rz);

    bool& Ux() { return flags[0]; }
    bool& Uy() { return flags[1]; }
    bool& Uz() { return flags[2]; }
    bool& Rx() { return flags[3]; }
    bool& Ry() { return flags[4]; }
    bool& Rz() { return flags[5]; }

    void FixAll();
    void PinFix();
    void ReleaseAll();
    
};

struct Node {
public:
    int id = -1;
    Point Location;
    Support Fix;

    Node(int id, double x, double y, double z)
        : id(id), Location(x, y, z) {};
    Node(double x, double y, double z) : Location(x, y, z) {};
    Node(Point p) : Location(p) {};

};

struct Material
{
public:
    int id = -1;
    double Young;
    double Poisson;
    double G() { return Young * 0.5 / (1 + Poisson); }

    Material(double young, double poisson) : Young(young), Poisson(poisson) {};

    friend std::ostream& operator<<(std::ostream& os, const Material& m);
};

std::ostream& operator<<(std::ostream& os, const Material& m);

#endif
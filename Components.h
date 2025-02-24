#ifndef _COMPONENTS_
#define _COMPONENTS_

#ifndef SWIGCSHARP
#include<iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
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
    double Dx() const { return displace[0]; }
    double Dy() const { return displace[1]; }
    double Dz() const { return displace[2]; }
    double Rx() const { return displace[3]; }
    double Ry() const { return displace[4]; }
    double Rz() const { return displace[5]; }

    friend std::ostream& operator<<(std::ostream& os, const Displacement& m);
    Displacement translate(Eigen::Matrix3d transmat);
    
    Displacement operator+(const Displacement& other) const {
        return Displacement(
            this->Dx() + other.Dx(),
            this->Dy() + other.Dy(),
            this->Dz() + other.Dz(),
            this->Rx() + other.Rx(),
            this->Ry() + other.Ry(),
            this->Rz() + other.Rz()
        );
    }
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

//struct Load {
//    Load() : Load(-1, 0, 0, 0, 0, 0, 0) {};
//    Load(int _id, double px, double py, double pz);
//    Load(int _id, double px, double py, double pz, double mx, double my, double mz);
//
//    int id = -1;
//    double loads[6];
//    double& Px() { return loads[0]; }
//    double& Py() { return loads[1]; }
//    double& Pz() { return loads[2]; }
//    double& Mx() { return loads[3]; }
//    double& My() { return loads[4]; }
//    double& Mz() { return loads[5]; }
//};

//__declspec(deprecated("Use NodeFix class instead."))
//struct __declspec(deprecated("Use NodeFix class instead.")) Support {

// True if Fixed
//struct Support {
class Support {
public:
    static const bool Fix = true;
    static const bool Free = false;
    static const bool Unlock = false;
    static const bool Lock = true;

    bool lockflags[6];
    bool fixflags[6];
    Support() : Support(Free, Free, Free, Free, Free, Free) {};
    Support(bool ux, bool uy, bool uz, bool rx, bool ry, bool rz);

    bool& Ux() { return fixflags[0]; }
    bool& Uy() { return fixflags[1]; }
    bool& Uz() { return fixflags[2]; }
    bool& Rx() { return fixflags[3]; }
    bool& Ry() { return fixflags[4]; }
    bool& Rz() { return fixflags[5]; }

    std::array<bool, 6> isdof_fixed() {
        std::array<bool, 6> fixed;
        fixed[0] = fixflags[0] || lockflags[0];
        fixed[1] = fixflags[1] || lockflags[1];
        fixed[2] = fixflags[2] || lockflags[2];
        fixed[3] = fixflags[3] || lockflags[3];
        fixed[4] = fixflags[4] || lockflags[4];
        fixed[5] = fixflags[5] || lockflags[5];

        return fixed;
    }

    void FixAll();
    void PinFix();
    void ReleaseAll();

    void UnlockAllRot() {
        lockflags[3] = Unlock;
        lockflags[4] = Unlock;
        lockflags[5] = Unlock;
    }

    bool IsAnyFix() {
        for each (bool f in isdof_fixed()) {
            if (f) return true;
        }
        return false;
        //return std::any_of(fixflags, fixflags + 6, [](bool x) { return x; });
    };
};

//struct NodeFix {
//    int id;
//    bool flags[6];
//    NodeFix() : NodeFix(false, false, false, false, false, false) {};
//    NodeFix(bool ux, bool uy, bool uz, bool rx, bool ry, bool rz);
//    NodeFix(int _id, bool ux, bool uy, bool uz, bool rx, bool ry, bool rz)
//        : NodeFix(ux, uy, uz, rx, ry, rz) {id = _id;};
//
//    bool& Ux() { return flags[0]; }
//    bool& Uy() { return flags[1]; }
//    bool& Uz() { return flags[2]; }
//    bool& Rx() { return flags[3]; }
//    bool& Ry() { return flags[4]; }
//    bool& Rz() { return flags[5]; }
//
//    void FixAll();
//    void PinFix();
//    void ReleaseAll();
//    bool IsAnyFix() {
//        return std::any_of(flags, flags + 6, [](bool x) { return x; });
//    };
//};

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
    double dense = 0;

    Material() : Young(205e3), Poisson(0.2) {};
    Material(double young, double poisson) : Young(young), Poisson(poisson) {};
    Material(double young, double poisson, double dense) : 
        Young(young), Poisson(poisson), dense(dense) {};

    friend std::ostream& operator<<(std::ostream& os, const Material& m);
};

std::ostream& operator<<(std::ostream& os, const Material& m);

#endif
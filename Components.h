#ifndef _COMPONENTS_
#define _COMPONENTS_

#ifndef SWIGCSHARP
#include<iostream>
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
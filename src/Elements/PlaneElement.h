#ifndef _PLANE_ELEMENT_
#define _PLANE_ELEMENT_

#include "ElementBase.h"

Eigen::Vector<double, 4> ShapeFunction4(double xi, double eta);
Eigen::Vector<double, 6> ShapeFunctionTriangle6(double xi, double eta);
Eigen::Vector<double, 8> ShapeFunctionSerendipity8(double xi, double eta);

struct MembraneStressData
{
public:
    double sigx, sigy, sigxy;

    MembraneStressData(double sigx, double sigy, double sigxy) : sigx(sigx), sigy(sigy), sigxy(sigxy) {};

    friend std::ostream &operator<<(std::ostream &os, const MembraneStressData &strs);
};

std::ostream &operator<<(std::ostream &os, const MembraneStressData &strs);

class PlaneElementBase : public ElementBase
{
public:
    double Beta = 0;
    Plane plane;
    Thickness thickness;

    virtual std::vector<NodeLoadData> AreaForceToNodeLoadData(std::vector<Vector> load_vecs) = 0;
};

struct PlateStressData
{
public:
    double Mx = 0, My = 0, Mxy = 0, Qx = 0, Qy = 0;
    double Nx = 0, Ny = 0, Qxy = 0;

    PlateStressData() {};
    PlateStressData(double mx, double my, double mxy,
                    double qx, double qy, double nx, double ny, double qxy) : Mx(mx), My(my), Mxy(mxy), Qx(qx), Qy(qy), Nx(nx), Ny(ny), Qxy(qxy) {};

    // 加算演算子
    friend PlateStressData operator+(const PlateStressData &lhs, const PlateStressData &rhs)
    {
        return PlateStressData(
            lhs.Mx + rhs.Mx,
            lhs.My + rhs.My,
            lhs.Mxy + rhs.Mxy,
            lhs.Qx + rhs.Qx,
            lhs.Qy + rhs.Qy,
            lhs.Nx + rhs.Nx,
            lhs.Ny + rhs.Ny,
            lhs.Qxy + rhs.Qxy);
    }

    // 加算代入演算子
    PlateStressData &operator+=(const PlateStressData &rhs)
    {
        Mx += rhs.Mx;
        My += rhs.My;
        Mxy += rhs.Mxy;
        Qx += rhs.Qx;
        Qy += rhs.Qy;
        Nx += rhs.Nx;
        Ny += rhs.Ny;
        Qxy += rhs.Qxy;
        return *this;
    }

    // 減算演算子
    friend PlateStressData operator-(const PlateStressData &lhs, const PlateStressData &rhs)
    {
        return PlateStressData(
            lhs.Mx - rhs.Mx,
            lhs.My - rhs.My,
            lhs.Mxy - rhs.Mxy,
            lhs.Qx - rhs.Qx,
            lhs.Qy - rhs.Qy,
            lhs.Nx - rhs.Nx,
            lhs.Ny - rhs.Ny,
            lhs.Qxy - rhs.Qxy);
    }

    // 減算代入演算子
    PlateStressData &operator-=(const PlateStressData &rhs)
    {
        Mx -= rhs.Mx;
        My -= rhs.My;
        Mxy -= rhs.Mxy;
        Qx -= rhs.Qx;
        Qy -= rhs.Qy;
        Nx -= rhs.Nx;
        Ny -= rhs.Ny;
        Qxy -= rhs.Qxy;
        return *this;
    }

    // スカラー倍演算子
    friend PlateStressData operator*(const PlateStressData &lhs, double scalar)
    {
        return PlateStressData(
            lhs.Mx * scalar,
            lhs.My * scalar,
            lhs.Mxy * scalar,
            lhs.Qx * scalar,
            lhs.Qy * scalar,
            lhs.Nx * scalar,
            lhs.Ny * scalar,
            lhs.Qxy * scalar);
    }

    friend PlateStressData operator*(double scalar, const PlateStressData &rhs)
    {
        return rhs * scalar;
    }

    // スカラー倍代入演算子
    PlateStressData &operator*=(double scalar)
    {
        Mx *= scalar;
        My *= scalar;
        Mxy *= scalar;
        Qx *= scalar;
        Qy *= scalar;
        Nx *= scalar;
        Ny *= scalar;
        Qxy *= scalar;
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const PlateStressData &strs);
};

std::ostream &operator<<(std::ostream &os, const PlateStressData &strs);

struct PlatePrincipalStressData
{
public:
    double M1 = 0, M2 = 0, theta_m1 = 0;
    double N1 = 0, N2 = 0, theta_n1 = 0;

    PlatePrincipalStressData(double m1, double m2, double theta_m1,
                             double n1, double n2, double theta_n1)
        : M1(m1), M2(m2), theta_m1(theta_m1), N1(n1), N2(n2), theta_n1(theta_n1) {};

    PlatePrincipalStressData(PlateStressData psd);
};



// struct PlateStress {
//     MembraneStressData plane_stress;
//     PlateStressData plate_stress;
//
//     PlateStress(Displacement d0, Displacement d1, Displacement d2);
// };





#endif
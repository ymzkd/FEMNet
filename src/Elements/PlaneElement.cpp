#include "PlaneElement.h"

std::ostream &operator<<(std::ostream &os, const MembraneStressData &strs)
{
    os << "Membrane Stress  " << "sig x: " << strs.sigx << "sig y: " << strs.sigy
       << "sig xy: " << strs.sigxy;
    return os;
}

std::ostream &operator<<(std::ostream &os, const PlateStressData &strs)
{
    os << "Plate Stress Mx: " << strs.Mx << ", My: " << strs.My
       << ", Mxy: " << strs.Mxy << ", Qx: " << strs.Qx << ", Qy: " << strs.Qy;
    return os;
}

Eigen::Vector<double, 4> ShapeFunction4(double xi, double eta)
{
    Eigen::Vector<double, 4> Nfunc;
    Nfunc << 0.25 * (1 - xi) * (1 - eta),
        0.25 * (1 + xi) * (1 - eta),
        0.25 * (1 + xi) * (1 + eta),
        0.25 * (1 - xi) * (1 + eta);
    return Nfunc;
}

Eigen::Vector<double, 6> ShapeFunctionTriangle6(double xi, double eta)
{
    double gamma = 1.0 - xi - eta;
    Eigen::Vector<double, 6> Nfunc;
    Nfunc << 2.0 * xi * xi - xi,
        2.0 * eta * eta - eta,
        2.0 * gamma * gamma - gamma,
        4.0 * eta * gamma,
        4.0 * gamma * xi,
        4.0 * xi * eta;
    return Nfunc;
}

Eigen::Vector<double, 8> ShapeFunctionSerendipity8(double xi, double eta)
{
    Eigen::Vector<double, 8> Nfunc;
    Nfunc << 0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1),
        0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1),
        0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1),
        0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1),
        0.5 * (1 - xi * xi) * (1 - eta),
        0.5 * (1 + xi) * (1 - eta * eta),
        0.5 * (1 - xi * xi) * (1 + eta),
        0.5 * (1 - xi) * (1 - eta * eta);
    return Nfunc;
}

PlatePrincipalStressData::PlatePrincipalStressData(PlateStressData psd)
{
    double dm = sqrt(pow(((psd.Mx - psd.My) * 0.5), 2) + psd.Mxy * psd.Mxy);
    double m1 = (psd.Mx + psd.My) * 0.5 + dm;
    double m2 = (psd.Mx + psd.My) * 0.5 - dm;

    double theta_m = 0;
    if (abs(psd.Mx - psd.My) > 0.00001)
    {
        double t2m = 2 * psd.Mxy / (psd.Mx - psd.My);
        theta_m = atan2(2 * psd.Mxy, (psd.Mx - psd.My)) / 2;
        // theta_m = atan(t2m) / 2;
    }

    double dn = sqrt(pow(((psd.Nx - psd.Ny) * 0.5), 2) + psd.Qxy * psd.Qxy);
    double n1 = (psd.Nx + psd.Ny) * 0.5 + dn;
    double n2 = (psd.Nx + psd.Ny) * 0.5 - dn;

    double theta_n = 0;
    if (abs(psd.Nx - psd.Ny) > 0.00001)
    {
        double t2n = 2 * psd.Qxy / (psd.Nx - psd.Ny);
        theta_n = atan2(2 * psd.Qxy, (psd.Nx - psd.Ny)) / 2;
        // theta_n = atan(t2n) / 2;
    }

    M1 = m1;
    M2 = m2;
    theta_m1 = theta_m;

    N1 = n1;
    N2 = n2;
    theta_n1 = theta_n;
}


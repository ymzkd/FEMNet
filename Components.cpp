#include <Eigen/Dense>

#include "Components.h"

//std::ostream& Material::operator<<(std::ostream& os) {
//    os << "id: " << id << std::endl;
//    os << "Young: " << Young << std::endl;
//    os << "Poisson: " << Poisson << std::endl;
//    return os;
//}

std::ostream& operator<<(std::ostream& os, const Vector& m)
{
    os << m.x << ", " << m.y << ", " << m.z;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Point& m)
{
    os << "(" << m.x << ", " << m.y << ", " << m.z << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Displacement& m)
{
    os << "dx: " << m.displace[0] << ", dy: " << m.displace[1]
        << ", dz: " << m.displace[2] << ", rx: " << m.displace[3]
        << ", ry: " << m.displace[4] << ", rz: " << m.displace[5] << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Material& m) {
    os << "id: " << m.id << std::endl;
    os << "Young: " << m.Young << std::endl;
    os << "Poisson: " << m.Poisson << std::endl;
    return os;
}

//double Vector::norm()
//{
//    return sqrt(multiply(*this, *this));
//}

Vector Vector::add(const Vector v0, const Vector v1)
{
    return Vector(v0.x+v1.x, v0.y+v1.y, v0.z + v1.z);
}

Vector Vector::subtract(const Vector v0, const Vector v1)
{
    return Vector(v0.x - v1.x, v0.y - v1.y, v0.z - v1.z);
}

double Vector::multiply(const Vector v0, const Vector v1)
{
    return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

Vector Vector::multiply(const Vector v0, const double v1)
{
    return Vector(v0.x * v1, v0.y * v1, v0.z * v1);
}

Vector Vector::cross(const Vector v0, const Vector v1)
{
    double x = v0.y * v1.z - v0.z * v1.y;
    double y = v0.z * v1.x - v0.x * v1.z;
    double z = v0.x * v1.y - v0.y * v1.x;
    return Vector(x, y, z);
}

double Point::distance_to(Point p1)
{
    return subtract(*this, p1).norm();
}

Vector Point::subtract(const Point p0, const Point p1)
{
    return Vector(p0.x - p1.x, p0.y - p1.y, p0.z - p1.z);
}

Point Point::add(const Point p0, const Point p1)
{
    return Point(p0.x + p1.x, p0.y + p1.y, p0.z + p1.z);
}

Point Point::divide(const Point p0, const double t)
{
    return Point(p0.x / t, p0.y / t, p0.z / t);
}

Support::Support(bool ux, bool uy, bool uz, bool rx, bool ry, bool rz)
{
    flags[0] = ux;
    flags[1] = uy;
    flags[2] = uz;
    flags[3] = rx;
    flags[4] = ry;
    flags[5] = rz;
}

// Constrain translational movement and rotation.
void Support::FixAll()
{
    Ux() = true;
    Uy() = true;
    Uz() = true;
    Rx() = true;
    Ry() = true;
    Rz() = true;
}

void Support::PinFix()
{
    Ux() = true;
    Uy() = true;
    Uz() = true;
    Rx() = false;
    Ry() = false;
    Rz() = false;
}

void Support::ReleaseAll()
{
    Ux() = false;
    Uy() = false;
    Uz() = false;
    Rx() = false;
    Ry() = false;
    Rz() = false;
}

Displacement::Displacement(double dx, double dy, double dz, double rx, double ry, double rz)
{
    displace[0] = dx;
    displace[1] = dy;
    displace[2] = dz;
    displace[3] = rx;
    displace[4] = ry;
    displace[5] = rz;
}

Displacement Displacement::translate(Eigen::Matrix3d transmat)
{
    Eigen::Vector3d dvec(displace[0], displace[1], displace[2]);
    Eigen::Vector3d rvec(displace[3], displace[4], displace[5]);
    
    dvec = transmat * dvec;
    rvec = transmat * rvec;
    return Displacement(dvec(0), dvec(1), dvec(2), rvec(0), rvec(1), rvec(2));
}

Load::Load(int _id, double px, double py, double pz)
{
    id = _id;
    loads[0] = px;
    loads[1] = py;
    loads[2] = pz;
    loads[3] = 0;
    loads[4] = 0;
    loads[5] = 0;
}

Load::Load(int _id, double px, double py, double pz, double mx, double my, double mz)
{
    id = _id;
    loads[0] = px;
    loads[1] = py;
    loads[2] = pz;
    loads[3] = mx;
    loads[4] = my;
    loads[5] = mz;
}

Point Plane::PointToCoord(Point p)
{
    Vector v0(p - origin);
    double x = v0 * ex;
    double y = v0 * ey;
    double z = v0 * ez;
    return Point(x, y, z);
}

double Plane::DistanceTo(Point p)
{
    return (p - origin)* ez;
}

Plane Plane::CreateFromPoints(Point p0, Point p1, Point p2)
{
    Vector v01 = p1 - p0;
    Vector v02 = p2 - p0;
    Vector vn = Vector::cross(v01, v02);

    Vector ex = v01 * (1.0 / v01.norm());
    
    // Vector dv(v02.x * vn.x, v02.y * vn.y, v02.z * vn.z);
    // v02 = v02 - dv;
    Vector ey = v02 - ex * (v02 * ex);
    ey = ey * (1.0 / ey.norm());
    Vector ez = vn * (1.0 / (vn.norm()));

    return Plane(p0, ex, ey, ez);
}

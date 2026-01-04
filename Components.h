#ifndef _COMPONENTS_
#define _COMPONENTS_

#ifndef SWIG
#include<iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

#include <Eigen/Dense>
#endif

//#define PI 3.141592653589793238462643
#define NODE_DOF 6

/// <summary>
/// 6自由度（Ux, Uy, Uz, Rx, Ry, Rz）のフラグを管理する共通データ型
/// </summary>
struct DOFFlags {
    bool flags[6];

    // デフォルトコンストラクタ（全false）
    DOFFlags() {
        for (int i = 0; i < 6; i++) flags[i] = false;
    }

    // 6個のboolを受け取るコンストラクタ
    DOFFlags(bool ux, bool uy, bool uz, bool rx, bool ry, bool rz) {
        flags[0] = ux; flags[1] = uy; flags[2] = uz;
        flags[3] = rx; flags[4] = ry; flags[5] = rz;
    }

    // インデックスベースアクセサ（範囲チェック付き）
    bool Get(int index) const {
        if (index < 0 || index >= 6)
            throw std::out_of_range("DOFFlags index out of range [0-5]");
        return flags[index];
    }

    void Set(int index, bool value) {
        if (index < 0 || index >= 6)
            throw std::out_of_range("DOFFlags index out of range [0-5]");
        flags[index] = value;
    }

    // DOF名ベースアクセサ（Supportクラスパターンと一貫性）
    bool& Ux() { return flags[0]; }
    bool& Uy() { return flags[1]; }
    bool& Uz() { return flags[2]; }
    bool& Rx() { return flags[3]; }
    bool& Ry() { return flags[4]; }
    bool& Rz() { return flags[5]; }

    const bool& Ux() const { return flags[0]; }
    const bool& Uy() const { return flags[1]; }
    const bool& Uz() const { return flags[2]; }
    const bool& Rx() const { return flags[3]; }
    const bool& Ry() const { return flags[4]; }
    const bool& Rz() const { return flags[5]; }

    // ヘルパーメソッド
    void SetAll(bool value) {
        for (int i = 0; i < 6; i++) flags[i] = value;
    }

    bool IsAny() const {
        for (int i = 0; i < 6; i++)
            if (flags[i]) return true;
        return false;
    }

    bool IsAll() const {
        for (int i = 0; i < 6; i++)
            if (!flags[i]) return false;
        return true;
    }
};

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
    friend Vector operator*(const double v1, const Vector v0) { return multiply(v0, v1); }

	Eigen::Vector3d toEigen() const {
		return Eigen::Vector3d(x, y, z);
	}

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
    Plane(const Plane& other)
        : origin(other.origin), ex(other.ex), ey(other.ey), ez(other.ez) {}

    Point PointToCoord(Point p);
    double DistanceTo(Point p);

    /// @brief 平面を指定された軸周りに回転させる
    /// @param angle 回転角度（ラジアン）
    /// @param axis 回転軸ベクトル
    void Rotate(double angle, Vector axis);

    /**
     * @brief Creates a plane from three given points.
     *
     * This function constructs a plane by defining its coordinate system based on the three input points.
     * The x-axis is determined by the p0 -> p1 direction, the y-axis is the vector from the p0->p2 direction 
     * excluding the component parallel to the x-axis, and the z-axis is the direction orthogonal to the two axes. 
     * The origin is assumed to be p0.
     *
     * @param p0 The first point on the plane, which will be used as the origin of the plane's coordinate system.
     * @param p1 The second point on the plane, used to define the x-axis direction.
     * @param p2 The third point on the plane, used to define the y-axis direction.
     *
     * @return A Plane object representing the plane defined by the three points, including its origin and orthonormal basis vectors.
     *
     * @note The function assumes that the three points are not collinear. If they are, the resulting plane may be undefined.
     */
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

    // 減算演算子
    Displacement operator-(const Displacement& other) const {
        return Displacement(
            this->Dx() - other.Dx(),
            this->Dy() - other.Dy(),
            this->Dz() - other.Dz(),
            this->Rx() - other.Rx(),
            this->Ry() - other.Ry(),
            this->Rz() - other.Rz()
        );
    }

    // スカラー乗算
    Displacement operator*(double scalar) const {
        return Displacement(
            this->Dx() * scalar,
            this->Dy() * scalar,
            this->Dz() * scalar,
            this->Rx() * scalar,
            this->Ry() * scalar,
            this->Rz() * scalar
        );
    }

    // スカラー除算
    Displacement operator/(double scalar) const {
        return Displacement(
            this->Dx() / scalar,
            this->Dy() / scalar,
            this->Dz() / scalar,
            this->Rx() / scalar,
            this->Ry() / scalar,
            this->Rz() / scalar
        );
    }

    // 単項マイナス（符号反転）
    Displacement operator-() const {
        return Displacement(
            -this->Dx(),
            -this->Dy(),
            -this->Dz(),
            -this->Rx(),
            -this->Ry(),
            -this->Rz()
        );
    }

    // 複合代入演算子
    Displacement& operator+=(const Displacement& other) {
        displace[0] += other.Dx();
        displace[1] += other.Dy();
        displace[2] += other.Dz();
        displace[3] += other.Rx();
        displace[4] += other.Ry();
        displace[5] += other.Rz();
        return *this;
    }

    Displacement& operator-=(const Displacement& other) {
        displace[0] -= other.Dx();
        displace[1] -= other.Dy();
        displace[2] -= other.Dz();
        displace[3] -= other.Rx();
        displace[4] -= other.Ry();
        displace[5] -= other.Rz();
        return *this;
    }

    Displacement& operator*=(double scalar) {
        for (int i = 0; i < 6; ++i) {
            displace[i] *= scalar;
        }
        return *this;
    }

    Displacement& operator/=(double scalar) {
        for (int i = 0; i < 6; ++i) {
            displace[i] /= scalar;
        }
        return *this;
    }

    // インデックスアクセス演算子
    double& operator[](int index) {
        return displace[index];
    }

    const double& operator[](int index) const {
        return displace[index];
    }

    // スカラー × Displacement の順序での乗算（friend関数）
    friend Displacement operator*(double scalar, const Displacement& disp) {
        return disp * scalar;
    }

	Vector Translation() const {
		return Vector(displace[0], displace[1], displace[2]);
	}

	Vector Rotation() const {
		return Vector(displace[3], displace[4], displace[5]);
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
    double Iyz = 0;
    double K;

    Section() : A(0), Iy(0), Iz(0), K(0) {};
    Section(double A, double Iy, double Iz, double K)
        : A(A), Iy(Iy), Iz(Iz), K(K){};
};

// True if Fixed
class Support : public DOFFlags {
public:
    static const bool Fix = true;
    static const bool Free = false;
    static const bool Unlock = false;
    static const bool Lock = true;

    DOFFlags lockflags;
    Support() : DOFFlags(Free, Free, Free, Free, Free, Free), lockflags() {}
    Support(bool ux, bool uy, bool uz, bool rx, bool ry, bool rz);

    // Ux(), Uy() 等は DOFFlags から継承されるため定義不要

    std::array<bool, 6> isdof_fixed() {
        std::array<bool, 6> fixed;
        for (int i = 0; i < 6; i++)
            fixed[i] = flags[i] || lockflags.flags[i];
        return fixed;
    }

    void FixAll();
    void PinFix();
    void ReleaseAll();

    void UnlockAllRot() {
        lockflags.Rx() = Unlock;
        lockflags.Ry() = Unlock;
        lockflags.Rz() = Unlock;
    }

    bool IsAnyFix() {
        auto fixed = isdof_fixed();
        for (bool f : fixed) {
            if (f) return true;
        }
        return false;
    }
};

struct NodeMass {
public:
    double Mass = 0;
    double LoadMass = 0;
    double ElementMass = 0;

    NodeMass() {};
    NodeMass(double node_mass, double load_mass, double element_mass)
        : Mass(node_mass), LoadMass(load_mass), ElementMass(element_mass) {
    };

    double SumMass() const {
        return Mass + LoadMass + ElementMass;
    }
};

struct Node {
public:
    int id = -1;
    Point Location;
    Support Fix;
	//double Mass = 0;
	NodeMass MassData;

    Node() : id(-1), Location(0, 0, 0) {};
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

struct Thickness {
public:
    double plane_thick = 0;
    double plate_thick = 0;
    double weight_thick = 0;

    
    Thickness(double thick)
     : plane_thick(thick), plate_thick(thick), weight_thick(thick) {};
    
    Thickness(double plane, double plate)
     : plane_thick(plane), plate_thick(plate), weight_thick(plane) {};

    Thickness(double plane, double plate, double weight)
     : plane_thick(plane), plate_thick(plate), weight_thick(weight) {};

    Thickness() : plane_thick(0), plate_thick(0), weight_thick(0) {};
};

struct NodeLoadData {
public:
    NodeLoadData() : NodeLoadData(-1, 0, 0, 0, 0, 0, 0) {};
    NodeLoadData(int _id, double px, double py, double pz);
    NodeLoadData(int _id, double px, double py, double pz, double mx, double my, double mz);

    int id = -1;
    double loads[6];
    double& Px() { return loads[0]; }
    double& Py() { return loads[1]; }
    double& Pz() { return loads[2]; }
    double& Mx() { return loads[3]; }
    double& My() { return loads[4]; }
    double& Mz() { return loads[5]; }

    double Px() const { return loads[0]; }
    double Py() const { return loads[1]; }
    double Pz() const { return loads[2]; }
    double Mx() const { return loads[3]; }
    double My() const { return loads[4]; }
    double Mz() const { return loads[5]; }

    // 二項演算子
    NodeLoadData operator+(const NodeLoadData& other) const;
    NodeLoadData operator-(const NodeLoadData& other) const;
    NodeLoadData operator*(double scalar) const;
    NodeLoadData operator/(double scalar) const;
    
    // 単項演算子
    NodeLoadData operator-() const;
    
    // 複合代入演算子
    NodeLoadData& operator+=(const NodeLoadData& other);
    NodeLoadData& operator-=(const NodeLoadData& other);
    NodeLoadData& operator*=(double scalar);
    NodeLoadData& operator/=(double scalar);
    
    // friend関数
    friend NodeLoadData operator*(double scalar, const NodeLoadData& load);
    friend std::ostream& operator<<(std::ostream& os, const NodeLoadData& load);

private:
    // IDを決定するヘルパー関数
    static int determineId(int id1, int id2);
};

std::ostream& operator<<(std::ostream& os, const NodeLoadData& load);

#endif
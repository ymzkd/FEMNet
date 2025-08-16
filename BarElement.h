#ifndef _BAR_ELEMENT_
#define _BAR_ELEMENT_

#ifndef SWIGCSHARP

// #ifdef USE_MKL
// #define EIGEN_USE_MKL_ALL
// #endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif

#include "Element.h"

struct BeamStressData
{
public:
    double Nx, Qy, Qz, Mx, My, Mz;

    BeamStressData() : Nx(0), Qy(0), Qz(0), Mx(0), My(0), Mz(0) {};
    BeamStressData(double nx, double qy, double qz, double mx, double my, double mz) : Nx(nx), Qy(qy), Qz(qz), Mx(mx), My(my), Mz(mz) {};

    // 加算演算子
    friend BeamStressData operator+(const BeamStressData &lhs, const BeamStressData &rhs)
    {
        return BeamStressData(
            lhs.Nx + rhs.Nx,
            lhs.Qy + rhs.Qy,
            lhs.Qz + rhs.Qz,
            lhs.Mx + rhs.Mx,
            lhs.My + rhs.My,
            lhs.Mz + rhs.Mz);
    }

    // 加算代入演算子
    BeamStressData &operator+=(const BeamStressData &rhs)
    {
        Nx += rhs.Nx;
        Qy += rhs.Qy;
        Qz += rhs.Qz;
        Mx += rhs.Mx;
        My += rhs.My;
        Mz += rhs.Mz;
        return *this;
    }

    // 減算演算子
    friend BeamStressData operator-(const BeamStressData &lhs, const BeamStressData &rhs)
    {
        return BeamStressData(
            lhs.Nx - rhs.Nx,
            lhs.Qy - rhs.Qy,
            lhs.Qz - rhs.Qz,
            lhs.Mx - rhs.Mx,
            lhs.My - rhs.My,
            lhs.Mz - rhs.Mz);
    }

    // 減算代入演算子
    BeamStressData &operator-=(const BeamStressData &rhs)
    {
        Nx -= rhs.Nx;
        Qy -= rhs.Qy;
        Qz -= rhs.Qz;
        Mx -= rhs.Mx;
        My -= rhs.My;
        Mz -= rhs.Mz;
        return *this;
    }

    // スカラー倍演算子
    friend BeamStressData operator*(const BeamStressData &lhs, double scalar)
    {
        return BeamStressData(
            lhs.Nx * scalar,
            lhs.Qy * scalar,
            lhs.Qz * scalar,
            lhs.Mx * scalar,
            lhs.My * scalar,
            lhs.Mz * scalar);
    }

    friend BeamStressData operator*(double scalar, const BeamStressData &rhs)
    {
        return rhs * scalar;
    }

    // スカラー倍代入演算子
    BeamStressData &operator*=(double scalar)
    {
        Nx *= scalar;
        Qy *= scalar;
        Qz *= scalar;
        Mx *= scalar;
        My *= scalar;
        Mz *= scalar;
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const BeamStressData &bsd);
};

std::ostream &operator<<(std::ostream &os, const BeamStressData &bsd);

// 端部応力を格納。分布荷重とかも考える。
struct BeamStress
{
public:
    BeamStressData S0, S1;

    BeamStress(BeamStressData s0, BeamStressData s1) : S0(s0), S1(s1) {};

    /// <summary>
    /// 両端で定義された梁要素応力を補間して中間応力計算
    /// </summary>
    /// <param name="t">0~1の位置パラメータ</param>
    /// <returns></returns>
    BeamStressData Interpolate(double t)
    {
        return BeamStressData(
            S0.Nx + (S1.Nx - S0.Nx) * t,
            S0.Qy + (S1.Qy - S0.Qy) * t,
            S0.Qz + (S1.Qz - S0.Qz) * t,
            S0.Mx + (S1.Mx - S0.Mx) * t,
            S0.My + (S1.My - S0.My) * t,
            S0.Mz + (S1.Mz - S0.Mz) * t);
    };

    friend std::ostream &operator<<(std::ostream &os, const BeamStress &bsd);
};

std::ostream &operator<<(std::ostream &os, const BeamStress &bsd);

class BarElementBase : public ElementBase
{
protected:
    static const int node_num = 2;

public:
    Node *Nodes[2];

    Section *Sec;

    std::vector<Node *> NodesList() override { return std::vector<Node *>{Nodes[0], Nodes[1]}; }
    Node *ni() { return Nodes[0]; }
    Node *nj() { return Nodes[1]; }
    double length();
    Eigen::VectorXd NodeLumpedMass()
    {
        return Eigen::VectorXd::Constant(node_num, Sec->A * length() * Mat.dense / node_num);
    }

    void AssembleMassMatrix(Eigen::SparseMatrix<double> &mat);
    virtual BeamStress stress(Displacement d0, Displacement d1) = 0;
};

class TrussElement : public BarElementBase
{
private:
    // static constexpr ElementType type = ElementType::Truss;
    static const int total_dof = 6;
    static const int node_num = 2;
    Eigen::Matrix<double, node_num, total_dof> trans_matrix();
    Eigen::MatrixXd stiffness_matrix_local();

    // トラス要素の幾何剛性行列(6x6)
    Eigen::MatrixXd geometric_local_stiffness_matrix(const std::vector<Displacement> &disp) override;
    // double element_length();
public:
    // Node* Nodes[2];
    // Section* Sec;

    TrussElement() {}
    TrussElement(Node *n0, Node *n1, Section *sec, Material mat);
    TrussElement(int _id, Node *n0, Node *n1, Section *sec, Material mat)
        : TrussElement(n0, n1, sec, mat)
    {
        id = _id;
    };

    // トラス要素の幾何剛性行列(6x6)
    Eigen::MatrixXd StiffnessMatrix();

    bool hasRotate() { return false; }

    int NodeNum() override { return node_num; }

    ElementType Type() { return ElementType::Truss; }
    // ElementType Type() { return type; }

    int TotalDof() { return total_dof; }

    void AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K);
    void AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat) override;
    void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp) override;

    Eigen::MatrixXd NodeConsistentMass() override;

    std::vector<NodeLoadData> InertialForceToNodeLoadData(Eigen::Vector3d accel_vec) override;

    BeamStress stress(Displacement d0, Displacement d1);
};

class BeamElement : public BarElementBase
{
protected:
    static const int total_dof = 12;
    double element_length();
    Eigen::MatrixXd trans_matrix();
    virtual Eigen::MatrixXd stiffness_matrix_local();

    Eigen::MatrixXd geometric_local_stiffness_matrix(const std::vector<Displacement> &disp) override;

private:
    // static constexpr ElementType type = ElementType::Beam;

public:
    double Beta;
    // Node* Nodes[2];
    // Section* Sec;
    Vector XAxis()
    {
        Eigen::Matrix3d tr0 = trans_matrix3(Nodes[0]->Location, Nodes[1]->Location, Beta);
        Eigen::Vector3d vec = tr0.row(0);
        return Vector(vec.x(), vec.y(), vec.z());
    }

    Vector YAxis()
    {
        Eigen::Matrix3d tr0 = trans_matrix3(Nodes[0]->Location, Nodes[1]->Location, Beta);
        Eigen::Vector3d vec = tr0.row(1);
        return Vector(vec.x(), vec.y(), vec.z());
    }

    Vector ZAxis()
    {
        Eigen::Matrix3d tr0 = trans_matrix3(Nodes[0]->Location, Nodes[1]->Location, Beta);
        Eigen::Vector3d vec = tr0.row(2);
        return Vector(vec.x(), vec.y(), vec.z());
    }

    BeamElement() {}
    // BeamElement() : ElementBase(ElementType::Beam) {}
    BeamElement(Node *n0, Node *n1, Section *sec, Material mat, double beta = 0);
    BeamElement(int _id, Node *n0, Node *n1, Section *sec, Material mat, double beta = 0)
        : BeamElement(n0, n1, sec, mat, beta)
    {
        id = _id;
    };

    virtual Eigen::MatrixXd StiffnessMatrix();
    bool hasRotate() { return true; }
    int NodeNum() override { return node_num; }
    // ElementType Type() { return type; }
    ElementType Type() { return ElementType::Beam; }
    int TotalDof() { return total_dof; }

    void AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K);
    void AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat) override;
    void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp) override;

    std::vector<NodeLoadData> InertialForceToNodeLoadData(Eigen::Vector3d accel_vec) override;
    Eigen::MatrixXd NodeConsistentMass() override;

    // double length();
    BeamStress stress(Displacement d0, Displacement d1);
    Displacement DisplaceAt(Displacement d0, Displacement d1, double p);
};

class ComplexBeamElement : public BeamElement
{
private:
    // static const int total_dof = 12;
    Eigen::MatrixXd stiffness_matrix_local() override;
    // double element_length();
    // Eigen::MatrixXd trans_matrix();
public:
    double Lambda_bz, Lambda_bz_, Lambda_sy, Lambda_sy_;
    double Lambda_by, Lambda_by_, Lambda_sz, Lambda_sz_;
    double lzi, lzj, lyi, lyj;

    ComplexBeamElement();
    ComplexBeamElement(Node *n0, Node *n1, Section *sec, Material mat, double beta = 0);
    ComplexBeamElement(int _id, Node *n0, Node *n1, Section *sec, Material mat, double beta = 0);

    Eigen::MatrixXd StiffnessMatrix() override;

    // Eigen::MatrixXd NodeConsistentMass();
};

#endif
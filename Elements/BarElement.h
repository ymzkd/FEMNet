#ifndef _BAR_ELEMENT_
#define _BAR_ELEMENT_

#include "ElementBase.h"

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



#endif
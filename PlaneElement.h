#ifndef _PLANE_ELEMENT_
#define _PLANE_ELEMENT_

#include "Element.h"


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

class TriPlaneElement : public PlaneElementBase
{
    friend class TriPlateElement;

private:
    static constexpr int node_num = 3;
    static constexpr int node_dof = 3;
    static constexpr int total_dof = 9;

    static constexpr ElementType type = ElementType::Membrane;
    Eigen::MatrixXd BMatrix();
    Eigen::Matrix3d DMatrix();
    Eigen::MatrixXd localStiffnessMatrix();

    Eigen::MatrixXd geometric_local_stiffness_matrix(const std::vector<Displacement> &disp) override;

    Eigen::MatrixXd trans_matrix();

public:
    Node *Nodes[3];
    std::vector<Node *> NodesList() override { return std::vector<Node *>{Nodes[0], Nodes[1], Nodes[2]}; }

    TriPlaneElement() {};
    TriPlaneElement(Node *n0, Node *n1, Node *n2, double t, Material mat, double beta = 0);
    TriPlaneElement(Node *n0, Node *n1, Node *n2, Thickness t, Material mat, double beta = 0);
    TriPlaneElement(int _id, Node *n0, Node *n1, Node *n2, double t, Material mat, double beta = 0)
        : TriPlaneElement(n0, n1, n2, t, mat, beta)
    {
        id = _id;
    };
    TriPlaneElement(int _id, Node *n0, Node *n1, Node *n2, Thickness t, Material mat, double beta = 0)
        : TriPlaneElement(n0, n1, n2, t, mat, beta)
    {
        id = _id;
    };

    /**
     * Calculates the node lumped mass for an element.
     *
     * @return A vector containing the lumped mass values for each node.
     */
    Eigen::VectorXd NodeLumpedMass()
    {
        return Eigen::VectorXd::Constant(node_num, Area() * thickness.weight_thick * Mat.dense / node_num);
    }

    std::vector<NodeLoadData> AreaForceToNodeLoadData(std::vector<Vector> load_vecs) override;

    std::vector<NodeLoadData> InertialForceToNodeLoadData(Eigen::Vector3d accel_vec) override;
    Eigen::MatrixXd NodeConsistentMass() override;

    double Area();
    ElementType Type() { return type; }
    int NodeNum() override { return node_num; }
    int TotalDof() { return total_dof; }
    Eigen::MatrixXd StiffnessMatrix();

    void AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K);
    void AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat) override;
    void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp) override;
    void AssembleMassMatrix(Eigen::SparseMatrix<double> &mat);

    MembraneStressData stress(Displacement d0, Displacement d1, Displacement d2);
};

class QuadPlaneElement : public PlaneElementBase
{
    friend class QuadPlateElement;

private:
    static constexpr int node_num = 4;
    static constexpr int node_dof_local = 2;
    static constexpr int node_dof = 3;
    static constexpr int total_dof = node_num * node_dof;
    static constexpr int total_dof_local = node_num * node_dof_local;

    static constexpr ElementType type = ElementType::Membrane;
    Eigen::Matrix2d JMatrix(double xi, double eta);
    Eigen::MatrixXd BMatrix(double xi, double eta);
    Eigen::Matrix3d DMatrix();
    Eigen::MatrixXd localStiffnessMatrix();
    Eigen::MatrixXd geometric_local_stiffness_matrix(const std::vector<Displacement> &disp) override;

    Eigen::MatrixXd trans_matrix();

public:
    Node *Nodes[4];
    std::vector<Node *> NodesList() override { return std::vector<Node *>{Nodes[0], Nodes[1], Nodes[2], Nodes[3]}; }
    // Plane plane;
    // Thickness thickness;

    QuadPlaneElement() {};
    QuadPlaneElement(Node *n0, Node *n1, Node *n2, Node *n3, double t, Material mat);
    QuadPlaneElement(Node *n0, Node *n1, Node *n2, Node *n3, Thickness t, Material mat);
    QuadPlaneElement(int _id, Node *n0, Node *n1, Node *n2, Node *n3, double t, Material mat) : QuadPlaneElement(n0, n1, n2, n3, t, mat)
    {
        id = _id;
    };
    QuadPlaneElement(int _id, Node *n0, Node *n1, Node *n2, Node *n3, Thickness t, Material mat) : QuadPlaneElement(n0, n1, n2, n3, t, mat)
    {
        id = _id;
    };

    Eigen::VectorXd NodeLumpedMass()
    {
        return Eigen::VectorXd::Constant(node_num, Area() * thickness.weight_thick * Mat.dense / node_num);
    }

    Eigen::MatrixXd NodeConsistentMass();

    std::vector<NodeLoadData> InertialForceToNodeLoadData(Eigen::Vector3d accel_vec);
    std::vector<NodeLoadData> AreaForceToNodeLoadData(std::vector<Vector> load_vecs) override;

    double Area();
    ElementType Type() { return type; }
    int NodeNum() override { return node_num; }
    int TotalDof() { return total_dof; }
    Eigen::MatrixXd StiffnessMatrix();

    void AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K);
    void AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat) override;
    void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp) override;
    void AssembleMassMatrix(Eigen::SparseMatrix<double> &mat);

    MembraneStressData stress(Displacement d0, Displacement d1,
                              Displacement d2, Displacement d3, double xi, double eta);
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

class TriPlateElement : public PlaneElementBase
{
private:
    static constexpr int node_num = 3;
    static constexpr int total_dof = 18;
    static constexpr int node_dof = 6;
    static constexpr ElementType type = ElementType::DKT;
    using LocalMatrixd = Eigen::Matrix<double, total_dof, total_dof>;

    Eigen::MatrixXd HVecs(Eigen::Vector<double, 6> shape_funcs);

    Eigen::MatrixXd BMatrix(double L2, double L3);

    Eigen::Matrix3d DMatrix();

    Eigen::MatrixXd localStiffnessMatrix();

    Eigen::MatrixXd geometric_local_stiffness_matrix(const std::vector<Displacement> &disp) override;

    LocalMatrixd trans_matrix();

public:
    Node *Nodes[3];
    std::vector<Node *> NodesList() override { return std::vector<Node *>{Nodes[0], Nodes[1], Nodes[2]}; }
    // Plane plane;
    // Thickness thickness;

    TriPlaneElement plane_element;

    TriPlateElement() {};
    TriPlateElement(Node *n0, Node *n1, Node *n2, Thickness t, Material mat);
    TriPlateElement(Node *n0, Node *n1, Node *n2, double t, Material mat);
    TriPlateElement(int _id, Node *n0, Node *n1, Node *n2, Thickness t, Material mat)
        : TriPlateElement(n0, n1, n2, t, mat)
    {
        id = _id;
    };
    TriPlateElement(int _id, Node *n0, Node *n1, Node *n2, double t, Material mat)
        : TriPlateElement(n0, n1, n2, t, mat)
    {
        id = _id;
    };

    Eigen::VectorXd NodeLumpedMass()
    {
        return Eigen::VectorXd::Constant(node_num, Area() * thickness.weight_thick * Mat.dense / node_num);
        // return Eigen::VectorXd::Constant(node_num, Area() * thickness * Mat.dense / node_num);
    }

    std::vector<NodeLoadData> InertialForceToNodeLoadData(Eigen::Vector3d accel_vec) override;
    std::vector<NodeLoadData> AreaForceToNodeLoadData(std::vector<Vector> load_vecs) override;

    Eigen::MatrixXd NodeConsistentMass() override;

    double Area();
    ElementType Type() { return type; }
    int NodeNum() override { return node_num; }
    int TotalDof() { return total_dof; }
    Eigen::MatrixXd StiffnessMatrix();
    Eigen::MatrixXd GeometricStiffnessMatrix(const std::vector<Displacement> &disp);

    void AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K);
    void AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat) override;
    void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp) override;
    void AssembleMassMatrix(Eigen::SparseMatrix<double> &mat);

    PlateStressData stress(
        Displacement d0, Displacement d1, Displacement d2, double xi, double eta);

    //  PlateStressData stress_save(
    // Displacement d0, Displacement d1, Displacement d2, double xi, double eta);
    // void shearstress(Displacement d0, Displacement d1, Displacement d2);
};

class QuadPlateElement : public PlaneElementBase
{
public:
    // 下のpublicの範囲でこの定数を使っているからこちらもpublicにしておく必要があるらしい。
    static constexpr int node_num = 4;

private:
    static constexpr int node_dof = 6;
    static constexpr int node_dof_local = 3;
    static constexpr int total_dof = node_dof * node_num;
    static constexpr int total_dof_local = node_dof_local * node_num;

    using LocalMatrixd = Eigen::Matrix<double, total_dof, total_dof>;

    static constexpr ElementType type = ElementType::DKQ;

    Eigen::Matrix2d JMatrix(double xi, double eta);

    Eigen::Matrix2d dJinv_dxi(double xi, double eta);
    Eigen::Matrix2d dJinv_deta(double xi, double eta);

    Eigen::MatrixXd HVecs(Eigen::VectorXd shape_funcs);
    Eigen::MatrixXd BMatrix(double xi, double eta);

    Eigen::Matrix3d DMatrix();
    Eigen::MatrixXd localStiffnessMatrix();

    Eigen::MatrixXd geometric_local_stiffness_matrix(const std::vector<Displacement> &disp) override;

    LocalMatrixd trans_matrix();

    // Nastran方式のエッジ補正行列
    Eigen::MatrixXd WarpCorrectMatrix1a();
    // エネルギー原理によるエッジ補正行列
    Eigen::MatrixXd WarpCorrectMatrix1b();
    // 法線方向モーメント補正行列
    Eigen::MatrixXd WarpCorrectMatrix2();

public:
    // static constexpr int node_num = 4;
    Node *Nodes[node_num];
    std::vector<Node *> NodesList() override { return std::vector<Node *>{Nodes[0], Nodes[1], Nodes[2], Nodes[3]}; }
    // Plane plane;
    // Thickness thickness;
    QuadPlaneElement plane_element;

    QuadPlateElement() {};
    QuadPlateElement(Node *n0, Node *n1, Node *n2, Node *n3, double t, Material mat);
    QuadPlateElement(Node *n0, Node *n1, Node *n2, Node *n3, Thickness t, Material mat);
    QuadPlateElement(int _id, Node *n0, Node *n1, Node *n2, Node *n3, double t, Material mat)
        : QuadPlateElement(n0, n1, n2, n3, t, mat)
    {
        id = _id;
    };
    QuadPlateElement(int _id, Node *n0, Node *n1, Node *n2, Node *n3, Thickness t, Material mat)
        : QuadPlateElement(n0, n1, n2, n3, t, mat)
    {
        id = _id;
    };

    /// <summary>
    /// 集中質量マトリクスを計算
    /// </summary>
    Eigen::VectorXd NodeLumpedMass()
    {
        return Eigen::VectorXd::Constant(node_num, Area() * thickness.weight_thick * Mat.dense / node_num);
    }

    Eigen::MatrixXd NodeConsistentMass();
    // Eigen::MatrixXd NodeConsistentMass2();

    std::vector<NodeLoadData> InertialForceToNodeLoadData(Eigen::Vector3d accel_vec);
    std::vector<NodeLoadData> AreaForceToNodeLoadData(std::vector<Vector> load_vecs) override;

    double Area();
    ElementType Type() { return type; }
    int NodeNum() override { return node_num; }
    int TotalDof() { return total_dof; }
    Eigen::MatrixXd StiffnessMatrix();
    Eigen::MatrixXd GeometricStiffnessMatrix(const std::vector<Displacement> &disp);

    // 剛性行列を組み込む
    void AssembleMatrix(Eigen::SparseMatrix<double> &mat, Eigen::MatrixXd K);

    // 剛性行列を組み込む
    void AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat) override;

    // 幾何剛性行列を組み込む
    void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp) override;

    // 集中質量行列を組み込む
    void AssembleMassMatrix(Eigen::SparseMatrix<double> &mat);

    // 応力を計算
    PlateStressData stress(
        Displacement d0, Displacement d1, Displacement d2, Displacement d3, double xi, double eta);
    // void shearstress(Displacement d0, Displacement d1,
    //     Displacement d2, Displacement d3, double xi, double eta);
};




#endif
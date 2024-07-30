#ifndef _ELEMENT_
#define _ELEMENT_

#ifndef SWIGCSHARP
#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif

#include "Components.h"
// #include "Model.h"

enum ElementType
{
    None, Beam, Truss, Membrane, Plate
};

class ElementBase
{
private:
    //static constexpr ElementType type = ElementType::None;
public:
    Material* Mat;
    // ElementDataBase *data;
    // int mid;
    // int TotalDOF = 0;
    int id = -1;
    ElementBase() {};
    //ElementBase(ElementType t) : type(t) {};
    // ElementBase(int mid, SSModel *model) : mid(mid), Model(model){};

    // Material *get_material();
    // virtual void check_element() = 0;
    // virtual void AssembleMatrix(double *matrix, int *dof_map, int dof_num) = 0;
    // virtual int DOFIdx(int i) = 0;
    virtual bool hasRotate() { return false; }
    //virtual ElementType Type() { return type; }
    virtual ElementType Type() { return ElementType::None; }
    virtual int TotalDof() = 0;
    virtual Eigen::MatrixXd StiffnessMatrix() = 0;
    virtual void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat) = 0;
};

struct BeamStressData {
public:
    double Nx, Qy, Qz, Mx, My, Mz;

    BeamStressData(double nx, double qy, double qz, double mx, double my, double mz) :
        Nx(nx), Qy(qy), Qz(qz), Mx(mx), My(my), Mz(mz) {};

    friend std::ostream& operator<<(std::ostream& os, const BeamStressData& bsd);
};

std::ostream& operator<<(std::ostream& os, const BeamStressData& bsd);


// 端部応力を格納。分布荷重とかも考える。
struct BeamStress {
public:
    BeamStressData S0, S1;

    BeamStress(BeamStressData s0, BeamStressData s1) :
        S0(s0), S1(s1) {};

    friend std::ostream& operator<<(std::ostream& os, const BeamStress& bsd);
};

std::ostream& operator<<(std::ostream& os, const BeamStress& bsd);

class TrussElement : public ElementBase
{
private:
    //static constexpr ElementType type = ElementType::Truss;
    static const int total_dof = 6;
    static const int node_num = 2;
    Eigen::Matrix<double, node_num, total_dof> trans_matrix();
    Eigen::MatrixXd stiffness_matrix_local();
    //double element_length();
public:
    Node* Nodes[2];
    Section* Sec;

    TrussElement() {}
    TrussElement(Node* n0, Node* n1, Section* sec, Material* mat);
    TrussElement(int _id, Node* n0, Node* n1, Section* sec, Material* mat)
        : TrussElement(n0, n1, sec, mat) {
        id = _id;
    };

    Eigen::MatrixXd StiffnessMatrix();
    bool hasRotate() { return false; }
    ElementType Type() { return ElementType::Truss; }
    //ElementType Type() { return type; }
    int TotalDof() { return total_dof; }
    void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);

     BeamStress stress(Displacement d0, Displacement d1);
};

class BeamElement : public ElementBase
{
private:
    //static constexpr ElementType type = ElementType::Beam;
    static const int total_dof = 12;
    Eigen::MatrixXd trans_matrix();
    Eigen::MatrixXd stiffness_matrix_local();
    double element_length();
public:
    double Beta;
    Node* Nodes[2];
    Section* Sec;
    
    BeamElement() {}
    //BeamElement() : ElementBase(ElementType::Beam) {}
    BeamElement(Node *n0, Node *n1, Section* sec, Material* mat, double beta=0);
    BeamElement(int _id, Node* n0, Node* n1, Section* sec, Material* mat, double beta = 0)
        :BeamElement(n0, n1, sec, mat, beta) {
        id = _id;
    };

    Eigen::MatrixXd StiffnessMatrix();
    bool hasRotate() { return true; }
    //ElementType Type() { return type; }
    ElementType Type() { return ElementType::Beam; }
    int TotalDof() { return total_dof; }
    void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);

    BeamStress stress(Displacement d0, Displacement d1);
};


struct MembraneStressData {
public:
    double sigx, sigy, sigxy;

    MembraneStressData(double sigx, double sigy, double sigxy) :
        sigx(sigx), sigy(sigy), sigxy(sigxy) {};

    friend std::ostream& operator<<(std::ostream& os, const MembraneStressData& strs);
};

std::ostream& operator<<(std::ostream& os, const MembraneStressData& strs);

class TriPlaneElement : public ElementBase
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

    Eigen::MatrixXd trans_matrix();

public:
    Node* Nodes[3];
    Plane plane;
    double thickness;

    TriPlaneElement() {};
    TriPlaneElement(Node* n0, Node* n1, Node* n2, double t, Material* mat);
    TriPlaneElement(int _id, Node* n0, Node* n1, Node* n2, double t, Material* mat)
        :TriPlaneElement(n0, n1, n2, t, mat) {
        id = _id;
    };

    double Area();
    ElementType Type() { return type; }
    int TotalDof() { return total_dof; }
    Eigen::MatrixXd StiffnessMatrix(); 
    void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);

    MembraneStressData stress(Displacement d0, Displacement d1, Displacement d2);
};

class QuadPlaneElement : public ElementBase
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
    Eigen::MatrixXd trans_matrix();

public:
    Node* Nodes[4];
    Plane plane;
    double thickness;

    QuadPlaneElement() {};
    QuadPlaneElement(Node* n0, Node* n1, Node* n2, Node* n3, double t, Material* mat);
    QuadPlaneElement(int _id, Node* n0, Node* n1, Node* n2, Node* n3, double t, Material* mat) :
        QuadPlaneElement(n0, n1, n2, n3, t, mat) {
        id = _id;
    };

    // double Area();
    ElementType Type() { return type; }
    int TotalDof() { return total_dof; }
    Eigen::MatrixXd StiffnessMatrix();
    void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);

    MembraneStressData stress(Displacement d0, Displacement d1, 
        Displacement d2, Displacement d3, double xi, double eta);
};

struct PlateStressData {
public:
    double Mx, My, Mxy, Qx, Qy;

    PlateStressData(double mx, double my, double mxy, double qx, double qy) :
        Mx(mx), My(my), Mxy(mxy), Qx(qx), Qy(qy) {};

    friend std::ostream& operator<<(std::ostream& os, const PlateStressData& strs);
};

std::ostream& operator<<(std::ostream& os, const PlateStressData& strs);

//struct PlateStress {
//    MembraneStressData plane_stress;
//    PlateStressData plate_stress;
//
//    PlateStress(Displacement d0, Displacement d1, Displacement d2);
//};

class TriPlateElement : public ElementBase
{
private:
    static constexpr int total_dof = 18;
    static constexpr int node_dof = 6;
    static constexpr ElementType type = ElementType::Plate;
    using LocalMatrixd = Eigen::Matrix<double, total_dof, total_dof>;

    Eigen::MatrixXd BMatrix(double xi, double eta);
    
    Eigen::Matrix3d DMatrix();
    Eigen::MatrixXd localStiffnessMatrix();
    LocalMatrixd trans_matrix();

public:
    Node* Nodes[3];
    Plane plane;
    double thickness;
    TriPlaneElement plane_element;

    TriPlateElement() {};
    TriPlateElement(Node* n0, Node* n1, Node* n2, double t, Material* mat);
    TriPlateElement(int _id, Node* n0, Node* n1, Node* n2, double t, Material* mat)
        :TriPlateElement(n0, n1, n2, t, mat) {
        id = _id;
    };

    double Area();
    ElementType Type() { return type; }
    int TotalDof() { return total_dof; }
    Eigen::MatrixXd StiffnessMatrix();
    void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
    PlateStressData stress(
        Displacement d0, Displacement d1, Displacement d2, double xi, double eta);
    //void shearstress(Displacement d0, Displacement d1, Displacement d2);
};


class QuadPlateElement : public ElementBase
{
public:
    // 下のpublicの範囲でこの定数を使っているからこちらもpublicにしておく必要があるらしい。
    static constexpr int node_num = 4;
private:

    static constexpr int node_dof = 6;
    // static constexpr int node_num = 4;
    static constexpr int node_dof_local = 3;
    static constexpr int total_dof = node_dof * node_num;
    static constexpr int total_dof_local = node_dof_local * node_num;

    using LocalMatrixd = Eigen::Matrix<double, total_dof, total_dof>;

    static constexpr ElementType type = ElementType::Plate;

    Eigen::Matrix2d JMatrix(double xi, double eta);
    
    Eigen::Matrix2d dJinv_dxi(double xi, double eta);
    Eigen::Matrix2d dJinv_deta(double xi, double eta);

    Eigen::MatrixXd HVecs(Eigen::VectorXd shape_funcs);
    Eigen::MatrixXd BMatrix(double xi, double eta);

    Eigen::Matrix3d DMatrix();
    Eigen::MatrixXd localStiffnessMatrix();

    LocalMatrixd trans_matrix();

public:
    // static constexpr int node_num = 4;
    Node* Nodes[node_num];
    Plane plane;
    double thickness;
    QuadPlaneElement plane_element;

    QuadPlateElement() {};
    QuadPlateElement(Node* n0, Node* n1, Node* n2, Node* n3, double t, Material* mat);
    QuadPlateElement(int _id, Node* n0, Node* n1, Node* n2, Node* n3, double t, Material* mat)
        :QuadPlateElement(n0, n1, n2, n3, t, mat) {
        id = _id;
    };

    // double Area();
    ElementType Type() { return type; }
    int TotalDof() { return total_dof; }
    Eigen::MatrixXd StiffnessMatrix();
    void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
    PlateStressData stress(
        Displacement d0, Displacement d1, Displacement d2, Displacement d3, double xi, double eta);
    //void shearstress(Displacement d0, Displacement d1, 
    //    Displacement d2, Displacement d3, double xi, double eta);

};

//class TrussElement : public ElementBase {
//private:
//    const ElementType type = ElementType::Truss;
//    static const int total_dof = 6;
//    Eigen::MatrixXd stiffness_matrix_local();
//
//public:
//    Node* Nodes[2];
//    Section* Sec;
//
//    TrussElement() {};
//    TrussElement(Node* n0, Node* n1, Section* sec, Material* mat);
//
//    Eigen::MatrixXd StiffnessMatrix();
//    bool hasRotate() { return false; }
//    ElementType Type() { return type; }
//    int TotalDof() { return total_dof; }
//    void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
//};

#endif
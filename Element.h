#ifndef _ELEMENT_
#define _ELEMENT_

#ifndef SWIGCSHARP

//#ifdef USE_MKL
//#define EIGEN_USE_MKL_ALL
//#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif

#include "Components.h"
//#include "LoadComponent.h"
// #include "Model.h"

Eigen::Matrix3d trans_matrix3(const Point p0, const Point p1, const double beta);
Eigen::Matrix3d trans_matrix3(const Plane plane);

enum class ElementType
{
    None, Beam, Truss, Membrane, Plate, DKT, DKQ
};

class ElementBase
{
protected:
    // static Eigen::Matrix3d trans_matrix3(const Point p0, const Point p1, const double beta);
    // static Eigen::Matrix3d trans_matrix3(const Plane plane);

private:
    //static constexpr ElementType type = ElementType::None;
    virtual Eigen::MatrixXd geometric_local_stiffness_matrix(const std::vector<Displacement> &disp) = 0;

public:
    Material Mat;
    // ElementDataBase *data;
    // int mid;
    // int TotalDOF = 0;
    int id = -1;
    ElementBase() {};
    virtual int NodeNum() = 0;
    virtual std::vector<Node*> NodesList() = 0;
    //ElementBase(ElementType t) : type(t) {};
    // ElementBase(int mid, SSModel *model) : mid(mid), Model(model){};

    // Material *get_material();
    // virtual void check_element() = 0;
    // virtual void AssembleMatrix(double *matrix, int *dof_map, int dof_num) = 0;
    // virtual int DOFIdx(int i) = 0;

	virtual std::vector<NodeLoadData> InertialForceToNodeLoadData(Eigen::Vector3d accel_vec) = 0;
    virtual Eigen::VectorXd NodeLumpedMass() = 0;
    virtual Eigen::MatrixXd NodeConsistentMass() = 0;
    virtual bool hasRotate() { return false; }
    //virtual ElementType Type() { return type; }
    virtual ElementType Type() { return ElementType::None; }
    virtual int TotalDof() = 0;
    virtual Eigen::MatrixXd StiffnessMatrix() = 0;
    virtual void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat) = 0;
    virtual void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double>& mat, const std::vector<Displacement>& disp) = 0;
    virtual void AssembleMassMatrix(Eigen::SparseMatrix<double>& mat) = 0;
};


#endif
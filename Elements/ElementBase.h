#ifndef _ELEMENTBASE_
#define _ELEMENTBASE_

#ifndef SWIG
#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif

#include "../Components.h"

Eigen::Matrix3d trans_matrix3(const Point p0, const Point p1, const double beta);
Eigen::Matrix3d trans_matrix3(const Plane plane);

enum class ElementType
{
    None = 0,
    Beam = 1,
    Truss = 2,
    Membrane = 3, // 総称（新規要素では未使用＝レガシー扱い）
    Plate = 4,    // 未使用(予約)
    DKT = 5,      // TriPlateElement
    DKQ = 6,      // QuadPlateElement
    // --- 追加（末尾に追記し既存値は変更しない） ---
    ComplexBeam = 7,  // ComplexBeamElement
    TriMembrane = 8,  // TriPlaneElement
    QuadMembrane = 9, // QuadPlaneElement
};

// 要素種別の分類判定（分類の定義をここ1箇所に集約する）
inline bool IsBarType(ElementType t)
{
    return t == ElementType::Beam || t == ElementType::Truss || t == ElementType::ComplexBeam;
}

inline bool IsPlaneType(ElementType t)
{
    return t == ElementType::Membrane || t == ElementType::TriMembrane || t == ElementType::QuadMembrane || t == ElementType::Plate || t == ElementType::DKT || t == ElementType::DKQ;
}

class ElementBase
{

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
    [[deprecated("Use GetStiffnessTriplets() instead. This coeffRef-based method is less efficient and will be removed in a future version.")]]
    virtual void AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat) = 0;
    [[deprecated("Use GetGeometricStiffnessTriplets() instead. This coeffRef-based method is less efficient and will be removed in a future version.")]]
    virtual void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double>& mat, const std::vector<Displacement>& disp) = 0;
    virtual void AssembleMassMatrix(Eigen::SparseMatrix<double>& mat) = 0;

    // Triplet方式での行列組立用メソッド
    virtual void GetStiffnessTriplets(std::vector<Eigen::Triplet<double>>& triplets) = 0;
    virtual void GetGeometricStiffnessTriplets(const std::vector<Displacement>& disp,
                                                std::vector<Eigen::Triplet<double>>& triplets) = 0;
};


#endif
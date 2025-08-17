#ifndef _TRUSS_ELEMENT_
#define _TRUSS_ELEMENT_

#include "BarElement.h"

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

#endif
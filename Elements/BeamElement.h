#ifndef _BEAM_ELEMENT_
#define _BEAM_ELEMENT_

#include "BarElement.h"

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

#endif
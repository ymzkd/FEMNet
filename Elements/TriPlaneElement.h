#ifndef _TRI_PLANE_ELEMENT_
#define _TRI_PLANE_ELEMENT_

#include "PlaneElement.h"

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

#endif
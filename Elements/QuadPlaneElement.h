#ifndef _QUAD_PLANE_ELEMENT_
#define _QUAD_PLANE_ELEMENT_

#include "PlaneElement.h"

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

#endif
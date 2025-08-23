#ifndef _TRI_PLATE_ELEMENT_
#define _TRI_PLATE_ELEMENT_

#include "PlaneElement.h"
#include "TriPlaneElement.h"

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
    TriPlateElement(Node *n0, Node *n1, Node *n2, Thickness t, Material mat, double beta = 0);
    TriPlateElement(Node *n0, Node *n1, Node *n2, double t, Material mat, double beta = 0);
    TriPlateElement(int _id, Node *n0, Node *n1, Node *n2, Thickness t, Material mat, double beta = 0)
        : TriPlateElement(n0, n1, n2, t, mat, beta)
    {
        id = _id;
    };
    TriPlateElement(int _id, Node *n0, Node *n1, Node *n2, double t, Material mat, double beta = 0)
        : TriPlateElement(n0, n1, n2, t, mat, beta)
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

#endif
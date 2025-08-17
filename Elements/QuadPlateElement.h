#ifndef _QUAD_PLATE_ELEMENT_
#define _QUAD_PLATE_ELEMENT_

#include "PlaneElement.h"
#include "QuadPlaneElement.h"   

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
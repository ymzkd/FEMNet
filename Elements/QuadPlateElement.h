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
    QuadPlateElement(Node *n0, Node *n1, Node *n2, Node *n3, double t, Material mat, double beta = 0);
    QuadPlateElement(Node *n0, Node *n1, Node *n2, Node *n3, Thickness t, Material mat, double beta = 0);
    QuadPlateElement(int _id, Node *n0, Node *n1, Node *n2, Node *n3, double t, Material mat, double beta = 0)
        : QuadPlateElement(n0, n1, n2, n3, t, mat, beta)
    {
        id = _id;
    };
    QuadPlateElement(int _id, Node *n0, Node *n1, Node *n2, Node *n3, Thickness t, Material mat, double beta = 0)
        : QuadPlateElement(n0, n1, n2, n3, t, mat, beta)
    {
        id = _id;
    };

    /// <summary>
    /// 集中質量マトリクスを計算
    /// </summary>
    Eigen::VectorXd NodeLumpedMass()
    {
        return Eigen::VectorXd::Constant(node_num, Area() * thickness.WeightThickness() * Mat.dense / node_num);
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
    [[deprecated("Use GetStiffnessTriplets() instead. This coeffRef-based method is less efficient and will be removed in a future version.")]]
    void AssembleStiffMatrix(Eigen::SparseMatrix<double> &mat) override;

    // 幾何剛性行列を組み込む
    [[deprecated("Use GetGeometricStiffnessTriplets() instead. This coeffRef-based method is less efficient and will be removed in a future version.")]]
    void AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double> &mat, const std::vector<Displacement> &disp) override;

    // 集中質量行列を組み込む
    void AssembleMassMatrix(Eigen::SparseMatrix<double> &mat);

    // Triplet方式での行列組立
    void GetStiffnessTriplets(std::vector<Eigen::Triplet<double>>& triplets) override;
    void GetGeometricStiffnessTriplets(const std::vector<Displacement>& disp,
                                        std::vector<Eigen::Triplet<double>>& triplets) override;

    // 応力を計算
    PlateStressData stress(
        Displacement d0, Displacement d1, Displacement d2, Displacement d3, double xi, double eta);
    // void shearstress(Displacement d0, Displacement d1,
    //     Displacement d2, Displacement d3, double xi, double eta);

    // 整合節点力 f = K_e * u_e を計算する。
    // d0..d3 はグローバル座標の節点変位。
    // local=true のとき各節点の力・モーメントを要素plane軸へ回転して返す。
    std::vector<NodeLoadData> NodalForces(
        Displacement d0, Displacement d1, Displacement d2, Displacement d3, bool local);
};

#endif
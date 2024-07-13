#ifndef _MODEL_
#define _MODEL_

#ifndef SWIGCSHARP
#include <vector>
#include <map>
#include <array>

#include <Eigen/Sparse>
// #include <Eigen/Dense>
#endif

class ElementHandle {
public:
    ElementBase* data;
    ElementHandle(ElementBase* data) : data(data) {};

    // operator-> のオーバーロード
    ElementBase* operator->() {
        return data;
    }

    const ElementBase* operator->() const {
        return data;
    }
};


class SSModel
{
private:
    std::vector<double> displacement;
    // Eigen::MatrixXd stiff_matrix_xd(BeamElementData data);
    // Eigen::MatrixXd stiff_matrix_xd(TrussElementData data);
    std::vector<int> FreeIndices();
    Eigen::SparseMatrix<double> AssembleMatrix();
    static Eigen::SparseMatrix<double> extractSubMatrix(
        const Eigen::SparseMatrix<double>& mat,
        const std::vector<int>& rowIndices,
        const std::vector<int>& colIndices);

public:
    int NodeNum = 0;
    std::vector<Node> Nodes;
    std::vector<Material> Materials;
    std::vector<Section> Sections;
    std::vector<ElementHandle> Elements;

    void add_element(ElementBase* data);
    void add_element(BeamElement *data);
    void add_element(TrussElement *data);
    void add_element(TriPlaneElement *data);
    void add_element(TriPlateElement *data);
    void add_element(QuadPlaneElement *data);
    void add_element(QuadPlateElement *data);

    std::vector<Displacement> Solve(std::list<Load> loads);

    // void sample_function1();
    // Eigen::MatrixXd sample_function2();
    // std::map<int, std::array<bool, 6>> Supports;
    // std::map<int, std::array<double, 6>> Loads;

    // std::vector<int> DOFMap;

    // int MatrixDim();

    // SSModel(){};
    // SSModel(int node_num);

    // void check_model();

    // void Analysis(double *result);
    // void AnalysisXd(double *result);
    // void AnalysisSparse(double *result);
    // void Solve();

    // void SolveSparse();

    // std::array<double, 6> getNodeDisp(int idx);
};

#endif
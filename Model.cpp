#include <unordered_map>

#include <Eigen/Sparse>

#include "Element.h"
#include "Components.h"
#include "Model.h"

Eigen::SparseMatrix<double> SSModel::extractSubMatrix(
    const Eigen::SparseMatrix<double>& mat,
    const std::vector<int>& rowIndices,
    const std::vector<int>& colIndices) {

    // 新しい疎行列を作成
    Eigen::SparseMatrix<double> subMat(rowIndices.size(), colIndices.size());

    // インデックス変換用のマップ
    std::unordered_map<int, int> rowMap;
    for (int i = 0; i < rowIndices.size(); ++i)
        rowMap[rowIndices[i]] = i;

    std::unordered_map<int, int> colMap;
    for (int j = 0; j < colIndices.size(); ++j)
        colMap[colIndices[j]] = j;

    // 指定された行と列の要素を新しい疎行列にコピー
    for (int k = 0; k < mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
            auto rowMapItem = rowMap.find(it.row());
            auto colMapItem = colMap.find(it.col());

            if (rowMapItem != rowMap.end() && colMapItem != colMap.end())
                subMat.insert(rowMapItem->second, colMapItem->second) = it.value();
        }
    }

    return subMat;
}

std::vector<int> SSModel::FreeIndices()
{
    std::vector<int> indices;
    int idx = 0;
    for each (Node n in Nodes)
        for each(bool f in n.Fix.flags) {
            if (!f) indices.push_back(idx);
            idx++;
        }
    return indices;
}

void SSModel::add_element(ElementBase* data) {
    Elements.push_back(ElementHandle(data));
    Elems.push_back(data);
}

//void SSModel::add_element(BeamElement *data){
//
//    Elements.push_back(ElementHandle(data));
//}
//
//void SSModel::add_element(TrussElement* data)
//{
//    Elements.push_back(ElementHandle(data));
//}
//
//void SSModel::add_element(TriPlaneElement* data)
//{
//    Elements.push_back(ElementHandle(data));
//}
//
//void SSModel::add_element(TriPlateElement* data)
//{
//    Elements.push_back(ElementHandle(data));
//}
//
//void SSModel::add_element(QuadPlaneElement* data)
//{
//    Elements.push_back(ElementHandle(data));
//}
//
//void SSModel::add_element(QuadPlateElement* data)
//{
//    Elements.push_back(ElementHandle(data));
//}

Eigen::SparseMatrix<double> SSModel::AssembleMatrix()
{
    int mat_size = Nodes.size() * 6;
    Eigen::SparseMatrix<double> mat(mat_size, mat_size);
    for each (ElementHandle eh in Elements)
        eh->AssembleStiffMatrix(mat);

    return mat;
}

std::vector<Displacement> SSModel::Solve(std::list<Load> loads)
{
    Eigen::VectorXd f(Nodes.size() * 6);
    f.setZero();
    for each (Load l in loads)
    {
        if (l.id < 0) continue;

        int pos = l.id * 6;
        f[pos] = l.Px();
        f[pos+1] = l.Py();
        f[pos+2] = l.Pz();
        f[pos+3] = l.Mx();
        f[pos+4] = l.My();
        f[pos+5] = l.Mz();
    }

    std::vector<int> indices = FreeIndices();
    Eigen::SparseMatrix<double> m = SSModel::extractSubMatrix(AssembleMatrix(), indices, indices);
    
    Eigen::VectorXd f_fix(indices.size());
    for (size_t i = 0; i < indices.size(); i++)
        f_fix(i) = f(indices[i]);
    
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
    solver.compute(m);
    Eigen::VectorXd d_fix = solver.solve(f_fix);
    
    Eigen::VectorXd d(Nodes.size() * 6);
    d.setZero();
    for (size_t i = 0; i < indices.size(); i++)
        d(indices[i]) = d_fix(i);

    std::vector<Displacement> disp;
    for (size_t i = 0; i < Nodes.size(); i++)
    {
        int pos = i * 6;
        disp.push_back(Displacement(d[pos], d[pos + 1], d[pos + 2], d[pos + 3], d[pos + 4], d[pos + 5]));
    }
    
    return disp;
}

std::vector<BeamElement*> SSModel::beam_elements()
{
    std::vector<BeamElement*> beams;
    for each (ElementBase* e in Elems)
    {
        if (e->Type() == ElementType::Beam)
            beams.push_back(dynamic_cast<BeamElement*>(e));
    }
    return beams;
}


//void SSModel::sample_function1()
//{
//    Eigen::MatrixXd m(2, 2);
//    m << 1, 2,
//        3, 4;
//    std::cout << m << std::endl;
//}

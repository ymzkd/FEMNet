#ifndef _MODEL_
#define _MODEL_

#ifndef SWIGCSHARP
#include <vector>
#include <map>
#include <array>


// #include <Eigen/Dense>

#ifdef EIGEN_USE_MKL_ALL
    //#define EIGEN_USE_MKL_ALL
    #include <Eigen/Sparse>
    #include <Eigen/PardisoSupport>
#else
    #include <Eigen/Sparse>
#endif

#endif
#include "LoadComponent.h"

class FEModel
{
private:
    const double GRAVACCEL = 9806.6;
    std::vector<int> FreeIndices();
    std::vector<int> FixIndices();
    
    /// <summary>
    /// 集中質量マトリクスの対象でない(並進以外の回転自由度等)または固定自由度の
    /// 全体自由度におけるインデックスを格納した配列を返す関数
    /// </summary>
    /// <returns></returns>
    std::vector<int> UnLumpedFixIndices();

    Eigen::SparseMatrix<double> AssembleStiffnessMatrix();
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="A"></param>
    /// <param name="fixed_indices"></param>
    /// <param name="free_matrix"></param>
    /// <param name="free_fixed_matrix"></param>
    /// <param name="fixed_matrix"></param>
    static void splitMatrixWithResize(
        const Eigen::SparseMatrix<double>& A, 
        const std::vector<int>& fixed_indices, 
        Eigen::SparseMatrix<double>& free_matrix, 
        Eigen::SparseMatrix<double>& free_fixed_matrix, 
        Eigen::SparseMatrix<double>& fixed_matrix);

    /// <summary>
    /// ┌      ┐
    /// │ A  B │
    /// │ C  D │
    /// └      ┘

    /// │ A  B │
    /// │ C  D │

    /// </summary>
    /// <param name="A"></param>
    /// <param name="fixed_indices"></param>
    /// <param name="free_matrix"></param>
    static void splitMatrixWithResize(
        const Eigen::SparseMatrix<double>& A, 
        const std::vector<int>& fixed_indices, 
        Eigen::SparseMatrix<double>& free_matrix);
    
    static Eigen::SparseMatrix<double> extractSubMatrix(
        const Eigen::SparseMatrix<double>& mat,
        const std::vector<int>& rowIndices,
        const std::vector<int>& colIndices);

public:
    

    int NodeNum() { return Nodes.size(); }
    int DOFNum() { return NodeNum() * 6; }

    std::vector<Node> Nodes;
    std::vector<Material> Materials;
    std::vector<Section> Sections;
    
    std::vector<std::shared_ptr<ElementBase>> Elements;
    
    // template<>
    // void add_element(ElementBase* data);
    //void add_element(ElementBase data);
    void add_element(BeamElement data);
    void add_element(ComplexBeamElement data);
    void add_element(TrussElement data);
    void add_element(TriPlaneElement data);
    void add_element(TriPlateElement data);
    void add_element(QuadPlaneElement data);
    void add_element(QuadPlateElement data);

    // index based element addition
    void add_truss_element(int id, int n1_id, int n2_id, int sec_id, int mat_id);
    void add_beam_element(int id, int n1_id, int n2_id, int sec_id, int mat_id, double beta);
    //void add_tri_plane_element(int id, int n1_id, int n2_id, int n3_id, int mat_id, int sec_id);
    void add_tri_plate_element(int id, int n1_id, int n2_id, int n3_id, double thickness, int mat_id);
    //void add_quad_plane_element(int id, int n1_id, int n2_id, int n3_id, int n4_id, int mat_id, int sec_id);
    void add_quad_plate_element(int id, int n1_id, int n2_id, int n3_id, int n4_id, double thickness, int mat_id);
    
    BarElementBase* GetBarElement(int id);
    BeamElement* GetBeamElement(int id);
    TrussElement* GetTrussElement(int id);

    [[deprecated("This function is deprecated. Please use SolveLinearStatic instead.")]]
    void Solve(
        std::vector<std::shared_ptr<LoadBase>>& loads,
        std::vector<Displacement>& disp,
        std::vector<NodeLoad>& react);

    void SolveLinearStatic(
        std::vector<std::shared_ptr<LoadBase>>& loads,
        std::vector<Displacement>& disp,
        std::vector<NodeLoad>& react);

    //void SolveVibrationTest();

    int SolveVibration(const int nev, std::vector<double>& eigen_values,
        std::vector<std::vector<Displacement>>& mode_vectors);

};

typedef FEModel FEModel;

// Static Solver Package
class FEStaticResult {
public:
    std::shared_ptr<FEModel> model;
    std::vector<std::shared_ptr<LoadBase>> loads;
    std::vector<Displacement> displace;
    std::vector<NodeLoad> react_force;

    FEStaticResult(
        std::shared_ptr<FEModel> model,
        std::vector<std::shared_ptr<LoadBase>> loads,
        std::vector<Displacement> displace,
        std::vector<NodeLoad> react_force)
        : model(model), loads(loads),
        displace(displace), react_force(react_force) {};

    BeamStressData GetBeamStress(int eid, double p);

    Displacement GetBeamDisplace(int eid, double p);
};

class FEVibrateResult {
public:
    std::shared_ptr<FEModel> model;
    
    int modes_num() { return eigs.size(); };
    std::vector<double> eigs;
    std::vector<std::vector<Displacement>> mode_vectors;

    FEVibrateResult(
        std::shared_ptr<FEModel> model,
        std::vector<std::vector<Displacement>> mode_vectors,
        std::vector<double> eigs)
        : model(model), mode_vectors(mode_vectors),
        eigs(eigs) {};

};

#endif
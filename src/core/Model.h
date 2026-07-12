#ifndef _MODEL_
#define _MODEL_

#ifndef SWIG
#include <vector>
#include <map>
#include <array>
#include <memory>
#include <string>

#include "SparseSolver.h"

#endif

#include "LoadComponent.h"
#include "RigidLink.h"
#ifndef SWIG
#include "SparseMatrixUtils.h"
#endif

#define PI 3.141592653589793238462643

// 前方宣言
class FEDynamicDampInitializer;
class FEDynamicStiffDampInitializer;
class DynamicAnalysis;

class IResponseSpectrum {
public:
    virtual double Acceleration(double t) = 0;
    virtual double Velocity(double t) = 0;
	virtual double Displacement(double t) = 0;
};

class FEModel
{
private:
    /// <summary>
    /// 集中質量マトリクスの対象でない(並進以外の回転自由度等)または固定自由度の
    /// 全体自由度におけるインデックスを格納した配列を返す関数
    /// </summary>
    /// <returns></returns>
    std::vector<int> UnLumpedFixIndices();

	// 剛性マトリクスの組み立て
    Eigen::SparseMatrix<double> AssembleStiffnessMatrix();

	// 質量マトリクスの組み立て
    Eigen::SparseMatrix<double> AssembleMassMatrix();

    Eigen::SparseMatrix<double> AssembleGeometricStiffnessMatrix(
        const std::vector<Displacement> &displacements);

	friend class DynamicAnalysis;
    friend class FEBucklingAnalysis;
	friend class FEVibrateResult;
    friend class ResponseSpectrumMethod;

public:
    FEModel();

    double GraityAccel = 9806.65;

    /// <summary>
    /// 非拘束自由度の全自由度におけるインデックスを格納した配列を返す関数
    /// </summary>
    std::vector<int> FreeIndices(bool rigid_link = false);
    

    /// <summary>
    /// 剛体連結されている自由度の全自由度におけるインデックスを格納した配列を返す関数
    /// </summary>
    std::vector<int> SlaveIndices();

    /// <summary>
    /// 拘束自由度の全自由度におけるインデックスを格納した配列を返す関数
    /// </summary>
    std::vector<int> FixIndices();

    int NodeNum() { return Nodes.size(); }
    int DOFNum() { return NodeNum() * 6; }
	int FreeDOFNum() { return FreeIndices().size(); }
	int FixedDOFNum() { return FixIndices().size(); }

    std::vector<Node> Nodes;
    std::vector<Material> Materials;
    std::vector<Section> Sections;
    
    std::vector<std::shared_ptr<ElementBase>> Elements;
    std::shared_ptr<RigidLinks> RigidLinkData;
    
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
    void add_tri_plate_element(int id, int n1_id, int n2_id, int n3_id, double thickness, int mat_id);
    void add_quad_plate_element(int id, int n1_id, int n2_id, int n3_id, int n4_id, double thickness, int mat_id);
    
    BarElementBase* GetBarElement(int id);
    BeamElement* GetBeamElement(int id);
    TrussElement* GetTrussElement(int id);
    QuadPlateElement* GetQuadPlateElement(int id);
    TriPlateElement* GetTriPlateElement(int id);

    // 慣性力を等価な節点荷重に変換
    std::vector<NodeLoadData> InnertialForceToNodeLoads(
		const InertialForce inertial_force);

    double SumNodeMass();

    // 要素質量に基づく節点質量を計算してセットアップ
    void ComputeElementNodeMass();

    void SolveLinearStatic(
        std::vector<std::shared_ptr<LoadBase>>& loads,
        std::vector<Displacement>& disp,
        std::vector<NodeLoad>& react);

    int SolveVibration(const int nev, std::vector<double>& eigen_values,
        std::vector<std::vector<Displacement>>& mode_vectors);

    // 構造モデルのテキスト形式ファイル入出力
    // Save: 現在のモデルを path に書き出す
    // Load: path からモデルを読み込み、現在のモデルを置き換える
    void Save(const std::string& path);
    void Load(const std::string& path);

};

// typedef FEModel FEModel;




#endif
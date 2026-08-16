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
    // 減衰の考慮
    std::shared_ptr<FEDynamicDampInitializer> DampInitializer = nullptr;

    // このスペクトルが表現している減衰比(基準減衰)。告示・指針の設計スペクトルは5%。
    // DampingCorrectionFactorの正規化基準であり、Fh補正無効時の実効減衰でもある。
    double BaseDampingFactor = 0.05;

    // 減衰の影響係数計算
    bool enable_damp_factor = false;
    virtual double DampingCorrectionFactor(const double t);

    // 周期tのスペクトル値が実際に対応する減衰比を返す。
    //   Fh補正無効時: 基準減衰(BaseDampingFactor)
    //   Fh補正有効時: DampInitializerのモード減衰(算定不能なら基準減衰にフォールバック)
    // CQCの相関係数など、スペクトル値と整合した減衰比が必要な箇所はこれを参照する。
    double effective_damping_rate(const double t);

    double acceleration_factored(double t);
    double velocity_factored(double t);
    double displacement_factored(double t);

    virtual double Acceleration(double t) = 0;
    virtual double Velocity(double t) = 0;
	virtual double Displacement(double t) = 0;
};

// 状態依存要素(非抗圧トラス等)の解析時状態(要素インデックス→有効フラグ)。
// 状態は FEModel の要素内ではなく各解析Operatorが所有し、
// 剛性組立や応力復元の際に明示的に渡す。未登録の要素は true(規定剛性で有効)。
class ElementStates
{
private:
    std::vector<std::pair<int, bool>> states;

public:
    // 指定要素の状態を返す(未登録は true)
    bool Get(int eid) const
    {
        for (const auto &s : states)
            if (s.first == eid)
                return s.second;
        return true;
    }

    void Set(int eid, bool active)
    {
        for (auto &s : states)
        {
            if (s.first == eid)
            {
                s.second = active;
                return;
            }
        }
        states.emplace_back(eid, active);
    }

    void Clear() { states.clear(); }
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

public:
    // === 行列・ベクトル組立サービス(解析Operator向け, SWIG非公開) ===
    // FEModelは構造データの保持と組立のみを担い、解析(ソルバー)は
    // 各解析Operator(FELinearStaticOp, FEVibrationAnalysis等)が実装する。

	// 剛性マトリクスの組み立て(上三角格納)
	// states: 状態依存要素の状態(nullptrなら全要素を規定剛性で組立)
    Eigen::SparseMatrix<double> AssembleStiffnessMatrix(const ElementStates *states = nullptr);

    // 荷重リストから全体節点荷重ベクトルを組み立てる(InertialForceは要素質量から展開)
    Eigen::VectorXd AssembleLoadVector(const std::vector<std::shared_ptr<LoadBase>> &loads);

	// 質量マトリクスの組み立て(上三角格納)
    Eigen::SparseMatrix<double> AssembleMassMatrix();

    // 幾何剛性マトリクスの組み立て(上三角格納)
    Eigen::SparseMatrix<double> AssembleGeometricStiffnessMatrix(
        const std::vector<Displacement> &displacements);

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
    void add_element(TensionTrussElement data);
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

    // 構造モデルのテキスト形式ファイル入出力
    // Save: 現在のモデルを path に書き出す
    // Load: path からモデルを読み込み、現在のモデルを置き換える
    void Save(const std::string& path);
    void Load(const std::string& path);

};

// typedef FEModel FEModel;




#endif
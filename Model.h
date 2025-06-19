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

#define PI 3.141592653589793238462643

// 前方宣言
class FEDynamicDampInitializer;
class FEDynamicStiffDampInitializer;
class DynamicAnalysis;

class IResponseSpectrum {
public:
    //IResponseSpectrum() {}
    //virtual ~IResponseSpectrum() {}
    virtual double Acceleration(double t) = 0;
    virtual double Velocity(double t) = 0;
	virtual double Displacement(double t) = 0;
};

class DASampler {

public:

    //friend class DynamicAnalysis;

    int step;
    std::vector<Displacement> velocity, displacement, acceleration;
    virtual void Sampling(DynamicAnalysis& analysis) = 0;
};

class DASampler_MaxDisplacement : public DASampler {
public:
    double max_displacement = 0.0;

    void Sampling(DynamicAnalysis& da) override;
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

	friend class DynamicAnalysis;

public:
    static constexpr double GRAVACCEL = 9806.6;\

    /// <summary>
    /// 非拘束自由度の全自由度におけるインデックスを格納した配列を返す関数
    /// </summary>
    std::vector<int> FreeIndices();

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
    QuadPlateElement* GetQuadPlateElement(int id);
    TriPlateElement* GetTriPlateElement(int id);

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

// 変形ケースを表す基底クラス
class FEDeformCase {
public:
    std::shared_ptr<FEModel> model;
    //std::vector<std::shared_ptr<LoadBase>> loads;
    //std::vector<Displacement> displace;
    //std::vector<NodeLoad> react_force;

    FEDeformCase() {};
	FEDeformCase(std::shared_ptr<FEModel> model) : model(model) {};

    virtual BeamStressData GetBeamStress(int eid, double p) = 0;
    virtual PlateStressData GetPlateStressData(int eid, double xi, double eta) = 0;
    virtual Displacement GetBeamDisplace(int eid, double p) = 0;
    //FEDeformCase(std::shared_ptr<FEModel> model, std::vector<std::shared_ptr<LoadBase>> loads, std::vector<Displacement> displace, std::vector<NodeLoad> react_force)
    virtual std::vector<Displacement> GetDisplacements() = 0;
    //virtual std::vector<Displacement> GetVelocities() = 0;
    //virtual std::vector<Displacement> GetAccelerations() = 0;
};

// Static Solver Package
class FEStaticResult : public FEDeformCase {
public:
    //std::shared_ptr<FEModel> model;
    std::vector<std::shared_ptr<LoadBase>> loads;
    std::vector<Displacement> displace;
    std::vector<NodeLoad> react_force;

    FEStaticResult(
        std::shared_ptr<FEModel> model,
        std::vector<std::shared_ptr<LoadBase>> loads,
        std::vector<Displacement> displace,
        std::vector<NodeLoad> react_force)
        : FEDeformCase(model), loads(loads),
        displace(displace), react_force(react_force) {};

    BeamStressData GetBeamStress(int eid, double p) override;

    /// <summary>
    /// Obtain stress data for plate elements
    /// </summary>
    /// <param name="eid">element index</param>
    /// <param name="xi">
    /// xi for square isoparametric elements and L1 for triangular element
    /// area coordinate system in the first parameter.
    /// </param>
    /// <param name="eta">
    /// eta for square isoparametric elements and L2 for triangular element
    /// area coordinate system in the second parameter.
    /// </param>
    /// <returns>Plate element stress data</returns>
    PlateStressData GetPlateStressData(int eid, double xi, double eta) override;


    Displacement GetBeamDisplace(int eid, double p) override;

    // FEDeformCase を介して継承されました
    std::vector<Displacement> GetDisplacements() override;
};

class FEVibrateResult {
public:
    std::shared_ptr<FEModel> model;
    
    int modes_num() { return eigs.size(); };
    std::vector<double> eigs; // Omegas
    std::vector<std::vector<Displacement>> mode_vectors;

	FEVibrateResult() {};

    FEVibrateResult(
        std::shared_ptr<FEModel> model,
        std::vector<std::vector<Displacement>> mode_vectors,
        std::vector<double> eigs)
        : model(model), mode_vectors(mode_vectors),
        eigs(eigs) {};

    std::vector<double> ParticipationFactors();
    std::vector<double> ParticipationFactors(Vector direction);
    std::vector<Displacement> ParticipationDirectedFactors();
    std::vector<double> EffectiveMassRates();
    std::vector<Displacement> EffectiveDirectedMassRates();

	std::vector<double> NaturalPeriods() {
		std::vector<double> periods;
		
        for (double eig : eigs)
			periods.push_back(2.0 * PI / eig);

        return periods;
	}
};

class DARecorder {
public:
    virtual void Record(DynamicAnalysis& da) = 0;
};

class DAEnergyRecorder: public DARecorder {
public:
	std::vector<double> kinetic_energy, potential_energy, damping_energy, input_energy;

    void Initialize();
    void RecordKineticEnergy(DynamicAnalysis& da);
    void RecordPotentialEnergy(DynamicAnalysis& da);
    void RecordDampingEnergy(DynamicAnalysis& da);
    void RecordInputEnergy(DynamicAnalysis& da);

    void Record(DynamicAnalysis& da) override;
		
};

/// <summary>
/// 時刻歴応答解析クラス
/// </summary>
class DynamicAnalysis : public FEDeformCase {
private:
#ifdef EIGEN_USE_MKL_ALL
    Eigen::PardisoLLT<Eigen::SparseMatrix<double>> solver;
#else
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
#endif
	Eigen::SparseMatrix<double> mass_mat;
	Eigen::SparseMatrix<double> stiffness_mat;
	Eigen::SparseMatrix<double> damping_mat;
	Eigen::SparseMatrix<double> compute_mat;
	// Eigen::SparseMatrix<double> m_sh, k_sh, d_sh;
	// Eigen::SparseMatrix<double> m_sha, m_shb, m_shc, m_shd;
	std::vector<int> free_indices, fixed_indices;
	//std::vector<int> shrink_indices, other_indices;

    Eigen::VectorXd current_disp, current_vel, current_accel;

    friend class FEDynamicDampInitializer;
    friend class FEDynamicStiffDampInitializer;
	friend class DAEnergyRecorder;

public:
	//std::shared_ptr<FEModel> model;
	DynamicAccelLoad accel_load;
	FEDynamicDampInitializer* damp_initializer = nullptr;

	int current_step = 0;
	//double damping_rate = 0.05; // 減衰比

    std::list<std::shared_ptr<DASampler>> samplers;
	DAEnergyRecorder energy_recorder;

	//const double tmp_w1 = 720;
	//const double tmp_w1 = 4200;
	//const double tmp_w1 = PI * 4.0;
	double beta = 0.25; // 平均加速度法
	//double beta = 1.0/6.0; // 線形加速度法(発散しがち)
	
    //DynamicAnalysis(FEModel* model, DynamicAccelLoad accel_load)
    //    : model(std::make_shared<FEModel>(model)), accel_load(accel_load) {
    //}

    DynamicAnalysis(std::shared_ptr<FEModel> model, DynamicAccelLoad accel_load, FEDynamicDampInitializer* damp = nullptr);
	//	: model(model), accel_load(accel_load) {
	//};

    bool Initialize();

	void Clear() {
		current_step = 0;
		//current_disp.clear();
		//current_vel.clear();
		//current_accel.clear();
		mass_mat.resize(0, 0);
		stiffness_mat.resize(0, 0);
		damping_mat.resize(0, 0);
		compute_mat.resize(0, 0);
	}

	// Newmarkのβ法による動的解析
    void ComputeStep();
	void ComputeSteps(int steps);

    std::vector<Displacement> GetDisplacements();
	std::vector<Displacement> GetVelocities();
	std::vector<Displacement> GetAccelerations();

    std::vector<BarElementBase*> GetBarElements() {
        std::vector<BarElementBase*> bar_elements;
        for (const auto& elem : model->Elements) {
            if (elem->Type() == ElementType::Beam || elem->Type() == ElementType::Truss) {
                if (auto* bar_elem = dynamic_cast<BarElementBase*>(elem.get())) {
                    bar_elements.push_back(bar_elem);
                }
            }
        }
        return bar_elements;
    }

    // FEDeformCase を介して継承されました
    BeamStressData GetBeamStress(int eid, double p) override;
    PlateStressData GetPlateStressData(int eid, double xi, double eta) override;
    Displacement GetBeamDisplace(int eid, double p) override;
};

/// <summary>
/// 時刻歴応答解析における比例減衰マトリクスを初期化する抽象基底クラス
/// </summary>
class FEDynamicDampInitializer {
private:
	DynamicAnalysis* analysis;
public:
    virtual bool Initialize(DynamicAnalysis* analysis) = 0;
};


/// <summary>
/// 時刻歴応答解析において比例減衰マトリクスを剛性比例として初期化するクラス
/// </summary>
class FEDynamicStiffDampInitializer : public FEDynamicDampInitializer {
    
public:
	double natural_angle_velocity = 0.0; // 自然角速度
	double damp_rate = 0.05; // 減衰比

	FEDynamicStiffDampInitializer(double damp_rate = 0.05)
		: damp_rate(damp_rate) { }
    
    bool Initialize(DynamicAnalysis* analysis) override;
};

enum ResponseSpectrumMethodType
{
	ABS, SRSS, CQC
};


class ResponseSpectrumMethod : public FEDeformCase {
private:
    std::vector<Displacement> calculate_displacementsCQC();
    std::vector<Displacement> calculate_displacementsSRSS();
    std::vector<Displacement> calculate_displacementsABS();
    std::vector<Displacement> calculate_displacements();
    std::vector<Displacement> displacements; // 解析結果の変位ベクトル

public:
    
    double damping_rate = 0.05; // 減衰比(CQC法の場合のみ計算に影響)
	Vector Direction = Vector(1.0, 1.0, 1.0); // 応答スペクトルの方向

    //std::shared_ptr<FEModel> model;
    FEVibrateResult VibrateResult;
    IResponseSpectrum* SpectrumFunction;
	ResponseSpectrumMethodType MethodType = ResponseSpectrumMethodType::ABS;

    ResponseSpectrumMethod() {}
    ResponseSpectrumMethod(std::shared_ptr<FEModel> model, FEVibrateResult vibrate_result, Vector direction,
        IResponseSpectrum* spectrum_function, ResponseSpectrumMethodType type);
		

    std::vector<Displacement> GetDisplacements() override;
    //std::vector<Displacement> GetVelocities();
    std::vector<Displacement> GetAccelerations();

    // FEDeformCase を介して継承されました
    BeamStressData GetBeamStress(int eid, double p) override;
    PlateStressData GetPlateStressData(int eid, double xi, double eta) override;
    Displacement GetBeamDisplace(int eid, double p) override;
};

#endif
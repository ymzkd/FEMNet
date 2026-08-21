// fem_common.i - Common SWIG definitions (language independent)

%include <std_vector.i>
%include <std_map.i>
%include <std_array.i>
%include <std_list.i>
%include <std_shared_ptr.i>
%include <std_string.i>
// 出力引数(double& 等)を各言語ネイティブの出力に変換するための %apply 用
%include <typemaps.i>

// Note: Default constructors have been added to BeamPolyLoad, Node, and Section
// so they can now be used in STL containers

// Shared pointer declarations for base classes
%shared_ptr(FEModel);
%shared_ptr(RigidLinks);
// Load Pointer
%shared_ptr(LoadBase);
%shared_ptr(PlateLoad);
%shared_ptr(BeamPolyLoad);
%shared_ptr(AxialPolyLoad);
%shared_ptr(BeamLoadBase);
// NodeLoad inherits from LoadBase, so it must also be marked as shared_ptr
// to avoid SWIG generating inconsistent destructor code (smartarg1 error)
%shared_ptr(NodeLoad);
%shared_ptr(InertialForce);
%shared_ptr(NodeBodyForce);

%{
    #include <Eigen/Dense>
    #include <Eigen/Sparse>

    #ifdef EIGEN_USE_MKL_ALL
    #include <Eigen/PardisoSupport>
    #endif

    #include "Components.h"
    #include "Elements/Elements.h"
    #include "LoadComponent.h"
    #include "RigidLink.h"
    #include "Model.h"
    #include "SeismicModule.h"
    #include "LoadIO.h"
%}

// Ignore non-const accessors that return references (use const versions instead)
%ignore NodeLoadData::Px();
%ignore NodeLoadData::Py();
%ignore NodeLoadData::Pz();
%ignore NodeLoadData::Mx();
%ignore NodeLoadData::My();
%ignore NodeLoadData::Mz();

// DOFFlags も同様に非const版(bool&)を隠し、const版(値返し)のみを公開する。
// 非const版が優先されると bool* の不透明ポインタになってしまう。
// 値の書き換えは Set(index, value) を使う。
%ignore DOFFlags::Ux();
%ignore DOFFlags::Uy();
%ignore DOFFlags::Uz();
%ignore DOFFlags::Rx();
%ignore DOFFlags::Ry();
%ignore DOFFlags::Rz();

// 生のC配列メンバは不透明ポインタ(double*, bool*)にしかならないため公開しない。
// 値はアクセサ(Px()/Dx() 等)経由で取得する。
%ignore DOFFlags::flags;
%ignore Displacement::displace;
%ignore NodeLoadData::loads;

// NodeLoad のアクセサは double&/double* を返すため不透明ポインタになる。
// 値返しの同名メソッドを %extend で差し替える(ファイル末尾の %extend NodeLoad 参照)。
%ignore NodeLoad::loads;
%ignore NodeLoad::Px();
%ignore NodeLoad::Py();
%ignore NodeLoad::Pz();
%ignore NodeLoad::Mx();
%ignore NodeLoad::My();
%ignore NodeLoad::Mz();

// Ignore Eigen types that cannot be wrapped
%ignore Eigen::SparseMatrix;
%ignore Eigen::MatrixXd;
%ignore Eigen::VectorXd;
%ignore Eigen::Vector3d;
%ignore trans_matrix3(const Point p0, const Point p1, const double beta);
%ignore trans_matrix3(const Plane plane);
%ignore Eigen::Matrix3d;
%ignore Displacement::translate(Eigen::Matrix3d transmat);
%ignore Vector::toEigen;

// Ignore RigidLink methods that return Eigen types
%ignore RigidLink::TransformationMatrix;
%ignore RigidLinks::TransformationMatrix;

// FEModelの組立サービス(Eigen型を返すためSWIG非公開)。
// 解析はOperator側(FELinearStaticOp, FEVibrationAnalysis等)から利用する。
%ignore FEModel::AssembleStiffnessMatrix;
%ignore FEModel::AssembleMassMatrix;
%ignore FEModel::AssembleGeometricStiffnessMatrix;
%ignore FEModel::AssembleLoadVector;

// Ignore pure virtual methods that use Eigen types
%ignore ElementBase::geometric_local_stiffness_matrix;
%ignore ElementBase::InertialForceToNodeLoadData;
%ignore ElementBase::NodeLumpedMass;
%ignore ElementBase::NodeConsistentMass;
%ignore ElementBase::StiffnessMatrix;
%ignore ElementBase::AssembleStiffMatrix;
%ignore ElementBase::AssembleGeometricStiffMatrix;
%ignore ElementBase::AssembleMassMatrix;

// ===================================================================
// STEP 1: Include Components.h first (defines basic types)
// ===================================================================
%include "Components.h"

// ===================================================================
// STEP 2: Include Elements module (depends on Components.h)
// ===================================================================
%include "Elements_common.i"

// ===================================================================
// STEP 3: Include LoadComponent.h (depends on Components.h and Elements)
// ===================================================================
%include "LoadComponent.h"

// ===================================================================
// STEP 4: STL templates (AFTER all classes are defined)
// ===================================================================
namespace std {
    // Basic types
    %template(VectorInt) std::vector<int>;
    %template(VectorDouble) std::vector<double>;

    // Support::isdof_fixed() の戻り値(未定義だと bool[6] の不透明ポインタになる)
    %template(ArrayBool6) std::array<bool, 6>;

    // Vector(3次元ベクトル)のコンテナ。PlateLoad::load_vecs と
    // PlaneElementBase::AreaForceToNodeLoadData() で使用する。
    %template(VectorVector) std::vector<Vector>;

    // Classes with default constructors
    %template(VectorDisp) std::vector<Displacement>;
    %template(VectorMode) std::vector<std::vector<Displacement>>;
    %template(VectorMaterial) std::vector<Material>;

    // Node and Section vectors (default constructors added)
    %template(VectorNode) std::vector<Node>;
    %template(VectorSection) std::vector<Section>;

    // Load-related vectors
    %template(VectorLoad) std::vector<std::shared_ptr<LoadBase>>;
    %template(VectorNodeLoad) std::vector<NodeLoad>;
    %template(VectorNodeBodyForce) std::vector<NodeBodyForce>;
    %template(VectorNodeLoadData) std::vector<NodeLoadData>;

    // StressData vectors
    %template(VectorBeamStressData) std::vector<BeamStressData>;
    %template(VectorPlateStressData) std::vector<PlateStressData>;

    // BeamPolyLoad vectors and lists (default constructors added)
    %template(VectorBeamPolyLoad) std::vector<BeamPolyLoad>;
    %template(ListBeamPolyLoad) std::list<BeamPolyLoad>;
    %template(VectorBeamPolyLoadPtr) std::vector<std::shared_ptr<BeamPolyLoad>>;

    // RigidLink vectors
    %template(VectorRigidLink) std::vector<RigidLink>;
}

// ===================================================================
// STEP 5: Include remaining modules
// ===================================================================
%include "Operator_common.i"
%include "SeismicModule_common.i"

// Include Model.h and SeismicModule.h
%include "RigidLink.h"
%include "Model.h"
%include "SeismicModule.h"
%include "LoadIO.h"

// ===================================================================
// Class extensions (AFTER all classes are fully defined)
// ===================================================================

// NodeLoad: double& を返す元のアクセサ(%ignore 済み)の代わりに
// 値返しのアクセサを提供する。実体は data(NodeLoadData)の const アクセサ。
%extend NodeLoad {
    double Px() const { return $self->data.Px(); }
    double Py() const { return $self->data.Py(); }
    double Pz() const { return $self->data.Pz(); }
    double Mx() const { return $self->data.Mx(); }
    double My() const { return $self->data.My(); }
    double Mz() const { return $self->data.Mz(); }
};

// Material extension (language independent)
%extend Material {
    std::string to_string() {
       std::ostringstream oss;
       oss << "Material value: " << *$self;
       return oss.str();
   }
};

// FEModel extension - helper methods for adding components
%extend FEModel {
    // Add a node with id and coordinates
    Node& AddNode(int id, double x, double y, double z) {
        $self->Nodes.push_back(Node(id, x, y, z));
        return $self->Nodes.back();
    }

    // Add a node with coordinates only (id will be -1)
    Node& AddNodeXYZ(double x, double y, double z) {
        $self->Nodes.push_back(Node(x, y, z));
        return $self->Nodes.back();
    }

    // Get node by index
    Node& GetNode(int index) {
        return $self->Nodes[index];
    }

    // Add a section
    Section& AddSection(double A, double Iy, double Iz, double K) {
        $self->Sections.push_back(Section(A, Iy, Iz, K));
        return $self->Sections.back();
    }

    // Get section by index
    Section& GetSection(int index) {
        return $self->Sections[index];
    }

    // Add a material with Young's modulus and Poisson's ratio
    Material& AddMaterial(double young, double poisson) {
        $self->Materials.push_back(Material(young, poisson));
        return $self->Materials.back();
    }

    // Add a material with density
    Material& AddMaterialWithDensity(double young, double poisson, double dense) {
        $self->Materials.push_back(Material(young, poisson, dense));
        return $self->Materials.back();
    }

    // Get material by index
    Material& GetMaterial(int index) {
        return $self->Materials[index];
    }
};

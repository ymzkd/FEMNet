// Operator_common.i - Common Analysis/Operator definitions (language independent)

// Analysis Pointer definitions
%shared_ptr(DASampler);
%shared_ptr(DASampler_MaxDisplacement);
%shared_ptr(DASampler_MaxVelocity);
%shared_ptr(DASampler_MaxAcceleration);
%shared_ptr(DASampler_MaxDispDirection);
%shared_ptr(DASampler_MaxVelocityDirection);
%shared_ptr(DASampler_MaxAccelDirection);
// レコーダ(DynamicAnalysis が複数保持するため shared_ptr で公開)
%shared_ptr(DARecorder);
%shared_ptr(DARecorder_KineticEnergy);
%shared_ptr(DARecorder_PotentialEnergy);
%shared_ptr(DARecorder_DampingEnergy);
%shared_ptr(DARecorder_InputEnergy);
%shared_ptr(FEDeformOperator);
%shared_ptr(FEModeOperator);
%shared_ptr(FELinearStaticOp);
%shared_ptr(LinearStaticCombinationOperator);
%shared_ptr(DynamicAnalysis);
// 減衰初期化子(DynamicAnalysis が共同所有するため shared_ptr で公開)
%shared_ptr(FEDynamicDampInitializer);
%shared_ptr(FEDynamicStiffDampInitializer);
%shared_ptr(FEDynamicMassDampInitializer);
%shared_ptr(FEDynamicRayleighDampInitializer);
// 時刻歴荷重(慣性力以外にも対応する DynamicLoad 階層)
%shared_ptr(DynamicLoad);
%shared_ptr(SeismicAccelLoad);
%shared_ptr(NodalDynamicLoad);
%shared_ptr(FEBucklingAnalysis);
%shared_ptr(FEVibrationAnalysis);
%shared_ptr(ResponseSpectrumMethod);

// Director feature for polymorphic classes
%feature("director") IResponseSpectrum;
%feature("director") DASampler;
%feature("director") DARecorder;
// 時刻歴荷重は .NET / Python 側で任意の時間関数として実装できるようにする
%feature("director") DynamicLoad;
// time_series() は内部データへの参照を返すため、.NET 側でオーバーライドすると
// 返却したオブジェクトが GC で回収され得る。この2つは director から除外し、
// .NET 側の派生クラスは load_vector / reference_value 等で実装する。
%feature("nodirector") DynamicLoad::has_time_series;
%feature("nodirector") DynamicLoad::time_series;

%{
    #include "FEAnalysis.h"
    #include "FELinearStaticOp.h"
    #include "FEDynamic.h"
    #include "FEBucklingAnalysis.h"
    #include "FEVibrationAnalysis.h"
    #include "ResponseSpectrumMethod.h"
%}

// STL templates for Analysis
namespace std {
    %template(VectorDASampler) std::vector<std::shared_ptr<DASampler>>;
    %template(VectorDARecorder) std::vector<std::shared_ptr<DARecorder>>;
    %template(VectorDynamicLoad) std::vector<std::shared_ptr<DynamicLoad>>;
    %template(LinearStaticDeformFactorVector) std::vector<LinearStaticDeformFactor>;
}

%include "FEAnalysis.h"
%include "FELinearStaticOp.h"
%include "FEDynamic.h"
%include "FEBucklingAnalysis.h"
%include "FEVibrationAnalysis.h"
%include "ResponseSpectrumMethod.h"

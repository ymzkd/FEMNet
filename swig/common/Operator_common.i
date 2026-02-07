// Operator_common.i - Common Analysis/Operator definitions (language independent)

// Note: DynamicAccelLoad now has a default constructor, so accel_load can be accessed directly

// Analysis Pointer definitions
%shared_ptr(DASampler);
%shared_ptr(DASampler_MaxDisplacement);
%shared_ptr(FEDeformOperator);
%shared_ptr(FEModeOperator);
%shared_ptr(FEStaticResult);
%shared_ptr(FELinearStaticOp);
%shared_ptr(LinearStaticCombinationOperator);
%shared_ptr(DynamicAnalysis);
%shared_ptr(FEBucklingAnalysis);
%shared_ptr(FEVibrateResult);
%shared_ptr(ResponseSpectrumMethod);

// Director feature for polymorphic classes
%feature("director") IResponseSpectrum;
%feature("director") DASampler;

%{
    #include "FEAnalysis.h"
    #include "FELinearStaticOp.h"
    #include "FEDynamic.h"
    #include "FEBucklingAnalysis.h"
    #include "FEVibrateResult.h"
    #include "ResponseSpectrumMethod.h"
%}

// STL templates for Analysis
namespace std {
    %template(VectorDASampler) std::vector<std::shared_ptr<DASampler>>;
    %template(LinearStaticDeformFactorVector) std::vector<LinearStaticDeformFactor>;
}

%include "FEAnalysis.h"
%include "FELinearStaticOp.h"
%include "FEDynamic.h"
%include "FEBucklingAnalysis.h"
%include "FEVibrateResult.h"
%include "ResponseSpectrumMethod.h"

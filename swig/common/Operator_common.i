// Operator_common.i - Common Analysis/Operator definitions (language independent)

// Ignore DynamicAccelLoad member to avoid default constructor issues
// (DynamicAccelLoad has no default constructor, and SWIG generates code that requires it)
%ignore DynamicAnalysis::accel_load;

// Provide access to accel_load through a method returning a pointer
%extend DynamicAnalysis {
    const DynamicAccelLoad* get_accel_load() const {
        return &($self->accel_load);
    }
}

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
}

%include "FEAnalysis.h"
%include "FELinearStaticOp.h"
%include "FEDynamic.h"
%include "FEBucklingAnalysis.h"
%include "FEVibrateResult.h"
%include "ResponseSpectrumMethod.h"

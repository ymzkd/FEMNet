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

// Director機能を有効化
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

// C#側で定義したIResponseSpectrumをGCに解放されないように
// 保持するためのメソッドを追加
%typemap(cscode) ResponseSpectrumMethod %{
    private IResponseSpectrum __sptectrumRefs;

    public void SetSpectrum(IResponseSpectrum spectrum) {
        __sptectrumRefs = spectrum;
        this.SpectrumFunction = spectrum;
    }
%}

// STL templates for Analysis
namespace std {
    %template(VectorDASampler) std::vector<std::shared_ptr<DASampler>>;
}

// FEDeformOperatorクラスにC#プロパティを追加
%typemap(cscode) FEDeformOperator %{
    /// <summary>
    /// Operation Name
    /// </summary>
    public virtual string OperationName { get; set; } = "";

    /// <summary>
    /// Operation Description 
    /// step summary, time, computation algorithm, etc.
    /// </summary>
    public virtual string OperationDescription { get; set; } = "";
%}


%include "FEAnalysis.h"
%include "FELinearStaticOp.h"
%include "FEDynamic.h"
%include "FEBucklingAnalysis.h"
%include "FEVibrateResult.h"
%include "ResponseSpectrumMethod.h"
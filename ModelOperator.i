// ModelOperator.i - Analysis module definitions for SWIG
// %include <std_shared_ptr.i>
// %include <std_vector.i>

// Analysis Pointer definitions
// %shared_ptr(FEModel);
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
// %feature("director") IResponseSpectrum;
// %feature("director") DASampler;

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

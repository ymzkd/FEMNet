%module (directors="1") FEMNet
// %include "SeismicModule.i"
%include <std_vector.i>
%include <std_map.i>
%include <std_array.i>
%include <std_list.i>
%include <std_shared_ptr.i>
%include <std_string.i>


%shared_ptr(FEModel);
// Load Pointer
%shared_ptr(LoadBase);
%shared_ptr(PlateLoad);
%shared_ptr(BeamPolyLoad);
%shared_ptr(AxialPolyLoad);
%shared_ptr(BeamLoadBase);
%shared_ptr(NodeLoad);
%shared_ptr(InertialForce);
%shared_ptr(NodeBodyForce);


// Director機能を有効化
%feature("director") IResponseSpectrum;
%feature("director") DASampler;

%{
    #include <Eigen/Dense>
    #include <Eigen/Sparse>
    
    #ifdef EIGEN_USE_MKL_ALL
    #include <Eigen/PardisoSupport>
    #endif
    
    #include "Elements/Elements.h"
    #include "Components.h"
    #include "Model.h"
    #include "FEAnalysis.h"
    #include "FELinearStaticOp.h"
    #include "FEDynamic.h"
    #include "FEBucklingAnalysis.h"
    #include "FEVibrateResult.h"
    #include "ResponseSpectrumMethod.h"
    #include "SeismicModule.h"
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

%include <std_string.i>
// %include <carrays.i>
// %array_class(double, doubleArray);
// %array_functions(double, doubleArr);

%ignore Eigen::SparseMatrix;
%ignore Eigen::MatrixXd;
%ignore extractSubMatrix();
%ignore trans_matrix3(const Point p0, const Point p1, const double beta);



%ignore Eigen::Matrix3d;
%ignore Displacement::translate(Eigen::Matrix3d transmat);

namespace std {
    // List
    // %template(ListLoad) std::list<Load>;
    %template(ListBeamPolyLoad) std::list<BeamPolyLoad>;
    
    // Vector
    %template(VectorNode) std::vector<Node>;
    %template(VectorDisp) std::vector<Displacement>;
    %template(VectorMode) std::vector<std::vector<Displacement>>;
    %template(VectorMaterial) std::vector<Material>;
	%template(VectorSection) std::vector<Section>;
	%template(VectorInt) std::vector<int>;
	%template(VectorDouble) std::vector<double>;

	%template(VectorBeamPolyLoad) std::vector<BeamPolyLoad>;
	%template(VectorLoad) std::vector<std::shared_ptr<LoadBase>>;
	%template(VectorNodeLoad) std::vector<NodeLoad>;
	%template(VectorNodeBodyForce) std::vector<NodeBodyForce>;
    %template(VectorNodeLoadData) std::vector<NodeLoadData>;
}

// カスタムメソッドを追加
%extend Material {
    std::string to_string() {
       std::ostringstream oss;
       oss << "Material value: " << *$self;
       return oss.str();
   }
};



// C#のカスタムコードを追加
%typemap(cscode) Material %{
    // ToStringメソッドをオーバーライド
    public override string ToString()
    {
        return to_string();
    }
%}

%ignore NodeLoadData::Px();
%ignore NodeLoadData::Py();
%ignore NodeLoadData::Pz();
%ignore NodeLoadData::Mx();
%ignore NodeLoadData::My();
%ignore NodeLoadData::Mz();

// %typemap(csout) double& %{
//     $result = $imcall;
// %}

// %ignore ElementHandle(ElementBase data);
// %ignore ElementHandle::ElementHandle(ElementBase data);
// %typemap(csin) ElementBaseSetPtr data "getCPtrAndAddReference($csinput)"
// %typemap(cscode) ElementHandle %{
//     private ElementBase elementReference;
//     private global::System.Runtime.InteropServices.HandleRef getCPtrAndAddReference(ElementBase element) {
//         elementReference = element;
//         return ElementBase.getCPtr(element);
//     }
//     public void SetElement(ElementBase element){
//         getCPtrAndAddReference(element);
//         set_element(element);
//     }
// %}

// Include module-specific SWIG files after C++ headers
%include "Elements/Elements.i"
%include "ModelOperator.i"


%include "Elements/Elements.h"
%include "Components.h"
%include "Model.h"
%include "FEAnalysis.h"
%include "FELinearStaticOp.h"
%include "FEDynamic.h"
%include "FEBucklingAnalysis.h"
%include "FEVibrateResult.h"
%include "ResponseSpectrumMethod.h"
%include "SeismicModule.h"

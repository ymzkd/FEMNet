%module (directors="1") FEMNet

%include <std_vector.i>
%include <std_map.i>
%include <std_array.i>
%include <std_list.i>
%include <std_shared_ptr.i>

%shared_ptr(FEModel);
// %shared_ptr(SSModel);

// Element Pointer
%shared_ptr(ElementBase);
%shared_ptr(BarElementBase);
%shared_ptr(TrussElement);
%shared_ptr(BeamElement);
%shared_ptr(ComplexBeamElement);
%shared_ptr(PlaneElementBase);
%shared_ptr(TriPlateElement);
%shared_ptr(QuadPlateElement);
%shared_ptr(TriPlaneElement);
%shared_ptr(QuadPlaneElement);

// Load Pointer
%shared_ptr(LoadBase);
%shared_ptr(PlateLoad);
%shared_ptr(BeamPolyLoad);
%shared_ptr(AxialPolyLoad);
%shared_ptr(BeamLoadBase);
%shared_ptr(NodeLoad);
%shared_ptr(InertialForce);

// Director機能を有効化
%feature("director") IResponseSpectrum;

%{
    #include <Eigen/Dense>
    #include <Eigen/Sparse>
    
    #ifdef EIGEN_USE_MKL_ALL
    #include <Eigen/PardisoSupport>
    #endif
    
    #include <Eigen/PardisoSupport>

    #include "Element.h"
    #include "Components.h"
    #include "Model.h"
%}

%include <std_string.i>
// %include <carrays.i>
// %array_class(double, doubleArray);
// %array_functions(double, doubleArr);

%ignore Eigen::SparseMatrix;
%ignore Eigen::MatrixXd;
%ignore extractSubMatrix();
%ignore trans_matrix3(const Point p0, const Point p1, const double beta);


%ignore BeamElement::StiffnessMatrix();
%ignore BeamElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore BeamElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore BeamElement::NodeLumpedMass();

%ignore TriPlaneElement::StiffnessMatrix();
%ignore TriPlaneElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TriPlaneElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TriPlaneElement::NodeLumpedMass();

%ignore TriPlateElement::StiffnessMatrix();
%ignore TriPlateElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TriPlateElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TriPlateElement::NodeLumpedMass();

%ignore QuadPlaneElement::StiffnessMatrix();
%ignore QuadPlaneElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore QuadPlaneElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore QuadPlaneElement::NodeLumpedMass();

%ignore QuadPlateElement::StiffnessMatrix();
%ignore QuadPlateElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore QuadPlateElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore QuadPlateElement::NodeLumpedMass();

%ignore ElementBase::StiffnessMatrix();
%ignore ElementBase::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore ElementBase::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore ElementBase::NodeLumpedMass();  

%ignore BarElementBase::StiffnessMatrix();
%ignore BarElementBase::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore BarElementBase::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore BarElementBase::NodeLumpedMass();

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
    %template(VectorElement) std::vector<std::shared_ptr<ElementBase>>;
    %template(VectorElem) std::vector<ElementBase*>;
    %template(VectorBars) std::vector<BarElementBase*>;
    %template(VectorBeams) std::vector<BeamElement*>;
    %template(VectorMaterial) std::vector<Material>;
	%template(VectorSection) std::vector<Section>;
	%template(VectorInt) std::vector<int>;
	%template(VectorDouble) std::vector<double>;
    

	%template(VectorBeamPolyLoad) std::vector<BeamPolyLoad>;
	%template(VectorLoad) std::vector<std::shared_ptr<LoadBase>>;
	%template(VectorNodeLoad) std::vector<NodeLoad>;
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

%extend BarElementBase {
    Node* getNodes(int index) {
        return $self->Nodes[index];
    }
}

%extend BeamElement {
    Node* getNodes(int index) {
        return $self->Nodes[index];
    }
}

%extend TrussElement {
    Node* getNodes(int index) {
        return $self->Nodes[index];
    }
}

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

#include "Element.h"
#include "Components.h"
#include "Model.h"


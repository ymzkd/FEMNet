%module FEMNet

%include <std_vector.i>
%include <std_map.i>
%include <std_array.i>
%include <std_list.i>

%{
    #include <Eigen/Dense>
    #include <Eigen/Sparse>
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

%ignore BeamElement::StiffnessMatrix();
%ignore BeamElement::AssembleStiffMatrix();

%ignore TriPlaneElement::StiffnessMatrix();
%ignore TriPlaneElement::AssembleStiffMatrix();

%ignore TriPlateElement::StiffnessMatrix();
%ignore TriPlateElement::AssembleStiffMatrix();

%ignore QuadPlaneElement::StiffnessMatrix();
%ignore QuadPlaneElement::AssembleStiffMatrix();

%ignore QuadPlateElement::StiffnessMatrix();
%ignore QuadPlateElement::AssembleStiffMatrix();

%ignore ElementBase::StiffnessMatrix();
%ignore ElementBase::AssembleStiffMatrix();

%ignore Eigen::Matrix3d;
%ignore Displacement::translate();

namespace std {
    // List
    %template(ListLoad) std::list<Load>;
    
    // Vector
    %template(VectorNode) std::vector<Node>;
    %template(VectorDisp) std::vector<Displacement>;
    %template(VectorElement) std::vector<ElementHandle>;
    %template(VectorElem) std::vector<ElementBase*>;
    %template(VectorBeams) std::vector<BeamElement*>;
    %template(VectorMaterial) std::vector<Material>;
	%template(VectorSection) std::vector<Section>;
	%template(VectorInt) std::vector<int>;
}

// カスタムメソッドを追加
%extend Material {
    std::string to_string() {
       std::ostringstream oss;
       oss << "Material value: " << *$self;
       return oss.str();
   }
};

%extend BeamElement {
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

%ignore ElementHandle(ElementBase data);
%ignore ElementHandle::ElementHandle(ElementBase data);
%typemap(csin) ElementBaseSetPtr data "getCPtrAndAddReference($csinput)"
%typemap(cscode) ElementHandle %{
    private ElementBase elementReference;
    private global::System.Runtime.InteropServices.HandleRef getCPtrAndAddReference(ElementBase element) {
        elementReference = element;
        return ElementBase.getCPtr(element);
    }
    public void SetElement(ElementBase element){
        getCPtrAndAddReference(element);
        set_element(element);
    }
%}

#include "Element.h"
#include "Components.h"
#include "Model.h"


// Elements.i - Elements module definitions for SWIG

// Element Pointer definitions
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

%{
    #include "Elements/Elements.h"
%}

%ignore ShapeFunction4(double xi, double eta);
%ignore ShapeFunctionTriangle6(double xi, double eta);
%ignore ShapeFunctionSerendipity8(double xi, double eta);

// Ignore Element methods that use Eigen matrices
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

// STL templates for Elements
namespace std {
    %template(VectorElement) std::vector<std::shared_ptr<ElementBase>>;
    %template(VectorElem) std::vector<ElementBase*>;
    %template(VectorBars) std::vector<BarElementBase*>;
    %template(VectorBeams) std::vector<BeamElement*>;
}

// Element extensions for accessing nodes
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

%include "Elements/Elements.h"
// Elements_common.i - Common Elements module definitions (language independent)

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

// Ignore shape functions
%ignore ShapeFunction4(double xi, double eta);
%ignore ShapeFunctionTriangle6(double xi, double eta);
%ignore ShapeFunctionSerendipity8(double xi, double eta);

// Ignore Element methods that use Eigen matrices
%ignore BeamElement::StiffnessMatrix();
%ignore BeamElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore BeamElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore BeamElement::NodeLumpedMass();

%ignore TrussElement::StiffnessMatrix();
%ignore TrussElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TrussElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TrussElement::NodeLumpedMass();

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
%ignore BarElementBase::stress;

%ignore ComplexBeamElement::StiffnessMatrix();
%ignore ComplexBeamElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore ComplexBeamElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore ComplexBeamElement::NodeLumpedMass();

%ignore PlaneElementBase::StiffnessMatrix();
%ignore PlaneElementBase::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore PlaneElementBase::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore PlaneElementBase::NodeLumpedMass();

// Explicitly include individual element headers for proper SWIG parsing
// IMPORTANT: These must come AFTER %ignore directives
%include "Elements/ElementBase.h"
%include "Elements/BarElement.h"
%include "Elements/TrussElement.h"
%include "Elements/BeamElement.h"
%include "Elements/ComplexBeamElement.h"
%include "Elements/PlaneElement.h"
%include "Elements/TriPlaneElement.h"
%include "Elements/QuadPlaneElement.h"
%include "Elements/TriPlateElement.h"
%include "Elements/QuadPlateElement.h"

// STL templates for Elements
namespace std {
    %template(VectorElement) std::vector<std::shared_ptr<ElementBase>>;
    %template(VectorElem) std::vector<ElementBase*>;
    %template(VectorBars) std::vector<BarElementBase*>;
    %template(VectorBeams) std::vector<BeamElement*>;
}

// Element extensions for accessing nodes (after includes)
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

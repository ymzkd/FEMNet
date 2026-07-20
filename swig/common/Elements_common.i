// Elements_common.i - Common Elements module definitions (language independent)

// Element Pointer definitions
%shared_ptr(ElementBase);
%shared_ptr(BarElementBase);
%shared_ptr(TrussElement);
%shared_ptr(TensionTrussElement);
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
%ignore BeamElement::AssembleGeometricStiffMatrix;
%ignore BeamElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore BeamElement::NodeLumpedMass();
%ignore BeamElement::GetStiffnessTriplets;
%ignore BeamElement::GetGeometricStiffnessTriplets;

%ignore TrussElement::StiffnessMatrix();
%ignore TrussElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TrussElement::AssembleGeometricStiffMatrix;
%ignore TrussElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TrussElement::NodeLumpedMass();
%ignore TrussElement::GetStiffnessTriplets;
%ignore TrussElement::GetGeometricStiffnessTriplets;

// IStateDependentElement は C++ 側の多重継承用の純粋インターフェース。
// C# は単一継承のため SWIG で公開しない (TensionTrussElement に直接メソッドを再露出する)
%ignore IStateDependentElement;
%ignore TensionTrussElement::GetTangentStiffnessTriplets;
%ignore TensionTrussElement::TangentStiffnessMatrix();

%ignore TriPlaneElement::StiffnessMatrix();
%ignore TriPlaneElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TriPlaneElement::AssembleGeometricStiffMatrix;
%ignore TriPlaneElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TriPlaneElement::NodeLumpedMass();
%ignore TriPlaneElement::GetStiffnessTriplets;
%ignore TriPlaneElement::GetGeometricStiffnessTriplets;

%ignore TriPlateElement::StiffnessMatrix();
%ignore TriPlateElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TriPlateElement::AssembleGeometricStiffMatrix;
%ignore TriPlateElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore TriPlateElement::NodeLumpedMass();
%ignore TriPlateElement::GetStiffnessTriplets;
%ignore TriPlateElement::GetGeometricStiffnessTriplets;

%ignore QuadPlaneElement::StiffnessMatrix();
%ignore QuadPlaneElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore QuadPlaneElement::AssembleGeometricStiffMatrix;
%ignore QuadPlaneElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore QuadPlaneElement::NodeLumpedMass();
%ignore QuadPlaneElement::GetStiffnessTriplets;
%ignore QuadPlaneElement::GetGeometricStiffnessTriplets;

%ignore QuadPlateElement::StiffnessMatrix();
%ignore QuadPlateElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore QuadPlateElement::AssembleGeometricStiffMatrix;
%ignore QuadPlateElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore QuadPlateElement::NodeLumpedMass();
%ignore QuadPlateElement::GetStiffnessTriplets;
%ignore QuadPlateElement::GetGeometricStiffnessTriplets;

%ignore ElementBase::StiffnessMatrix();
%ignore ElementBase::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore ElementBase::AssembleGeometricStiffMatrix(Eigen::SparseMatrix<double>& mat, const std::vector<Displacement>& disp);
%ignore ElementBase::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore ElementBase::NodeLumpedMass();
%ignore ElementBase::GetStiffnessTriplets;
%ignore ElementBase::GetGeometricStiffnessTriplets;

%ignore BarElementBase::StiffnessMatrix();
%ignore BarElementBase::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore BarElementBase::AssembleGeometricStiffMatrix;
%ignore BarElementBase::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore BarElementBase::NodeLumpedMass();
%ignore BarElementBase::stress;
%ignore BarElementBase::GetStiffnessTriplets;
%ignore BarElementBase::GetGeometricStiffnessTriplets;

%ignore ComplexBeamElement::StiffnessMatrix();
%ignore ComplexBeamElement::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore ComplexBeamElement::AssembleGeometricStiffMatrix;
%ignore ComplexBeamElement::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore ComplexBeamElement::NodeLumpedMass();
%ignore ComplexBeamElement::GetStiffnessTriplets;
%ignore ComplexBeamElement::GetGeometricStiffnessTriplets;

%ignore PlaneElementBase::StiffnessMatrix();
%ignore PlaneElementBase::AssembleStiffMatrix(Eigen::SparseMatrix<double>& mat);
%ignore PlaneElementBase::AssembleGeometricStiffMatrix;
%ignore PlaneElementBase::AssembleMassMatrix(Eigen::SparseMatrix<double>& mat);
%ignore PlaneElementBase::NodeLumpedMass();
%ignore PlaneElementBase::GetStiffnessTriplets;
%ignore PlaneElementBase::GetGeometricStiffnessTriplets;

// Explicitly include individual element headers for proper SWIG parsing
// IMPORTANT: These must come AFTER %ignore directives
%include "Elements/ElementBase.h"
%include "Elements/BarElement.h"
%include "Elements/TrussElement.h"
%include "Elements/TensionTrussElement.h"
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
    %template(VectorPlanes) std::vector<PlaneElementBase*>;
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

%extend TensionTrussElement {
    Node* getNodes(int index) {
        return $self->Nodes[index];
    }

    // IStateDependentElement::IsActive は SWIG で非公開のため、TensionTrussElement に
    // getter/setter として再露出する (C# からは tt.GetIsActive()/SetIsActive() でアクセス)
    bool GetIsActive() { return $self->IsActive; }
    void SetIsActive(bool val) { $self->IsActive = val; }
}

%extend PlaneElementBase {
    Node* getNodes(int index) {
        return $self->NodesList()[index];
    }
}

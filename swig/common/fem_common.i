// fem_common.i - Common SWIG definitions (language independent)

%include <std_vector.i>
%include <std_map.i>
%include <std_array.i>
%include <std_list.i>
%include <std_shared_ptr.i>
%include <std_string.i>

// Disable default constructor generation for classes without default constructors
%nodefaultctor BeamPolyLoad;
%nodefaultctor DynamicAccelLoad;
%nodefaultctor Node;
%nodefaultctor Section;

// Shared pointer declarations for base classes
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

%{
    #include <Eigen/Dense>
    #include <Eigen/Sparse>

    #ifdef EIGEN_USE_MKL_ALL
    #include <Eigen/PardisoSupport>
    #endif

    #include "Components.h"
    #include "LoadComponent.h"
    #include "Model.h"
    #include "SeismicModule.h"
%}

// Ignore Eigen types that cannot be wrapped
%ignore Eigen::SparseMatrix;
%ignore Eigen::MatrixXd;
%ignore extractSubMatrix();
%ignore trans_matrix3(const Point p0, const Point p1, const double beta);
%ignore Eigen::Matrix3d;
%ignore Displacement::translate(Eigen::Matrix3d transmat);

// Ignore NodeLoadData accessors
%ignore NodeLoadData::Px();
%ignore NodeLoadData::Py();
%ignore NodeLoadData::Pz();
%ignore NodeLoadData::Mx();
%ignore NodeLoadData::My();
%ignore NodeLoadData::Mz();

// STL templates for common types
// Note: Classes without default constructors (BeamPolyLoad, DynamicAccelLoad, Node, Section)
// cannot be used in STL container templates for Python bindings directly.
// Use pointer-based vectors or access through methods instead.
namespace std {
    // Basic types
    %template(VectorInt) std::vector<int>;
    %template(VectorDouble) std::vector<double>;

    // Classes with default constructors
    %template(VectorDisp) std::vector<Displacement>;
    %template(VectorMode) std::vector<std::vector<Displacement>>;
    %template(VectorMaterial) std::vector<Material>;

    // Load-related vectors
    %template(VectorLoad) std::vector<std::shared_ptr<LoadBase>>;
    %template(VectorNodeLoad) std::vector<NodeLoad>;
    %template(VectorNodeBodyForce) std::vector<NodeBodyForce>;
    %template(VectorNodeLoadData) std::vector<NodeLoadData>;

    // For classes without default constructors, use pointer vectors
    %template(VectorBeamPolyLoadPtr) std::vector<std::shared_ptr<BeamPolyLoad>>;
}

// Material extension (language independent)
%extend Material {
    std::string to_string() {
       std::ostringstream oss;
       oss << "Material value: " << *$self;
       return oss.str();
   }
};

// FEModel extension - helper methods for adding components
%extend FEModel {
    // Add a node with id and coordinates
    Node& AddNode(int id, double x, double y, double z) {
        $self->Nodes.push_back(Node(id, x, y, z));
        return $self->Nodes.back();
    }

    // Add a node with coordinates only (id will be -1)
    Node& AddNodeXYZ(double x, double y, double z) {
        $self->Nodes.push_back(Node(x, y, z));
        return $self->Nodes.back();
    }

    // Get node by index
    Node& GetNode(int index) {
        return $self->Nodes[index];
    }

    // Add a section
    Section& AddSection(double A, double Iy, double Iz, double K) {
        $self->Sections.push_back(Section(A, Iy, Iz, K));
        return $self->Sections.back();
    }

    // Get section by index
    Section& GetSection(int index) {
        return $self->Sections[index];
    }

    // Add a material with Young's modulus and Poisson's ratio
    Material& AddMaterial(double young, double poisson) {
        $self->Materials.push_back(Material(young, poisson));
        return $self->Materials.back();
    }

    // Add a material with density
    Material& AddMaterialWithDensity(double young, double poisson, double dense) {
        $self->Materials.push_back(Material(young, poisson, dense));
        return $self->Materials.back();
    }

    // Get material by index
    Material& GetMaterial(int index) {
        return $self->Materials[index];
    }
};

// Include sub-modules
%include "Elements_common.i"
%include "Operator_common.i"
%include "SeismicModule_common.i"

// Include headers
%include "Components.h"
%include "LoadComponent.h"
%include "Model.h"
%include "SeismicModule.h"

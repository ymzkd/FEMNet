// fem_common.i - Common SWIG definitions (language independent)

%include <std_vector.i>
%include <std_map.i>
%include <std_array.i>
%include <std_list.i>
%include <std_shared_ptr.i>
%include <std_string.i>

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
    #include "Model.h"
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
namespace std {
    // List
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

// Material extension (language independent)
%extend Material {
    std::string to_string() {
       std::ostringstream oss;
       oss << "Material value: " << *$self;
       return oss.str();
   }
};

// Include sub-modules
%include "Elements_common.i"
%include "Operator_common.i"
%include "SeismicModule_common.i"

// Include headers
%include "Components.h"
%include "Model.h"
%include "SeismicModule.h"

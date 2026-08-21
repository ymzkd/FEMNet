// fem_python.i - Python-specific SWIG interface file
%module(directors="1") femnet

// Enable automatic docstring generation
%feature("autodoc", "1");

// Exception handling: convert C++ exceptions to Python exceptions
%include <exception.i>
%exception {
    try {
        $action
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "Unknown exception");
    }
}

// Include common definitions
%include "fem_common.i"

// ===================================================================
// Fix: swig::from<NodeLoad> for vector-to-tuple conversion
//
// When %shared_ptr(NodeLoad) is active, SWIG's swig::from<NodeLoad>()
// uses traits_from_ptr which creates a raw NodeLoad* pointer wrapped as
// SWIGTYPE_p_NodeLoad. But the registered Python proxy expects
// SWIGTYPE_p_std__shared_ptrT_NodeLoad_t, so elements become SwigPyObject.
//
// This specialization wraps each NodeLoad in a shared_ptr before conversion,
// matching the %shared_ptr(NodeLoad) declaration. Placed after %include
// "fem_common.i" so it appears after swig::traits_from template is defined.
// ===================================================================
%{
namespace swig {
    template <>
    struct traits_from<NodeLoad> {
        static PyObject* from(const NodeLoad& val) {
            std::shared_ptr<NodeLoad>* smartresult = new std::shared_ptr<NodeLoad>(
                std::make_shared<NodeLoad>(val));
            // Use cached type lookup (SWIG_TypeQuery caches results internally)
            static swig_type_info* stype = nullptr;
            if (!stype) {
                stype = SWIG_TypeQuery("std::shared_ptr< NodeLoad > *");
            }
            return SWIG_NewPointerObj(
                SWIG_as_voidptr(smartresult),
                stype,
                SWIG_POINTER_OWN);
        }
    };
}
%}

// Python-specific extensions

// Material __str__ and __repr__ methods
%extend Material {
    %pythoncode %{
        def __str__(self):
            return self.to_string()

        def __repr__(self):
            return f"Material({self.to_string()})"
    %}
};

// ResponseSpectrumMethod Python helper
%extend ResponseSpectrumMethod {
    %pythoncode %{
        def set_spectrum(self, spectrum):
            """
            Set the response spectrum function.
            This method keeps a reference to prevent garbage collection.

            Args:
                spectrum: An IResponseSpectrum implementation
            """
            self._spectrum_ref = spectrum
            self.SpectrumFunction = spectrum
    %}
};

// FEDeformOperator Python properties
%extend FEDeformOperator {
    %pythoncode %{
        @property
        def operation_name(self):
            """Operation name"""
            if not hasattr(self, '_operation_name'):
                self._operation_name = ""
            return self._operation_name

        @operation_name.setter
        def operation_name(self, value):
            self._operation_name = value

        @property
        def operation_description(self):
            """Operation description (step summary, time, computation algorithm, etc.)"""
            if not hasattr(self, '_operation_description'):
                self._operation_description = ""
            return self._operation_description

        @operation_description.setter
        def operation_description(self, value):
            self._operation_description = value
    %}
};

// Add Python-specific module docstring
%pythoncode %{
"""
FEMNet - Finite Element Method Network Library

This module provides Python bindings for the FEMNet C++ library,
which implements various finite element analysis capabilities including:

- Static and dynamic analysis
- Buckling analysis
- Modal analysis (eigenvalue problems)
- Response spectrum analysis
- Various element types (truss, beam, plate, plane elements)

Example usage:
    import femnet

    # Create a model
    model = femnet.FEModel()

    # Add nodes
    model.AddNode(0, 0, 0)
    model.AddNode(1, 0, 0)

    # ... add elements, materials, loads, etc.
"""
%}

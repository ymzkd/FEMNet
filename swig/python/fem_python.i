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

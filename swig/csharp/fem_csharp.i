// fem_csharp.i - C#-specific SWIG interface file
%module(directors="1") FEMNet

// Include common definitions
%include "fem_common.i"

// C#-specific typemaps and extensions

// Material: Override ToString() method
%typemap(cscode) Material %{
    // Override ToString method
    public override string ToString()
    {
        return to_string();
    }
%}

// ResponseSpectrumMethod: Keep reference to prevent GC collection
%typemap(cscode) ResponseSpectrumMethod %{
    private IResponseSpectrum __spectrumRefs;

    /// <summary>
    /// Set the response spectrum function.
    /// This method keeps a reference to prevent garbage collection.
    /// </summary>
    /// <param name="spectrum">An IResponseSpectrum implementation</param>
    public void SetSpectrum(IResponseSpectrum spectrum) {
        __spectrumRefs = spectrum;
        this.SpectrumFunction = spectrum;
    }
%}

// FEDeformOperator: Add C# properties
%typemap(cscode) FEDeformOperator %{
    /// <summary>
    /// Operation Name
    /// </summary>
    public virtual string OperationName { get; set; } = "";

    /// <summary>
    /// Operation Description
    /// step summary, time, computation algorithm, etc.
    /// </summary>
    public virtual string OperationDescription { get; set; } = "";
%}

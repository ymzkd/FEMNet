cmake_minimum_required(VERSION 3.25)
project(FEMNet VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17) # C++17
set(MKL_INTERFACE lp64)
find_package(MKL CONFIG)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(Spectra CONFIG REQUIRED)

if (MKL_FOUND)
    message(STATUS "MKL found")
    add_definitions(-DUSE_MKL)
else()
    message(STATUS "MKL not found, proceeding without MKL")
endif()

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# === Library ===
# Set .Net project directory
if(NOT DEFINED NET_PROJECT_DIR)
    set(NET_PROJECT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../FEMNet/")
endif()

if(NOT DEFINED SWIG_OUTPUT_DIR)
    set(SWIG_OUTPUT_DIR "${NET_PROJECT_DIR}/")
endif()

# Add swig flags here
set(CMAKE_SWIG_FLAGS "")

# SWIGエントリ
set(SWIG_MODULE_fem_EXTRA_DEPS
    Elements/Elements.i
    Operator.i
    ModelOperator.i
)

set_property(SOURCE fem.i PROPERTY CPLUSPLUS ON)
set_source_files_properties(fem.i PROPERTIES
    SWIG_FLAGS "-includeall")

swig_add_library(fem
    TYPE SHARED
    LANGUAGE CSharp
    SOURCES 
        fem.i 
        Components.cpp 
        Model.cpp 
        Elements/ElementBase.cpp 
        Elements/BarElement.cpp 
        Elements/TrussElement.cpp 
        Elements/BeamElement.cpp 
        Elements/ComplexBeamElement.cpp 
        Elements/PlaneElement.cpp 
        Elements/TriPlaneElement.cpp 
        Elements/QuadPlaneElement.cpp 
        Elements/TriPlateElement.cpp 
        Elements/QuadPlateElement.cpp 
        LoadComponent.cpp 
        SeismicModule.cpp 
        FELinearStaticOp.cpp 
        FEDynamic.cpp 
        FEBucklingAnalysis.cpp 
        FEVibrateResult.cpp 
        ResponseSpectrumMethod.cpp
    OUTPUT_DIR ${SWIG_OUTPUT_DIR}
    OUTFILE_DIR ${SWIG_OUTPUT_DIR}
)

set_property(TARGET fem PROPERTY SWIG_COMPILE_OPTIONS -namespace ${PROJECT_NAME})

target_include_directories(fem PUBLIC ${CMAKE_SOURCE_DIR})
if (MKL_FOUND)
    target_compile_definitions(fem PUBLIC EIGEN_USE_MKL_ALL)
    target_compile_options(fem PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(fem PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(fem PUBLIC MKL::MKL Eigen3::Eigen Spectra::Spectra)
else()
    target_link_libraries(fem PUBLIC Eigen3::Eigen Spectra::Spectra)
endif()
# for copying fem.dll to .Net project dir
set_target_properties( fem
    PROPERTIES
    # These copy triangulate.dll on Windows to .Net project directory
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${NET_PROJECT_DIR}
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${NET_PROJECT_DIR}
    
    # This copies triangulate.so on Linux to .Net project directory
    LIBRARY_OUTPUT_DIRECTORY ${NET_PROJECT_DIR}
    # LIBRARY_OUTPUT_DIRECTORY ${NET_PROJECT_DIR}

    # Set address of C++ headers
    INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}
)

# === App ===
add_executable(sampleapp 
    app.cpp 
    Components.cpp 
    Model.cpp 
    Elements/ElementBase.cpp 
    Elements/BarElement.cpp 
    Elements/TrussElement.cpp 
    Elements/BeamElement.cpp 
    Elements/ComplexBeamElement.cpp 
    Elements/PlaneElement.cpp 
    Elements/TriPlaneElement.cpp 
    Elements/QuadPlaneElement.cpp 
    Elements/TriPlateElement.cpp 
    Elements/QuadPlateElement.cpp 
    LoadComponent.cpp 
    SeismicModule.cpp 
    FELinearStaticOp.cpp 
    FEDynamic.cpp 
    FEBucklingAnalysis.cpp 
    FEVibrateResult.cpp 
    ResponseSpectrumMethod.cpp
)

target_include_directories(sampleapp PUBLIC ${CMAKE_SOURCE_DIR})
if (MKL_FOUND)
    target_compile_definitions(sampleapp PUBLIC EIGEN_USE_MKL_ALL)
    target_compile_options(sampleapp PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(sampleapp PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(sampleapp PUBLIC MKL::MKL Eigen3::Eigen Spectra::Spectra)
else()
    target_link_libraries(sampleapp PUBLIC Eigen3::Eigen Spectra::Spectra)
endif()


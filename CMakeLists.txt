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

# === Build Options ===
option(BUILD_CSHARP "Build C# bindings" ON)
option(BUILD_PYTHON "Build Python bindings" OFF)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# Python settings
if(BUILD_PYTHON)
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
    message(STATUS "Python3 found: ${Python3_EXECUTABLE}")
    message(STATUS "Python3 include dirs: ${Python3_INCLUDE_DIRS}")
    message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
endif()

# === Common source files ===
set(FEMNET_SOURCES
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

# === C# Library ===
if(BUILD_CSHARP)
    # Set .Net project directory
    if(NOT DEFINED NET_PROJECT_DIR)
        set(NET_PROJECT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../FEMNet/")
    endif()

    if(NOT DEFINED SWIG_OUTPUT_DIR)
        set(SWIG_OUTPUT_DIR "${NET_PROJECT_DIR}/")
    endif()

    set_property(SOURCE swig/csharp/fem_csharp.i PROPERTY CPLUSPLUS ON)
    set_source_files_properties(swig/csharp/fem_csharp.i PROPERTIES
        SWIG_FLAGS "-includeall"
        INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/swig/common")

    swig_add_library(fem
        TYPE SHARED
        LANGUAGE CSharp
        SOURCES
            swig/csharp/fem_csharp.i
            ${FEMNET_SOURCES}
        OUTPUT_DIR ${SWIG_OUTPUT_DIR}
        OUTFILE_DIR ${SWIG_OUTPUT_DIR}
    )

    set_property(TARGET fem PROPERTY SWIG_COMPILE_OPTIONS -namespace ${PROJECT_NAME})
    target_include_directories(fem PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/swig/common)

    if (MKL_FOUND)
        target_compile_definitions(fem PUBLIC EIGEN_USE_MKL_ALL)
        target_compile_options(fem PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
        target_include_directories(fem PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(fem PUBLIC MKL::MKL Eigen3::Eigen Spectra::Spectra)
    else()
        target_link_libraries(fem PUBLIC Eigen3::Eigen Spectra::Spectra)
    endif()

    # for copying fem.dll to .Net project dir
    set_target_properties(fem
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY_RELEASE ${NET_PROJECT_DIR}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG ${NET_PROJECT_DIR}
        LIBRARY_OUTPUT_DIRECTORY ${NET_PROJECT_DIR}
        INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

# === Python Library ===
if(BUILD_PYTHON)
    if(NOT DEFINED PYTHON_OUTPUT_DIR)
        set(PYTHON_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/python/femnet/")
    endif()

    set_property(SOURCE swig/python/fem_python.i PROPERTY CPLUSPLUS ON)
    set_source_files_properties(swig/python/fem_python.i PROPERTIES
        SWIG_FLAGS "-py3"
        INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/swig/common")

    swig_add_library(femnet
        TYPE SHARED
        LANGUAGE python
        SOURCES
            swig/python/fem_python.i
            ${FEMNET_SOURCES}
        OUTPUT_DIR ${PYTHON_OUTPUT_DIR}
        OUTFILE_DIR ${PYTHON_OUTPUT_DIR}
    )

    target_include_directories(femnet PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/swig/common
        ${Python3_INCLUDE_DIRS})

    if (MKL_FOUND)
        target_compile_definitions(femnet PUBLIC EIGEN_USE_MKL_ALL)
        target_compile_options(femnet PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
        target_include_directories(femnet PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(femnet PUBLIC MKL::MKL Eigen3::Eigen Spectra::Spectra Python3::Python)
    else()
        target_link_libraries(femnet PUBLIC Eigen3::Eigen Spectra::Spectra Python3::Python)
    endif()

    set_target_properties(femnet
        PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PYTHON_OUTPUT_DIR}
        INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

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


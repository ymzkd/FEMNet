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

# === Optional CUDA support ===
option(USE_CUDA "Enable CUDA sparse solvers" OFF)
if(USE_CUDA)
    # MSVC 19.50+ (VS 2026+) is not yet supported by CUDA 12.x;
    # set --allow-unsupported-compiler BEFORE enable_language(CUDA)
    # so the compiler identification test can pass.
    if(MSVC AND MSVC_VERSION GREATER_EQUAL 1950)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --allow-unsupported-compiler")
    endif()

    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit REQUIRED)
        set(CUDA_FOUND TRUE)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_ARCHITECTURES 89)
        # Suppress Eigen warnings in nvcc and enable relaxed constexpr
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -diag-suppress=20012,20015 -Xcudafe --diag_suppress=20012 -Xcudafe --diag_suppress=20015 -w")
        message(STATUS "CUDA found: ${CMAKE_CUDA_COMPILER}")
        message(STATUS "CUDA Toolkit version: ${CUDAToolkit_VERSION}")
    else()
        set(CUDA_FOUND FALSE)
        message(STATUS "CUDA compiler not found, proceeding without GPU solvers")
    endif()
else()
    set(CUDA_FOUND FALSE)
    message(STATUS "CUDA disabled by user")
endif()

# CUDA solver sources
if(CUDA_FOUND)
    set(CUDA_SOLVER_SOURCES
        cuda/CudaSolverBase.cu
        cuda/CudaPCGSolver.cu
    )
endif()

# Helper function: configure CUDA for a target
function(configure_cuda_target target_name)
    if(CUDA_FOUND)
        target_compile_definitions(${target_name} PUBLIC USE_CUDA)
        target_link_libraries(${target_name} PUBLIC CUDA::cusolver CUDA::cusparse CUDA::cublas CUDA::cudart)
        set_target_properties(${target_name} PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
            CUDA_RESOLVE_DEVICE_SYMBOLS ON
        )
    endif()
endfunction()

# Helper function: configure MKL for a target
function(configure_mkl_target target_name)
    if(MKL_FOUND)
        target_compile_definitions(${target_name} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:EIGEN_USE_MKL_ALL>)
        target_compile_options(${target_name} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
        target_include_directories(${target_name} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_link_libraries(${target_name} PUBLIC MKL::MKL Eigen3::Eigen Spectra::Spectra)
    else()
        target_link_libraries(${target_name} PUBLIC Eigen3::Eigen Spectra::Spectra)
    endif()
endfunction()

# === Build Options ===
option(BUILD_CSHARP "Build C# bindings" ON)
option(BUILD_PYTHON "Build Python bindings" OFF)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

if(MSVC)
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:/utf-8>")
else()
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Python settings
if(BUILD_PYTHON)
    if(WIN32)
        # --- Windows: Search for Anaconda/Miniconda Python ---
        # Use system Python (Anaconda/user Python) instead of vcpkg's Python
        set(CONDA_SEARCH_PATHS
            "C:/ProgramData/anaconda3"
            "C:/ProgramData/miniconda3"
            "$ENV{USERPROFILE}/anaconda3"
            "$ENV{USERPROFILE}/miniconda3"
            "C:/tools/anaconda3"
            "C:/tools/miniconda3"
        )

        # Search for Anaconda Python with development files
        foreach(CONDA_PATH ${CONDA_SEARCH_PATHS})
            if(EXISTS "${CONDA_PATH}/python.exe" AND
               EXISTS "${CONDA_PATH}/include/Python.h" AND
               EXISTS "${CONDA_PATH}/libs")
                set(FOUND_PYTHON_ROOT "${CONDA_PATH}")
                break()
            endif()
        endforeach()

        if(FOUND_PYTHON_ROOT)
            message(STATUS "Found Anaconda/Miniconda at: ${FOUND_PYTHON_ROOT}")

            # Get Python version by running python --version
            execute_process(
                COMMAND "${FOUND_PYTHON_ROOT}/python.exe" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
                OUTPUT_VARIABLE PYTHON_VERSION_SHORT
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )

            execute_process(
                COMMAND "${FOUND_PYTHON_ROOT}/python.exe" -c "import sys; print(f'{sys.version_info.major}{sys.version_info.minor}')"
                OUTPUT_VARIABLE PYTHON_VERSION_NODOT
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )

            message(STATUS "Python version: ${PYTHON_VERSION_SHORT}")

            # Set Python variables to bypass vcpkg's wrapper
            set(Python3_EXECUTABLE "${FOUND_PYTHON_ROOT}/python.exe" CACHE FILEPATH "" FORCE)
            set(Python3_INCLUDE_DIRS "${FOUND_PYTHON_ROOT}/include" CACHE PATH "" FORCE)
            set(Python3_LIBRARIES "${FOUND_PYTHON_ROOT}/libs/python${PYTHON_VERSION_NODOT}.lib" CACHE FILEPATH "" FORCE)

            # Create imported target manually since vcpkg wrapper may override find_package
            add_library(Python3::Python SHARED IMPORTED)
            set_target_properties(Python3::Python PROPERTIES
                IMPORTED_LOCATION "${FOUND_PYTHON_ROOT}/python${PYTHON_VERSION_NODOT}.dll"
                IMPORTED_IMPLIB "${FOUND_PYTHON_ROOT}/libs/python${PYTHON_VERSION_NODOT}.lib"
                INTERFACE_INCLUDE_DIRECTORIES "${FOUND_PYTHON_ROOT}/include"
            )

            message(STATUS "Python3 found: ${Python3_EXECUTABLE}")
            message(STATUS "Python3 include dirs: ${Python3_INCLUDE_DIRS}")
            message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
        else()
            message(STATUS "No Anaconda/Miniconda found, falling back to find_package")
            find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
            message(STATUS "Python3 found: ${Python3_EXECUTABLE}")
            message(STATUS "Python3 include dirs: ${Python3_INCLUDE_DIRS}")
            message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
        endif()
    else()
        # --- Linux/macOS: Standard find_package ---
        find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
        message(STATUS "Python3 found: ${Python3_EXECUTABLE}")
        message(STATUS "Python3 include dirs: ${Python3_INCLUDE_DIRS}")
        message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
    endif()
endif()

# === Common source files ===
set(FEMNET_SOURCES
    Components.cpp
    Model.cpp
    SparseMatrixUtils.cpp
    RigidLink.cpp
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
    FEAnalysis.cpp
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
            $<$<BOOL:${CUDA_FOUND}>:${CUDA_SOLVER_SOURCES}>
        OUTPUT_DIR ${SWIG_OUTPUT_DIR}
        OUTFILE_DIR ${SWIG_OUTPUT_DIR}
    )

    set_property(TARGET fem PROPERTY SWIG_COMPILE_OPTIONS -namespace ${PROJECT_NAME})
    target_include_directories(fem PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/swig/common)

    configure_mkl_target(fem)
    configure_cuda_target(fem)

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
            $<$<BOOL:${CUDA_FOUND}>:${CUDA_SOLVER_SOURCES}>
        OUTPUT_DIR ${PYTHON_OUTPUT_DIR}
        OUTFILE_DIR ${PYTHON_OUTPUT_DIR}
    )

    target_include_directories(femnet PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/swig/common
        ${Python3_INCLUDE_DIRS})

    configure_mkl_target(femnet)
    target_link_libraries(femnet PUBLIC Python3::Python)
    configure_cuda_target(femnet)

    set_target_properties(femnet
        PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${PYTHON_OUTPUT_DIR}
        LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PYTHON_OUTPUT_DIR}
        LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PYTHON_OUTPUT_DIR}
        RUNTIME_OUTPUT_DIRECTORY ${PYTHON_OUTPUT_DIR}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PYTHON_OUTPUT_DIR}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PYTHON_OUTPUT_DIR}
        INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

# === App ===
add_executable(sampleapp
    app.cpp
    ${FEMNET_SOURCES}
    $<$<BOOL:${CUDA_FOUND}>:${CUDA_SOLVER_SOURCES}>
)

target_include_directories(sampleapp PUBLIC ${CMAKE_SOURCE_DIR})
configure_mkl_target(sampleapp)
configure_cuda_target(sampleapp)


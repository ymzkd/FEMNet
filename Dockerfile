FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install build tools and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    swig \
    python3-dev \
    python3-pip \
    libeigen3-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Spectra (header-only library, not available via apt)
RUN git clone --depth 1 --branch v1.0.1 https://github.com/yixuan/spectra.git /tmp/spectra \
    && cd /tmp/spectra \
    && cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local \
    && cmake --install build \
    && rm -rf /tmp/spectra

WORKDIR /femnet
COPY . .

# Step 1: Build and run C++ sample app
RUN cmake -B build -DBUILD_CSHARP=OFF -DBUILD_PYTHON=OFF \
    && cmake --build build \
    && ./build/sampleapp

# Step 2: Build Python bindings
RUN cmake -B build-python -DBUILD_CSHARP=OFF -DBUILD_PYTHON=ON \
    && cmake --build build-python

# Step 3: Install Python package
RUN pip install -e python/ --break-system-packages

# Step 4: Run all Python examples
RUN python3 python/examples/01_simple_truss.py \
    && python3 python/examples/02_cantilever_beam.py \
    && python3 python/examples/03_modal_analysis.py \
    && python3 python/examples/04_buckling_analysis.py \
    && python3 python/examples/05_dynamic_analysis.py \
    && python3 python/examples/06_plane_stress.py \
    && python3 python/examples/07_plate_element.py \
    && python3 python/examples/08_response_spectrum.py

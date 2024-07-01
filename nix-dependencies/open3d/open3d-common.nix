{
  lib,
  pkgs,
  buildPythonPackage,
  embree3
}:

buildPythonPackage {
  pname = "open3d-cpu";
  version = "0.18.0";
  pyproject = false;
  #src = fetchGit { url = "https://github.com/isl-org/Open3D.git"; rev = "0f06a149c4fb9406fd3e432a5cb0c024f38e2f0e"; };
  src = ./Open3D;
  build-system = with pkgs; [ cmake git openssl pkg-config ];
  dependencies = [
    pkgs.python3.pkgs.pybind11
    pkgs.python3.pkgs.python-lzf
  ];
  buildInputs = with pkgs; with xorg; [
    libX11
    libXrandr
    libXinerama
    vulkan-headers
    libXcursor
    libcxx
    curl
    assimp
    eigen
    fmt
    glew
    glfw
    imgui
    libjpeg
    jsoncpp
    msgpack-cxx
    nanoflann
    libpng
    qhull
    librealsense
    tinyobjloader
    vtk
    openssl_3_0
    boost
    directx-headers
    blas
    lapack
    tbb
    minizip
    fmt
    zeromq
    cppzmq
    embree3
  ];
  cmakeFlags = [ 
    "-DOPEN3D_USE_ONEAPI_PACKAGES=OFF"
    "-DUSE_BLAS=ON" 
    "-DBUNDLE_OPEN3D_ML=OFF"
    "-DBUILD_JUPYTER_EXTENSION=OFF"
    "-DBUILD_AZURE_KINECT=OFF"
    "-DBUILD_PYTORCH_OPS=OFF"
    "-DBUILD_CUDA_MODULE=OFF"
    "-DCMAKE_BUILD_TYPE=Release"
    "-DBUILD_UNIT_TESTS=OFF"
    "-DBUILD_BENCHMARKS=OFF"
    "-DBUILD_WEBRTC=OFF"
    "-DBUILD_SHARED_LIBS=OFF"
    "-DUSE_SYSTEM_BLAS=ON"
    "-DUSE_SYSTEM_ASSIMP=ON"
    "-DUSE_SYSTEM_CURL=ON"
    "-DUSE_SYSTEM_CUTLASS=ON"
    "-DUSE_SYSTEM_EIGEN3=ON"
    "-DUSE_SYSTEM_FILAMENT=OFF"
    "-DUSE_SYSTEM_FMT=ON"
    "-DUSE_SYSTEM_GLEW=ON"
    "-DUSE_SYSTEM_GLFW=ON"
    "-DUSE_SYSTEM_IMGUI=ON"
    "-DUSE_SYSTEM_JPEG=ON"
    "-DUSE_SYSTEM_JSONCPP=ON"
    "-DUSE_SYSTEM_LIBLZF=OFF"
    "-DUSE_SYSTEM_MSGPACK=ON"
    "-DUSE_SYSTEM_NANOFLANN=ON"
    "-DUSE_SYSTEM_OPENSSL=OFF"
    "-DUSE_SYSTEM_PNG=ON"
    "-DUSE_SYSTEM_PYBIND11=ON"
    "-DUSE_SYSTEM_QHULLCPP=ON"
    "-DUSE_SYSTEM_STDGPU=ON"
    "-DUSE_SYSTEM_TBB=ON"
    "-DUSE_SYSTEM_TINYOBJLOADER=ON"
    "-DUSE_SYSTEM_VTK=ON"
    "-DUSE_SYSTEM_ZEROMQ=ON"
    "-DWITH_MINIZIP=TRUE"
    "-DWITH_IPPICV=OFF"
    "-DBUILD_ISPC_MODULE=OFF"
    "-DUSE_SYSTEM_TBB=ON"
    "-DUSE_SYSTEM_FMT=ON"
    "-DUSE_SYSTEM_EMBREE=ON"
  ];
}


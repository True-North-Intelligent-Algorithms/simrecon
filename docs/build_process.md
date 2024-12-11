## C code Build Process  

The build uses CMake. The CMake file is [here](../src/CMakeLists.txt).  

The recommended way to start the build is from Visual Studio Code. For convenience, I have added the [.vscode](../.vscode) files to the GitHub repo. The directory locations need to be changed for local computers.

You will need to install the CMake and C/C++ extension. 

It is recommended to use the MinGW compiler.  

## Dependencies  

The dependencies required are listed in [settings.json](../.vscode/settings.json).  

Dependencies needed are [OpenBLAS](https://www.openblas.net/), [LAPACK](https://www.netlib.org/lapack/), and [FFTW](https://www.fftw.org/).  

IVE is also needed, but this is included in the source tree.  

```json
{
    "cmake.configureSettings": {
        "cmake.generator": "MinGW Makefiles",
        "CMAKE_BUILD_TYPE": "Debug",
        "CMAKE_INSTALL_PREFIX": "${env:CONDA_PREFIX}",
        "BLAS_LIBRARIES": "C:/Users/bnort/work/tools/OpenBLAS-0.3.27-x64/lib/libopenblas.lib",
        "LAPACK_LIBRARIES": "C:/Users/bnort/work/tools/LAPACKE/LAPACK.lib",
        "FFTW_INCLUDE_DIRS": "C:/Users/bnort/work/tools/fftw-3.3.5-dll64",
        "FFTW_FLOAT_LIB": "C:/Users/bnort/work/tools/fftw-3.3.5-dll64/libfftw3f-3.lib",
        "FFTW2_INCLUDE_DIRS": "D:/Janelia/scripts/SIMrecon_svn/usr_local/include",
        "FFTW2_FLOAT_LIB": "D:/Janelia/scripts/SIMrecon_svn/usr_local/lib/libsfftw.lib",
        "FFTW2_FLOAT_LIBR": "D:/Janelia/scripts/SIMrecon_svn/usr_local/lib/libsrfftw.lib"
    },
    "cmake.sourceDirectory": "${workspaceFolder}/src"
}
```


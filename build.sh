rm -rf cmake_build

mkdir cmake_build
cd cmake_build

#CXXFLAGS="$CXXFLAGS -Wfatal-errors -Wno-deprecated-declarations"
#echo $CMAKE_ARGS

cmake -G "Unix Makefiles" \
    -DFFTW_INCLUDE_DIRS="C:/Users/bnort/work/tools/fftw-3.3.5-dll64" \
    -DFFTW_FLOAT_LIB="C:/Users/bnort/work/tools/fftw-3.3.5-dll64/libfftw3f-3.lib" \
    -DBLAS_LIBRARIES="C:/Users/bnort/work/tools/OpenBLAS-0.3.27-x64/lib/libopenblas.lib" \
    -DLAPACK_LIBRARIES="C:/Users/bnort/work/tools/LAPACKE/LAPACK.lib" \
    ../src

#cmake -G"NMake Makefiles" \
#    ${CMAKE_ARGS} \
#    -DBUILD_MRC=ON \
#    -DBUILD_OTF_VIEWER=OFF \
#    -DCMAKE_BUILD_TYPE=Release \
#    -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
#    -DBLAS_LIBRARIES="C:/Users/bnort/work/tools/OpenBLAS-0.3.27-x64/lib/libopenblas.lib" \
#    -DLAPACK_LIBRARIES="C:/Users/bnort/work/tools/LAPACKE\LAPACK.lib" \
#    -DFFTW_INCLUDE_DIRS="C:/Users/bnort/work/tools/fftw-3.3.5-dll64" \
#    -DFFTW_FLOAT_LIB="C:/Users/bnort/work/tools/fftw-3.3.5-dll64/libfftw3f-3.lib" \
#    ../src

#nmake -j 2
#make install

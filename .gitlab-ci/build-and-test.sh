#!/bin/bash -x

[[ $OS_VERSION != "osx" ]] && ccache -s && ccache -M 5G

rm -rf build && mkdir -p build && (cd build && if [[ $BUILD_TYPE == "default" ]]; then cmake -G 'Unix Makefiles' .. -DCMAKE_BUILD_TYPE=Release -DNEKTAR_TEST_ALL=ON -DNEKTAR_BUILD_TIMINGS=ON -DNEKTAR_ERROR_ON_WARNINGS=OFF ..; elif [[ $BUILD_TYPE == "full" ]]; then cmake -DCMAKE_BUILD_TYPE:STRING=Debug -DNEKTAR_FULL_DEBUG:BOOL=ON -DNEKTAR_TEST_ALL:BOOL=ON -DNEKTAR_BUILD_TIMINGS:BOOL=ON -DNEKTAR_USE_ARPACK:BOOL=ON -DNEKTAR_USE_FFTW:BOOL=ON -DNEKTAR_USE_MPI:BOOL=ON -DNEKTAR_USE_SCOTCH:BOOL=ON -DNEKTAR_USE_PETSC:BOOL=ON -DNEKTAR_USE_HDF5:BOOL=ON -DNEKTAR_USE_MESHGEN:BOOL=ON -DNEKTAR_USE_CCM:BOOL=ON -DNEKTAR_BUILD_PYTHON:BOOL=ON -DNEKTAR_TEST_USE_HOSTFILE=ON -DNEKTAR_ERROR_ON_WARNINGS=OFF ..; fi) && \
    make -C build -j $NUM_CPUS all 2>&1 && \
    make -C build -j $NUM_CPUS install && \
    (cd build && ctest -j $NUM_CPUS --output-on-failure)

exit_code=$?
if [[ $exit_code -ne 0 ]]
then
    [[ $OS_VERSION != "osx" ]] && rm -rf build/dist
    exit $exit_code
fi

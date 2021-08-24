#!/bin/bash

set -e

# Create a build directory and compile against Nektar++
rm -rf build && mkdir build && cd build
cmake -DNektar++_DIR=$1 ..
make install

# Run test case in parallel.
cd ..
export OMPI_MCA_btl_vader_single_copy_mechanism=none
if [[ $EUID -eq 0 ]]; then
    mpiargs="-n $2 --allow-run-as-root"
else
    mpiargs="-n $2"
fi

test_output=`mpirun $mpiargs ./build/dist/ExampleSolver sample-laplace.xml | grep "L 2 error" | awk '{print ($7 < 1e-10)}'`
if [ "$test_output" -eq 1 ]; then
    echo "Test passed tolerance"
    exit 0
fi

# Test failed tolerance
echo "Test failed tolerance"
exit 1

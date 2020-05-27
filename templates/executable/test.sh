#!/bin/bash

set -e

# Create a build directory and compile against Nektar++
rm -rf build && mkdir build && cd build
cmake -DNektar++_DIR=$1 ..
make install

# Run test case in parallel.
cd ..
export OMPI_MCA_btl_vader_single_copy_mechanism=none
test_output=`mpirun -n $2 ./build/dist/ExampleSolver sample-laplace.xml | grep "L 2 error" | awk '{print ($7 < 1e-10)}'`
if [ "$test_output" -eq 1 ]; then
    echo "Test passed tolerance"
    exit 0
fi

# Test failed tolerance
echo "Test failed tolerance"
exit 1

#!/bin/bash

cd `dirname $0`

TARGET=nektar++-`cat VERSION`

if [ -d $TARGET ]; then
    rm -rf $TARGET
fi
mkdir -p $TARGET

echo "Generating file list..."
rsync -avqH --cvs-exclude --exclude-from dist-exclude * $TARGET

echo "Generating library distribution..."
tar -zc -f $TARGET.tar.gz $TARGET
rm -rf $TARGET

# Get documentation up to date
#cd docs/html/doxygen
#doxygen doxygen

#!/bin/bash
# @author Chris Cantwell
#
# This script generates Nektar++ distributions (nektar++-VERSION.tar.gz)
# The file VERSION contains the version number of the generated release.
#
# @requires rsync doxygen tar

BASE=`dirname $0`
cd $BASE

TARGET=nektar++-`cat VERSION`

# Make Code tree target
if [ -d $TARGET ]; then
    rm -rf $TARGET
fi
mkdir -p $TARGET

# Create code tree
echo "Generating code tree..."
rsync -avqH --cvs-exclude --exclude-from dist-exclude * $TARGET

# Package code tree
echo "Packaging code distribution..."
tar -zc -f $TARGET.tar.gz $TARGET

# Clean up
rm -rf $TARGET

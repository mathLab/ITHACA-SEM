#!/bin/bash
# @author Chris Cantwell
#
# This script generates Nektar++ distributions, comprising of the following:
# - A code distribution (nektar++-VERSION.tar.gz)
# - A website package (nektar++-web-VERSION.tar.gz)
# The file VERSION contains the version number of the generated release.
# A working copy of the ThirdParty repository should be placed in a subdirectory
# called ThirdParty, or sym-linked to such from elsewhere. If available, the
# ThirdParty distribution will also be compiled and included in the web package.
#
# @requires rsync doxygen tar

BASE=`dirname $0`
cd $BASE

TARGET=nektar++-`cat VERSION`

# Make Web tree target
if [ -d $TARGET-web ]; then
    rm -rf $TARGET-web
fi
mkdir -p $TARGET-web
mkdir -p $TARGET-web/downloads/

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

# Generate ThirdParty package if available
if [ -d ThirdParty -o -h ThirdParty ]; then
    if [ -f ThirdParty/distribute.sh ]; then
        ThirdParty/distribute.sh
        mv ThirdParty/ThirdParty-*.tar.gz $TARGET-web/downloads/
    else
        echo "ThirdParty directory exists, but without distribution script."
    fi
else
    echo "ThirdParty not available. Please package separately."
fi

# Generate documentation for distributed code
echo "Generating doxygen docs...this will take a while..."
cd $TARGET/docs/html/doxygen
doxygen doxygen > /dev/null 2>&1

cd ../../../../

# Create web tree
echo "Generating web tree..."
rsync -avqH --cvs-exclude --exclude='code' docs/html/* $TARGET-web
mv $TARGET/docs/html/code $TARGET-web/
mv $TARGET.tar.gz $TARGET-web/downloads/

# Package web tree
echo "Packaging web distribution..."
tar -zc -f $TARGET-web.tar.gz $TARGET-web
mv $TARGET-web/downloads/$TARGET.tar.gz .

# Clean up
rm -rf $TARGET-web
rm -rf $TARGET

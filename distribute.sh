#!/bin/bash

BASE=`dirname $0`
cd $BASE

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
cd $BASE
echo "Generating doxygen docs..."
cd docs/html/doxygen
doxygen doxygen

cd ../../../
if [ -d $TARGET-web ]; then
    rm -rf $TARGET-web
fi
mkdir -p $TARGET-web

echo "Generating website tree..."
rsync -avqH --cvs-exclude docs/html/* $TARGET-web
cp $TARGET.tar.gz $TARGET-web/downloads/

echo "Generating website tarball..."
tar -zc -f $TARGET-web.tar.gz $TARGET-web
rm -rf $TARGET-web

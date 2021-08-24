#/bin/bash
echo Patching all library links in OCC lib files and nektar libs/exec

path=$1
install_path=$2

for file in $path/*.dylib
do
    if [ -L $file ]; then
        continue
    fi
    echo Repairing: $file
    install_name_tool -id $install_path/`basename $file` $file
    DYLIBS=`otool -L $file | grep -v "/" | awk -F' ' '{ print $1 }'`
    for dylib in $DYLIBS
    do
        install_name_tool -change $dylib $install_path/`basename $dylib` $file
    done
done

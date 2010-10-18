# Generate list of available solvers
FILE(GLOB dir_list "*")
LIST(REMOVE_ITEM dir_list "${CMAKE_SOURCE_DIR}/PreProcessing/.svn")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})

SET(USERGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR})
SET(USERGUIDE ${CMAKE_BINARY_DIR}/docs/user-guide)

FILE(MAKE_DIRECTORY ${USERGUIDE}/html)

FIND_PROGRAM(PDFLATEX pdflatex)
MARK_AS_ADVANCED(PDFLATEX)
FIND_PROGRAM(BIBTEX bibtex)
MARK_AS_ADVANCED(BIBTEX)
FIND_PROGRAM(MAKEINDEX makeindex)
MARK_AS_ADVANCED(MAKEINDEX)
FIND_PROGRAM(HTLATEX htlatex)
MARK_AS_ADVANCED(HTLATEX)

ADD_CUSTOM_TARGET(user-guide-html
    export TEXINPUTS=${USERGUIDESRC}//: &&
    ${HTLATEX} ${USERGUIDESRC}/user-guide.tex
    "${USERGUIDESRC}/styling.cfg,html,3,next,NoFonts"
    COMMAND BIBINPUTS=${CMAKE_SOURCE_DIR}/docs ${BIBTEX} user-guide.aux
    COMMAND TEXINPUTS=${USERGUIDESRC}//:
    ${HTLATEX} ${USERGUIDESRC}/user-guide.tex 
    "${USERGUIDESRC}/styling.cfg,html,3,next,NoFonts"
    WORKING_DIRECTORY ${USERGUIDE}/html
)

# If tex4ht successful, create img dir and copy images across
FILE(GLOB_RECURSE imgfiles "img/*.png" "img/*.jpg" "*/img/*.png" "*/img/*.jpg")
ADD_CUSTOM_COMMAND(TARGET user-guide-html
    POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${USERGUIDE}/html/img)
FOREACH(img ${imgfiles})
    ADD_CUSTOM_COMMAND(TARGET user-guide-html
        POST_BUILD COMMAND 
            ${CMAKE_COMMAND} -E copy ${img} ${USERGUIDE}/html/img)
ENDFOREACH()

FILE(GLOB_RECURSE pdffiles "*/img/*.pdf")
FIND_PROGRAM(CONVERT convert)
FOREACH(pdf ${pdffiles})
    GET_FILENAME_COMPONENT(BASENAME ${pdf} NAME_WE)
    ADD_CUSTOM_COMMAND(TARGET user-guide-html
        POST_BUILD COMMAND
        ${CONVERT} ${pdf} ${USERGUIDE}/html/img/${BASENAME}.png)
ENDFOREACH()

ADD_CUSTOM_TARGET(user-guide-pdf
    ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
    COMMAND BIBINPUTS=${CMAKE_SOURCE_DIR}/docs TEXMFOUTPUT=${USERGUIDE}
    ${BIBTEX} ${USERGUIDE}/user-guide.aux
    COMMAND TEXMFOUTPUT=${USERGUIDE} ${MAKEINDEX} ${USERGUIDE}/user-guide.idx
    COMMAND ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
    COMMAND ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
    WORKING_DIRECTORY ${USERGUIDESRC}
)

# Install output, if it is generated.
INSTALL(CODE "IF (EXISTS ${USERGUIDE}/html)
    FILE(GLOB UGHTML ${USERGUIDE}/html/*)
    LIST(LENGTH UGHTML UGHTML_LEN)
    IF(UGHTML_LEN GREATER 0)
        FILE(INSTALL ${USERGUIDE}/html
             DESTINATION \${CMAKE_INSTALL_PREFIX}/${NEKTAR_DOC_DIR}/user-guide)
    ENDIF()
ENDIF()
IF (EXISTS ${USERGUIDE}/user-guide.pdf)
    FILE(INSTALL ${USERGUIDE}/user-guide.pdf
         DESTINATION \${CMAKE_INSTALL_PREFIX}/${NEKTAR_DOC_DIR}/user-guide)
ENDIF()")

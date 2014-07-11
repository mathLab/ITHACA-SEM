SET(USERGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR})
SET(USERGUIDE ${CMAKE_BINARY_DIR}/docs/user-guide)

FILE(MAKE_DIRECTORY ${USERGUIDE}/html)

FIND_PROGRAM(HTLATEX htlatex)
ADD_CUSTOM_TARGET(user-guide-html
    export TEXINPUTS=${USERGUIDESRC}//: && ${HTLATEX} ${USERGUIDESRC}/user-guide.tex 
        "${USERGUIDESRC}/styling.cfg,html,3,next"
    WORKING_DIRECTORY ${USERGUIDE}/html
)

FIND_PROGRAM(PDFLATEX pdflatex)
ADD_CUSTOM_TARGET(user-guide-pdf
    ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
    COMMAND ${PDFLATEX} --output-directory ${USERGUIDE}
        ${USERGUIDESRC}/user-guide.tex
    WORKING_DIRECTORY ${USERGUIDESRC}
)

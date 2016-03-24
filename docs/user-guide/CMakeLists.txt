SET(USERGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR})
SET(USERGUIDE ${CMAKE_BINARY_DIR}/docs/user-guide)

FILE(MAKE_DIRECTORY ${USERGUIDE}/html)

FIND_PROGRAM(HTLATEX htlatex)
MARK_AS_ADVANCED(HTLATEX)
ADD_CUSTOM_TARGET(user-guide-html
    export TEXINPUTS=${CMAKE_SOURCE_DIR}//:${USERGUIDESRC}//: &&
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


FIND_PROGRAM(PDFLATEX pdflatex)
MARK_AS_ADVANCED(PDFLATEX)
FIND_PROGRAM(BIBTEX bibtex)
MARK_AS_ADVANCED(BIBTEX)
FIND_PROGRAM(MAKEINDEX makeindex)
MARK_AS_ADVANCED(MAKEINDEX)

ADD_CUSTOM_TARGET(user-guide-pdf
    export TEXINPUTS=${CMAKE_SOURCE_DIR}//: && 
    ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
    COMMAND TEXMFOUTPUT=${USERGUIDE} ${BIBTEX} ${USERGUIDE}/user-guide.aux
    COMMAND TEXMFOUTPUT=${USERGUIDE} ${MAKEINDEX} ${USERGUIDE}/user-guide.idx
    COMMAND export TEXINPUTS=${CMAKE_SOURCE_DIR}//: && 
    ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
    WORKING_DIRECTORY ${USERGUIDESRC}
)

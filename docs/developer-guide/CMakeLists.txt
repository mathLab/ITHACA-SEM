SET(DEVGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR})
SET(DEVGUIDE ${CMAKE_BINARY_DIR}/docs/developer-guide)

FILE(MAKE_DIRECTORY ${DEVGUIDE}/html)

FIND_PROGRAM(HTLATEX htlatex)
ADD_CUSTOM_TARGET(developer-guide-html
    export TEXINPUTS=${CMAKE_SOURCE_DIR}//:${DEVGUIDESRC}//: && 
    ${HTLATEX} ${DEVGUIDESRC}/developer-guide.tex 
    "${DEVGUIDESRC}/styling.cfg,html,3,next,NoFonts"
    WORKING_DIRECTORY ${DEVGUIDE}/html
)

# If tex4ht successful, create img dir and copy images across
FILE(GLOB_RECURSE imgfiles "*/img/*.png" "*/img/*.jpg")
ADD_CUSTOM_COMMAND(TARGET developer-guide-html
    POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${DEVGUIDE}/html/img)
FOREACH(img ${imgfiles})
    ADD_CUSTOM_COMMAND(TARGET developer-guide-html
        POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${img} ${DEVGUIDE}/html/img)
ENDFOREACH()

FILE(GLOB_RECURSE pdffiles "*/img/*.pdf")
FIND_PROGRAM(CONVERT convert)
FOREACH(pdf ${pdffiles})
    GET_FILENAME_COMPONENT(BASENAME ${pdf} NAME_WE)
    ADD_CUSTOM_COMMAND(TARGET developer-guide-html
        POST_BUILD COMMAND
        ${CONVERT} ${pdf} ${DEVGUIDE}/html/img/${BASENAME}.png)
ENDFOREACH()

FIND_PROGRAM(PDFLATEX pdflatex)
FIND_PROGRAM(BIBTEX bibtex)
FIND_PROGRAM(MAKEINDEX makeindex)
ADD_CUSTOM_TARGET(developer-guide-pdf
    export TEXINPUTS=${CMAKE_SOURCE_DIR}//: &&
    ${PDFLATEX} --output-directory ${DEVGUIDE} ${DEVGUIDESRC}/developer-guide.tex
    COMMAND TEXMFOUTPUT=${DEVGUIDE} ${BIBTEX} ${DEVGUIDE}/developer-guide.aux
    COMMAND TEXMFOUTPUT=${DEVGUIDE} ${MAKEINDEX} ${DEVGUIDE}/developer-guide.idx
    COMMAND export TEXINPUTS=${CMAKE_SOURCE_DIR}//: &&
    ${PDFLATEX} --output-directory ${DEVGUIDE} ${DEVGUIDESRC}/developer-guide.tex
    WORKING_DIRECTORY ${DEVGUIDESRC}
)

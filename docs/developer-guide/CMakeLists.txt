SET(DEVGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR})
SET(DEVGUIDE ${CMAKE_BINARY_DIR}/docs/developer-guide)

FILE(MAKE_DIRECTORY ${DEVGUIDE}/html)

FIND_PROGRAM(PDFLATEX pdflatex)
MARK_AS_ADVANCED(PDFLATEX)
FIND_PROGRAM(BIBTEX bibtex)
MARK_AS_ADVANCED(BIBTEX)
FIND_PROGRAM(MAKEINDEX makeindex)
MARK_AS_ADVANCED(MAKEINDEX)
FIND_PROGRAM(HTLATEX htlatex)
MARK_AS_ADVANCED(HTLATEX)

ADD_CUSTOM_TARGET(developer-guide-html
    TEXINPUTS=${CMAKE_SOURCE_DIR}:${DEVGUIDESRC}//:
    ${HTLATEX} ${DEVGUIDESRC}/developers-guide.tex 
    "${DEVGUIDESRC}/styling.cfg,html,3,next,NoFonts"
    COMMAND BIBINPUTS=${DEVGUIDESRC} ${BIBTEX} developers-guide.aux
    COMMAND TEXINPUTS=${CMAKE_SOURCE_DIR}:${DEVGUIDESRC}//:
    ${HTLATEX} ${DEVGUIDESRC}/developers-guide.tex 
    "${DEVGUIDESRC}/styling.cfg,html,3,next,NoFonts"
    WORKING_DIRECTORY ${DEVGUIDE}/html
)

# If tex4ht successful, create img dir and copy images across
FILE(GLOB_RECURSE imgfiles "${DEVGUIDESRC}/*.png" "${DEVGUIDESRC}/*.jpg")
ADD_CUSTOM_COMMAND(TARGET developer-guide-html
    POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${DEVGUIDE}/html/img)
FOREACH(img ${imgfiles})
    ADD_CUSTOM_COMMAND(TARGET developer-guide-html
        POST_BUILD COMMAND 
            ${CMAKE_COMMAND} -E copy ${img} ${DEVGUIDE}/html/img)
ENDFOREACH()

FILE(GLOB_RECURSE pdffiles "${DEVGUIDESRC}/*.pdf")
FIND_PROGRAM(CONVERT convert)
FOREACH(pdf ${pdffiles})
    GET_FILENAME_COMPONENT(BASENAME ${pdf} NAME_WE)
    ADD_CUSTOM_COMMAND(TARGET developer-guide-html
        POST_BUILD COMMAND
        ${CONVERT} ${pdf} ${DEVGUIDE}/html/img/${BASENAME}.png)
ENDFOREACH()

ADD_CUSTOM_TARGET(developer-guide-pdf
    export TEXINPUTS=${CMAKE_SOURCE_DIR}//: && 
    ${PDFLATEX} --output-directory ${DEVGUIDE} ${DEVGUIDESRC}/developers-guide.tex
    COMMAND TEXMFOUTPUT=${DEVGUIDE} ${BIBTEX} ${DEVGUIDE}/developers-guide.aux
    COMMAND TEXMFOUTPUT=${DEVGUIDE} ${MAKEINDEX} ${DEVGUIDE}/developers-guide.idx
    COMMAND TEXINPUTS=${CMAKE_SOURCE_DIR}//:  
    ${PDFLATEX} --output-directory ${DEVGUIDE} ${DEVGUIDESRC}/developers-guide.tex
    COMMAND TEXINPUTS=${CMAKE_SOURCE_DIR}//:  
    ${PDFLATEX} --output-directory ${DEVGUIDE} ${DEVGUIDESRC}/developers-guide.tex
    WORKING_DIRECTORY ${DEVGUIDESRC}
)

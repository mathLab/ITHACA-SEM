SET(DEVGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR})
SET(DEVGUIDE ${CMAKE_BINARY_DIR}/docs/developer-guide)

FILE(MAKE_DIRECTORY ${DEVGUIDE}/html)

FIND_PROGRAM(HTLATEX htlatex)
ADD_CUSTOM_TARGET(developer-guide-html
    export TEXINPUTS=${DEVGUIDESRC}//: && ${HTLATEX} ${DEVGUIDESRC}/developer-guide.tex 
    "${DEVGUIDESRC}/styling.cfg,html,2,next,NoFonts"
    WORKING_DIRECTORY ${DEVGUIDE}/html
)

FIND_PROGRAM(PDFLATEX pdflatex)
FIND_PROGRAM(BIBTEX bibtex)
FIND_PROGRAM(MAKEINDEX makeindex)
ADD_CUSTOM_TARGET(developer-guide-pdf
    ${PDFLATEX} --output-directory ${DEVGUIDE} ${DEVGUIDESRC}/developer-guide.tex
    COMMAND TEXMFOUTPUT=${DEVGUIDE} ${BIBTEX} ${DEVGUIDE}/developer-guide.aux
    COMMAND TEXMFOUTPUT=${DEVGUIDE} ${MAKEINDEX} ${DEVGUIDE}/developer-guide.idx
    COMMAND ${PDFLATEX} --output-directory ${DEVGUIDE} 
            ${DEVGUIDESRC}/developer-guide.tex
    WORKING_DIRECTORY ${DEVGUIDESRC}
)

message(STATUS "downloading...
     src='http://www.nektar.info/thirdparty/modmetis-5.1.0_1.tar.bz2'
     dst='/tmp/ybao/nektar++/ThirdParty/modmetis-5.1.0_1.tar.bz2'
     timeout='none'")

file(DOWNLOAD
  "http://www.nektar.info/thirdparty/modmetis-5.1.0_1.tar.bz2"
  "/tmp/ybao/nektar++/ThirdParty/modmetis-5.1.0_1.tar.bz2"
  SHOW_PROGRESS
  EXPECTED_MD5;6c6816aea0f53db6c71b1d98ed4ad42b
  # no TIMEOUT
  STATUS status
  LOG log)

list(GET status 0 status_code)
list(GET status 1 status_string)

if(NOT status_code EQUAL 0)
  message(FATAL_ERROR "error: downloading 'http://www.nektar.info/thirdparty/modmetis-5.1.0_1.tar.bz2' failed
  status_code: ${status_code}
  status_string: ${status_string}
  log: ${log}
")
endif()

message(STATUS "downloading... done")

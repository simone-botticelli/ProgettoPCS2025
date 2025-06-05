list(APPEND src_txt_headers ${CMAKE_CURRENT_SOURCE_DIR}/ExportTxt.hpp)
list(APPEND src_txt_sources ${CMAKE_CURRENT_SOURCE_DIR}/ExportTxt.cpp)

list(APPEND src_txt_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_txt_sources ${src_txt_sources} PARENT_SCOPE)
set(src_txt_headers ${src_txt_headers} PARENT_SCOPE)
set(src_txt_includes ${src_txt_includes} PARENT_SCOPE)
